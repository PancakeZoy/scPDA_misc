# The search folder is too large to upload, please contact us if you need it.
require(dplyr)
require(Seurat)
require(h5read)
require(rhdf5)
require(mclust)
require(pROC)
require(glue)
require(purrr)
require(ggplot2)
require(gridExtra)
require(grid)
require(rlist)
"%ni%" <- Negate("%in%")

roc_auc = function(prob, prot, pos, types){
  label = if_else(types %in% pos, 1, 0)
  roc <- roc(label, prob[prot,], quiet=TRUE); auc = sprintf('%.2f', auc(roc))
  return(list(roc=roc, auc=auc))
}
################################# parameters #################################
layer1_layer2_hidden <- list(
  c(32, 16, 5), c(32, 16, 10),
  c(64, 32, 5), c(64, 32, 10), c(64, 32, 15),
  c(100, 50, 10), c(100, 50, 15), c(100, 50, 20),
  c(150, 75, 15), c(150, 75, 20),
  c(200, 100, 15), c(200, 100, 20)
)
batch_size <- c(256)
n_epochs <- c(500)
lr <- c(0.0005, 0.001, 0.005, 0.01)
gamma <- c(0.99)
kld_weight <- c(0.01, 0.05, 0.1, 0.25, 0.4)
alpha <- c(0.01, 0.03, 0.05, 0.1, 0.3, 0.5)

param_grid <- expand.grid(
  layer = seq_along(layer1_layer2_hidden),
  batch_size = batch_size,
  n_epochs = n_epochs,
  lr = lr,
  gamma = gamma,
  kld_weight = kld_weight,
  alpha = alpha,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
param_grid$layer <- purrr::map(param_grid$layer, ~ layer1_layer2_hidden[[.x]])
param_grid$layer <- purrr::map_chr(param_grid$layer, ~ paste0("(", paste(.x, collapse = ", "), ")"))
param_grid = param_grid %>% select(layer, lr, kld_weight, alpha)
param_grid$orig_index <- seq_len(nrow(param_grid))
param_grid = param_grid %>% filter(lr!=0.001) %>%
  filter(alpha %ni% c(0.01, 0.5)) %>%
  filter(kld_weight != 0.05) %>%
  filter(layer %in% c("(100, 50, 15)", "(100, 50, 20)", "(150, 75, 15)"))

################################# AUC #################################
tea <- readRDS('data/teaseq.rds')
empty = tea$empty; tea=tea$teaseq; tea=tea%>%subset(celltype != 'Unknown')
raw.counts = tea %>% GetAssayData(slot='counts') %>% as.matrix()
cell_log = raw.counts %>% log1p
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)
type_vec=tea$celltype

res = c()
marker.list <- list(
  list(prot='CD14', pos=c('CD14 Mono', 'Dendritic'), pos.name='DC/CD14 Mono'),
  list(prot='CD11c', pos=c('CD14 Mono', 'CD16 Mono', 'NK','Dendritic'), pos.name='DC/Mono/NK'),
  list(prot='CD123', pos=c('CD14 Mono', 'CD16 Mono', 'B','Dendritic'), pos.name='DC/Mono/B'),
  list(prot='CD141', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'),
  list(prot='CD16', pos=c('CD16 Mono', 'Dendritic', 'NK'), pos.name='DC/CD16 Mono/NK'),
  list(prot='CD172a', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'),
  list(prot='CD185', pos=c('B'), pos.name='B'),
  list(prot='CD19', pos=c('B'), pos.name='B'),
  list(prot='CD192', pos=c('CD14 Mono'), pos.name='CD14 Mono'),
  list(prot='CD21', pos=c('B'), pos.name='B'),
  list(prot='CD24', pos=c('B'), pos.name='B'),
  list(prot='CD27', pos=c('CD4 Memory', 'CD4 Naive', 'CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='T Cells'),
  list(prot='CD3', pos=c('CD4 Memory', 'CD4 Naive', 'CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='T Cells'),
  list(prot='CD40', pos=c('B', 'Dendritic', 'CD16 Mono', 'CD14 Mono'), pos.name='DC/B/Mono'),
  list(prot='CD56', pos=c('MAIT', 'NK'), pos.name='NK/MAIT'),
  list(prot='CD86', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'),
  list(prot='CD8a', pos=c('CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='CD8 T/MAIT'),
  list(prot='FceRI', pos=c('Dendritic'), pos.name='Dendritic'),
  list(prot='IgD', pos=c('B'), pos.name='B'),
  list(prot='IgM', pos=c('B'), pos.name='B')
)

for (i in (param_grid$orig_index-1)){
  path <- glue("/Volumes/SSD_1TB/search_scPDA/TEA/{i}.h5")
  pi = h5read(path, 'pi') %>% as.matrix()
  colnames(pi)=colnames(raw.counts); rownames(pi)=rownames(raw.counts)
  auc_scpda_ls = c();  auc_names=c()
  for (marker in marker.list){
    prot = marker$prot; pos=marker$pos
    auc_names = append(auc_names, prot)
    
    auc_scpda = roc_auc(pi, prot=prot, pos=pos, types=type_vec)$auc
    auc_scpda_ls = append(auc_scpda_ls, as.numeric(auc_scpda))
  }
  avg = mean(auc_scpda_ls)
  res = c(res, avg)
}
param_grid$auc = res
############################# NRC #####################################
sobj = readRDS('data/dsb.rds')
meta = sobj@meta.data
raw.counts = GetAssayData(sobj, assay='Protein', slot='counts') %>% as.matrix() %>% t %>% as.data.frame()
colnames(raw.counts)[grepl("sotype", colnames(raw.counts))] = paste0("Isotype", 1:4)

res = list(avg_unstain=c(), median_unstain=c(), avg_nonexpr=c(), median_nonexpr=c())
for (i in (param_grid$orig_index-1)){
  path = glue('/Volumes/SSD_1TB/search_scPDA/DSB/{i}.h5')
  pi = h5read(path, 'pi') %>% as.matrix() %>% t()
  colnames(pi)=colnames(raw.counts); rownames(pi)=rownames(raw.counts)
  scPDA.counts=round(raw.counts * (1-pi))
  # Neg Control Cells
  neg_idx = meta$Type1=='Neg Control'
  before = raw.counts[neg_idx,]
  after = scPDA.counts[neg_idx,]
  nrc.unstain = 1-colMeans(after)/colMeans(before)
  res$avg_unstain = c(res$avg_unstain, mean(nrc.unstain))
  res$median_unstain = c(res$median_unstain, median(nrc.unstain))
  # Non-expressive Proteins
  nonexpr_markers = c('CD117','CD137','CD138','CD206','CD223','CD273','CD80', 
                      "Isotype1","Isotype2","Isotype3","Isotype4")
  before = raw.counts %>% select(all_of(nonexpr_markers))
  after = scPDA.counts %>% select(all_of(nonexpr_markers))
  nrc.nonexpr = 1-colMeans(after)/colMeans(before)
  res$avg_nonexpr = c(res$avg_nonexpr, mean(nrc.nonexpr))
  res$median_nonexpr = c(res$median_nonexpr, median(nrc.nonexpr))
}
res = list.cbind(res)
param_grid = cbind(param_grid, res)


################################ Plot ################################
plist1 = list(); plist2 = list(); plist3 = list()
cols <- c("layer", "lr", "kld_weight", "alpha")
names(cols) = c("(1st layer, 2nd layer, bottle neck)", "learning rate", "KLD weight", "alpha")

layer_levels <- c("(100, 50, 15)", "(100, 50, 20)", "(150, 75, 15)")
param_grid$layer <- factor(param_grid$layer, levels = layer_levels)
param_grid$lr <- factor(param_grid$lr, levels = sort(unique(param_grid$lr)))
param_grid$kld_weight <- factor(param_grid$kld_weight, levels = sort(unique(param_grid$kld_weight)))
param_grid$alpha <- factor(param_grid$alpha, levels = sort(unique(param_grid$alpha)))

for (i in 1:length(cols)) {
  col = cols[i]
  
  p1 <- ggplot(param_grid, aes(x = .data[[col]], y = auc, fill = .data[[col]])) +
    geom_boxplot(outlier.shape = NA) + 
    labs(x = names(cols)[i], y = "AUC") +
    theme_minimal() + NoLegend() + ylim(0.8, 1) + 
    geom_jitter(width = 0.2, alpha = 1, size = 0.1)
  plist1 = c(plist1, list(p1))
  
  p2 <- ggplot(param_grid, aes(x = .data[[col]], y = avg_unstain, fill = .data[[col]])) +
    geom_boxplot(outlier.shape = NA) + 
    labs(x = names(cols)[i], y = "NRC for unstained cells") +
    theme_minimal() + NoLegend() + ylim(0.8, 1) + 
    geom_jitter(width = 0.2, alpha = 1, size = 0.1)
  plist2 = c(plist2, list(p2))
  
  p3 <- ggplot(param_grid, aes(x = .data[[col]], y = avg_nonexpr, fill = .data[[col]])) +
    geom_boxplot(outlier.shape = NA) + 
    labs(x = names(cols)[i], y = "NRC for nonexpressive proteins") +
    theme_minimal() + NoLegend() + ylim(0.8, 1) + 
    geom_jitter(width = 0.2, alpha = 1, size = 0.1)
  plist3 = c(plist3, list(p3))
}
plot1 = grid.arrange(grobs = plist1, nrow = 1, ncol = 4, 
                     top = textGrob("Sensitivity of AUC to Hyperparameters", 
                                    gp = gpar(fontface = "bold", fontsize=20)))
plot2 = grid.arrange(grobs = plist2, nrow = 1, ncol = 4, 
                     top = textGrob("Sensitivity of NRC for Unstained Cells to Hyperparameters", 
                                    gp = gpar(fontface = "bold", fontsize=20)))
plot3 = grid.arrange(grobs = plist3, nrow = 1, ncol = 4, 
                     top = textGrob("Sensitivity of NRC for nonexpressive proteins to Hyperparameters", 
                                    gp = gpar(fontface = "bold", fontsize=20)))


full_plot <- arrangeGrob(
  grobs = list(
    arrangeGrob(plot2, top = textGrob("a", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                      gp = gpar(fontsize = 30, fontface = "bold"))),
    grid::nullGrob(),
    arrangeGrob(plot3, top = textGrob("b", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                      gp = gpar(fontsize = 30, fontface = "bold"))),
    grid::nullGrob(),
    arrangeGrob(plot1, top = textGrob("c", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"),
                                      gp = gpar(fontsize = 30, fontface = "bold")))
  ),
  nrow = 5,
  heights = c(1, 0.1, 1, 0.1, 1)
)
grid.newpage()
grid.draw(full_plot)
ggsave('fig_reprod/supp/Tuning.png', plot=full_plot, width = 19, height = 10, dpi = 200)
