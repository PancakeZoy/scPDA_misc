library(here)
setwd(here())
source('code/tools.R')
theme_PlainTitle = theme(plot.title = element_text(face = "plain", size=10))
title.center = theme(plot.title = element_text(hjust = 0.5, face='plain'))
box = annotation_custom(grob = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 2)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

ridge = function(df, levels=NA, min_h, bold_labels){
  df_long <- df %>% gather(key="Variable", value="Value")
  
  if (!any(is.na(levels))){
    df_long$Variable <- factor(df_long$Variable, levels = levels)
  }
  
  # Create formatted labels
  labels <- levels(df_long$Variable)
  labels_formatted <- sapply(labels, function(l) {
    if (!is.null(bold_labels) && l %in% bold_labels) {
      styled_label <- sprintf("**<span style='color:red; font-style:italic; text-decoration:underline;'>%s</span>**", l)
    } else {
      l
    }
  })
  
  ggplot(df_long, aes(x = Value, y = Variable, fill = Variable)) +
    geom_density_ridges(scale = 3, rel_min_height = min_h) + 
    geom_vline(xintercept=0, linetype='dashed') +
    theme_minimal() +
    theme(
      legend.position = "none",  
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = ggtext::element_markdown()# Enable Markdown rendering
    ) +
    scale_y_discrete(labels = labels_formatted)
}
#################### Raw Counts ####################
wnn25 <- readRDS('data/wnn25.rds')
DefaultAssay(wnn25) <- 'ADT'
raw.counts = wnn25 %>% GetAssayData(slot='counts') %>% as.matrix()
wnn25 = wnn25 %>% NormalizeData(normalization.method='CLR', margin=2)
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)
all.types=unique(wnn25$celltype.l2)
#################### scPDA #########################
pi = h5read('results/wnn25/wnn25_scPDA.h5', 'pi') %>% as.matrix()
scPDA.counts = round(raw.counts * (1-pi))
wnn25.scPDA = SetAssayData(wnn25, slot='counts', new.data=scPDA.counts) %>% 
  NormalizeData(normalization.method='CLR', margin=2)
#################### scAR #########################
scAR.counts = h5read('results/wnn25/wnn25_scAR.h5', 'count')
colnames(scAR.counts)=colnames(raw.counts); rownames(scAR.counts)=rownames(raw.counts)
wnn25.scAR = SetAssayData(wnn25, slot='counts', new.data=scAR.counts) %>% 
  NormalizeData(normalization.method='CLR', margin=2)
#################### DSB #########################
dsb.counts =  readRDS('results/wnn25/wnn25_DSB.rds')
wnn25.dsb = SetAssayData(wnn25, slot='counts', new.data = (exp(dsb.counts)-0.5)) %>% 
  NormalizeData(normalization.method='CLR', margin=2)
#################### GMM #########################
gmm.counts = readRDS('results/wnn25/wnn25_GMM.rds')
wnn25.gmm = SetAssayData(wnn25, slot='counts', new.data=gmm.counts) %>% 
  NormalizeData(normalization.method='CLR', margin=2)
#################### DecontPro ###################
decont.counts = readr::read_csv('results/wnn25/wnn25_decontpro.csv')
decont.counts = column_to_rownames(decont.counts, var=colnames(decont.counts)[1]) %>% as.matrix()
rownames(decont.counts)=rownames(raw.counts); colnames(decont.counts)=colnames(raw.counts)
wnn25.decont = SetAssayData(wnn25, slot='counts', new.data=decont.counts) %>% 
  NormalizeData(normalization.method='CLR', margin=2)

###############################################################################
########################## Negative Score Compare Plot ########################
den.score = function(type, raw, denoise, pda, n.compare=8){
  preprocess = function(obj){
    GetAssayData(obj %>% subset(celltype.l2 %in% type),slot='data') %>% 
      t %>% as.data.frame() %>% rename(CD197=`CD197-CCR7`, CD278=`CD278-ICOS`, CD127=`CD127-IL7Ra`)
  }
  raw = preprocess(raw); denoise=preprocess(denoise); pda=preprocess(pda)
  levels = names(sort(colMeans(pda)))
  neg.prot = pda %>% colMeans() %>% sort() %>% head(n.compare) %>% names()
  d.score = mean((colMeans(raw)[neg.prot] - colMeans(denoise)[neg.prot])/colMeans(raw)[neg.prot])
  return(d.score)
}
ds.df = data.frame(DSB=NA, scAR=NA, GMM=rep(NA, length(all.types)), DecontPro=NA, scPDA=NA,row.names = all.types)
for (tp in all.types){
  ds.df[tp, 'scPDA'] = den.score(type=tp, raw=wnn25, denoise=wnn25.scPDA, pda=wnn25.scPDA)
  ds.df[tp, 'scAR'] = den.score(type=tp, raw=wnn25, denoise=wnn25.scAR, pda=wnn25.scPDA)
  ds.df[tp, 'DecontPro'] = den.score(type=tp, raw=wnn25, denoise=wnn25.decont, pda=wnn25.scPDA)
  ds.df[tp, 'GMM'] = den.score(type=tp, raw=wnn25, denoise=wnn25.gmm, pda=wnn25.scPDA)
  ds.df[tp, 'DSB'] = den.score(type=tp, raw=wnn25, denoise=wnn25.dsb, pda=wnn25.scPDA)
}

df_long.ds <- tidyr::gather(ds.df %>% rownames_to_column(var='CellType'), key = "method", value = "Negative Score", -CellType)
avg=apply(ds.df, 2, mean); med=apply(ds.df, 2, median)
df_long.ds = df_long.ds %>% mutate(avg.med=NA)
for (m in c('GMM', 'DSB', 'scAR', 'scPDA', 'DecontPro')){
  df_long.ds[df_long.ds$method==m, "avg.med"]= sprintf('(Avg: %.0f%%, Med: %.0f%%)', avg[m]*100, med[m]*100)
}
df_long.ds$method = paste(df_long.ds$method, df_long.ds$avg.med)
method.name = unique(df_long.ds$method)
df_long.ds$method = factor(df_long.ds$method, levels=method.name)
colors <- setNames(rgb_alpha(c("grey", "steelblue", "green4", "orange1", "red3"),0.6), method.name)
p.dscore = ggplot(df_long.ds, aes(x = CellType, y = `Negative Score`, color = method)) + geom_point(size = 2) + 
  theme_minimal() + box + labs(title='NRC for all cell types', x=NULL, y=NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.justification = 0.5, 
        legend.title = element_blank(), legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face='plain'))+ 
  scale_y_continuous(labels = label_percent(scale = 100)) + coord_cartesian(ylim = c(-1, 1)) +
  scale_color_manual(values = colors)

ggsave('fig_reprod/Figure2_g.png', plot=p.dscore, width=12, height=3, dpi=300)
###############################################################################
################################ Ridge Plots ##################################
cluster_ridge = function(type, raw, pda, title, n.compare=8, min_h=0.018, bold_labels = NULL){
  raw = GetAssayData(raw %>% subset(celltype.l2 %in% type),slot='data') %>% 
    t %>% as.data.frame() %>% rename(CD197=`CD197-CCR7`, CD278=`CD278-ICOS`, CD127=`CD127-IL7Ra`)
  pda = GetAssayData(pda %>% subset(celltype.l2 %in% type) ,slot='data') %>% 
    t %>% as.data.frame() %>% rename(CD197=`CD197-CCR7`, CD278=`CD278-ICOS`, CD127=`CD127-IL7Ra`)
  levels = names(sort(colMeans(pda)))
  p1 = ridge(raw, levels, min_h, bold_labels)+ggtitle('Raw Counts')+xlab(NULL)+ylab(NULL)+theme_PlainTitle
  p2=ridge(pda, levels, min_h, bold_labels)+ggtitle('scPDA denoised ')+xlab(NULL)+ylab(NULL)+theme_PlainTitle+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  neg.prot = pda %>% colMeans() %>% sort() %>% head(n.compare) %>% names()
  d.score = mean((colMeans(raw)[neg.prot] - colMeans(pda)[neg.prot])/colMeans(raw)[neg.prot])
  ridge_cluster = p1+p2+plot_layout(ncol = 2)+
    plot_annotation(title = sprintf("%s (NRC: %.2f)", title, d.score),
                    caption = "CLR Expression Levels") &
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  return(ridge_cluster)
}
type=c('NK')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='NK', bold_labels=c('CD161','CD3'))
ggsave("fig_reprod/Ridge_NK.png", ridge_cluster, width = 4, height = 4, dpi = 500)

type=c('CD4 Memory')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='CD4 Memory', bold_labels=c('CD45RA','CD45RO'))
ggsave("fig_reprod/Ridge_MemoryCD4.png", ridge_cluster, width = 4, height = 4, dpi = 500)

type=c('MAIT')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='MAIT')
ggsave("fig_reprod/Ridge_MAIT.png", ridge_cluster, width = 4, height = 4, dpi = 500)

type=c('CD14 Mono')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='CD14 Mono')
ggsave("fig_reprod/Ridge_Mono14.png", ridge_cluster, width = 4, height = 4, dpi = 500)

type=c('CD16 Mono')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='CD16 Mono')
ggsave("fig_reprod/Ridge_Mono16.png", ridge_cluster, width = 4, height = 4, dpi = 500)

type=c('pDC')
ridge_cluster=cluster_ridge(type, wnn25, wnn25.scPDA, title='plasmacytoid DC')
ggsave("fig_reprod/Ridge_pDC.png", ridge_cluster, width = 4, height = 4, dpi = 500)

##################### Combine Plots #####################
img1 <- image_read("fig_reprod/Ridge_NK.png")
img2 <- image_read("fig_reprod/Ridge_MemoryCD4.png")
img3 <- image_read("fig_reprod/Ridge_MAIT.png")
img4 <- image_read("fig_reprod/Ridge_Mono14.png")
img5 <- image_read("fig_reprod/Ridge_Mono16.png")
img6 <- image_read("fig_reprod/Ridge_pDC.png")
img7 <- image_read("fig_reprod/Figure2_g.png") %>% image_background(color='white')

img_list <- lapply(list(img1, img2, img3, img4, img5, img6), function(img) {
  image_resize(img, "800x800")
})
row1 <- image_append(c(img_list[[1]], img_list[[2]], img_list[[3]]))
row2 <- image_append(c(img_list[[4]], img_list[[5]], img_list[[6]]))
row3 <- image_resize(img7, "2400x")
combined_image <- image_append(c(row1, row2, row3), stack = TRUE)
image_write(combined_image, path = "fig_reprod/Figure2.jpg", format='jpeg', quality=60)
file.remove(
  "fig_reprod/Ridge_NK.png",
  "fig_reprod/Ridge_MemoryCD4.png",
  "fig_reprod/Ridge_MAIT.png",
  "fig_reprod/Ridge_Mono14.png",
  "fig_reprod/Ridge_Mono16.png",
  "fig_reprod/Ridge_pDC.png",
  "fig_reprod/Figure2_g.png"
)