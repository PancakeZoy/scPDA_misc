library(here)
setwd(here())
source('code/tools.R')
set.seed(1)
box = annotation_custom(grob = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 2)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
######################### Ground Truth Noise #########################
dsb = readRDS('data/dsb.rds');
raw.counts = dsb %>% GetAssayData(slot='counts') %>% as.matrix()
cell.counts = dsb %>% subset(Type1 != 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
unstain.counts = dsb %>% subset(Type1 == 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); unstain.names=colnames(unstain.counts); cell.names=colnames(cell.counts)
true.noise = rowMeans(unstain.counts)

######################### scPDA Noise #########################
pi = h5read('results/dsb/dsb_scPDA.h5', 'pi') %>% as.matrix()
colnames(pi)=colnames(raw.counts); rownames(pi)=rownames(raw.counts)
scPDA.counts=round(raw.counts * (1-pi)); scPDA.unstain=scPDA.counts[, unstain.names]
scPDA.noise = rowMeans(scPDA.unstain)
scPDA.prop = (true.noise-scPDA.noise)/true.noise
######################### DSB Noise ###########################
dsb.counts = readRDS('results/dsb/dsb_DSB.rds')
dsb.unstain = dsb.counts[, unstain.names] %>% expm1()
dsb.noise = rowMeans(dsb.unstain)
dsb.prop = (true.noise-dsb.noise)/true.noise
######################### scAR Noise ###########################
scAR.counts = h5read('results/dsb/dsb_scAR.h5', 'count') %>% as.matrix()
colnames(scAR.counts)=colnames(raw.counts); rownames(scAR.counts)=rownames(raw.counts)
scAR.unstain=scAR.counts[, unstain.names]
scAR.noise = rowMeans(scAR.unstain)
scAR.prop = (true.noise-scAR.noise)/true.noise
######################### GMM Noise ###########################
gmm.counts = readRDS('results/dsb/dsb_GMM.rds')
gmm.unstain=gmm.counts[, unstain.names]
gmm.noise = rowMeans(gmm.unstain)
gmm.prop = (true.noise-gmm.noise)/true.noise

######################### DecontPro Noise ###########################
decont.counts = readr::read_csv('results/dsb/dsb_DecontPro.csv')
decont.counts = column_to_rownames(decont.counts, var=colnames(decont.counts)[1]) %>% as.matrix()
rownames(decont.counts)=rownames(raw.counts); colnames(decont.counts)=colnames(raw.counts)
decont.unstain=decont.counts[, unstain.names]
decont.noise = rowMeans(decont.unstain)
decont.prop = (true.noise-decont.noise)/true.noise

######################### Denoising Proportion For Neg Control ###########################
prop.df = data.frame(GMM=gmm.prop, DSB=dsb.prop, scAR=scAR.prop, DecontPro=decont.prop, scPDA=scPDA.prop)
df_long.prop <- tidyr::gather(prop.df %>% rownames_to_column(var='protein'), key = "method", value = "NRC", -protein)
avg=apply(prop.df, 2, mean); med=apply(prop.df, 2, median)
df_long.prop = df_long.prop %>% mutate(avg.med=NA)
for (m in c('GMM', 'DSB','scAR', 'DecontPro', 'scPDA')){
  df_long.prop[df_long.prop$method==m, "avg.med"]= sprintf('(Avg: %.0f%%, Med: %.0f%%)', avg[m]*100, med[m]*100)
}
df_long.prop$method = paste(df_long.prop$method, df_long.prop$avg.med)
method.name = unique(df_long.prop$method)
df_long.prop$method = factor(df_long.prop$method, levels=method.name)
colors <- setNames(rgb_alpha(c("grey", "steelblue","orange1", "pink3", "red3"), 0.6), method.name)
df_long.prop$protein[grep('type', df_long.prop$protein)] = paste0('Isotype', rep(1:4,length(method.name)))
p.prop = ggplot(df_long.prop, aes(x = protein, y = `NRC`, color = method)) + geom_point(size = 2) + 
  theme_minimal() + box + labs(title='NRC for all proteins across unstained cells', x=NULL, y=NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.justification = 0.5, 
        legend.title = element_blank(), legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face='plain'))+ 
  scale_y_continuous(labels = label_percent(scale = 100)) + coord_cartesian(ylim = c(-1, 1)) +
  scale_color_manual(values = colors)

#################### Denoising Proportion For Non-expressive Proteins ##################
prot.mean = function(prot){
  raw.mean=raw.counts[prot, cell.names] %>% mean()
  gmm.mean=gmm.counts[prot, cell.names] %>% mean()
  dsb.mean=dsb.counts[prot, cell.names] %>% mean()
  scAR.mean=scAR.counts[prot, cell.names] %>% mean()
  decont.mean=decont.counts[prot, cell.names] %>% mean()
  scPDA.mean=scPDA.counts[prot, cell.names] %>% mean()
  return(data.frame(GMM=(raw.mean-gmm.mean)/raw.mean, 
                    DSB=(raw.mean-dsb.mean)/raw.mean,
                    scAR=(raw.mean-scAR.mean)/raw.mean, 
                    DecontPro = (raw.mean-decont.mean)/raw.mean,
                    scPDA=(raw.mean-scPDA.mean)/raw.mean, 
                    row.names = prot))
}
neg.prots = c('CD117', 'CD137', 'CD138', 'CD206', 'CD223',
              'CD273', 'CD80', "MouseIgG1kappaisotype",
              "MouseIgG2akappaisotype", "Mouse IgG2bkIsotype", "RatIgG2bkIsotype")
homo = data.frame()
for (p in neg.prots){homo = rbind(homo, prot.mean(p))}
rownames(homo)[grep('type', neg.prots)] = paste0('Isotype', 1:4)
df_long.homo <- tidyr::gather(homo %>% rownames_to_column(var='protein'), key = "method", value = "Expression", -protein)
avg=apply(homo, 2, mean); med=apply(homo, 2, median)
df_long.homo = df_long.homo %>% mutate(avg.med=NA)
for (m in c('GMM', 'DSB', 'scAR', 'DecontPro', 'scPDA')){
  df_long.homo[df_long.homo$method==m, "avg.med"]= sprintf('(Avg: %.0f%%, Med: %.0f%%)', avg[m]*100, med[m]*100)
}
df_long.homo$method = paste(df_long.homo$method, df_long.homo$avg.med)
method.name = unique(df_long.homo$method)
df_long.homo$method = factor(df_long.homo$method, levels=method.name)
colors <- setNames(rgb_alpha(c('green4',"steelblue", "orange1", "pink3", "red3"),0.6), method.name)
p.homo = ggplot(df_long.homo, aes(x = protein, y = Expression, color = method)) + geom_point(size = 2) + 
  theme_minimal() + box + labs(title='NRC for non-expressive proteins', x=NULL, y=NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.justification = 0.5, 
        legend.title = element_blank(), legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face='plain'))+ 
  scale_y_continuous(labels = label_percent(scale = 100)) + coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = colors)

plots_list = list()
plots_list[['prop']] = p.prop; 
plots_list[['homo']] = p.homo
combined_plot <- grid.arrange( grobs = plots_list, ncol = 1)
ggsave("fig_reprod/Figure1_j-k.png", combined_plot, width = 12, height = 6, dpi = 200)