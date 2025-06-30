library(here)
setwd(here())
source('code/tools.R')

theme.border = theme(panel.background = element_rect(colour = "black", fill = NA, linewidth = 1))
title.center = theme(plot.title = element_text(hjust = 0.5, face='plain'))
box = annotation_custom(grob = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 2)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
spacer <- plot_spacer() + theme_void() + theme(plot.background = element_blank())

tea <- readRDS('data/teaseq.rds')
empty = tea$empty; tea=tea$teaseq; tea=tea%>%subset(celltype != 'Unknown')
raw.counts = tea %>% GetAssayData(slot='counts') %>% as.matrix()
cell_log = raw.counts %>% log1p
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)
type_vec=tea$celltype
################################# scPDA #################################
pi = h5read('results/tea/tea_scPDA.h5', 'pi') %>% as.matrix()
colnames(pi)=colnames(raw.counts); rownames(pi)=rownames(raw.counts)
scPDA.counts = round(raw.counts * (1-pi))
################################# DSB #################################
dsb.counts = readRDS('results/tea/tea_dsb.rds')
################################# scAR #################################
scar.counts = h5read('results/tea/tea_scAR.h5', 'count') %>% as.matrix()
rownames(scar.counts) = prot.names; colnames(scar.counts) = cell.names
################################# DecontPro #################################
decont.counts = readr::read_csv('results/tea/tea_decontpro.csv')
decont.counts = column_to_rownames(decont.counts, var=colnames(decont.counts)[1]) %>% as.matrix()
rownames(decont.counts)=rownames(raw.counts); colnames(decont.counts)=colnames(raw.counts)

###########################################################################
################################ Figure a-c ################################
p_list = list()
prot='CD3'; pos=c('CD4 Memory', 'CD4 Naive', 'CD8 Memory', 'CD8 Naive', 'MAIT'); pos.name='T Cells'
roc.list = lapply(list(Raw=raw.counts,DSB=dsb.counts, DecontPro=decont.counts,scAR=scar.counts), 
                  function(x) assign_prob(counts=x, prot=prot,pos=pos,types=type_vec))
roc.list[['scPDA']] = roc_auc(pi, prot=prot,pos=pos,types=type_vec)
roc_vec = sapply(roc.list, function(x) as.numeric(x$auc))
p1 = overlap.gg(dsb.counts, prot=prot, pos=pos, title='DSB', pos.name=pos.name, types=type_vec, xlab='DSB Value', AUC=roc_vec['DSB']) + NoLegend() + labs(y=NULL)
p2 = overlap.gg(decont.counts, prot=prot, pos=pos, title='DecontPro', pos.name=pos.name, types=type_vec, AUC=roc_vec['DecontPro']) + NoLegend() + labs(y=NULL)
p3 = overlap.gg(scar.counts, prot=prot, pos=pos, title='scAR', pos.name=pos.name, types=type_vec, AUC=roc_vec['scAR']) + NoLegend() + labs(y=NULL)
p4 = overlap.gg(scPDA.counts, prot=prot, pos=pos, title='scPDA', pos.name=pos.name, types=type_vec, AUC=roc_vec['scPDA']) + NoLegend() + labs(y=NULL)
p1_with_leg = overlap.gg(raw.counts, prot=prot, pos=pos, title='Raw', pos.name=pos.name, types=type_vec, AUC=roc_vec['RAW']) + theme(legend.position = "top",  legend.justification = "center", legend.box = "horizontal",legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8)) +guides(fill=guide_legend(title=NULL))
g <- ggplotGrob(p1_with_leg); legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]
hist.4 = grid.arrange(p1, p2, p3, p4, ncol=2)
CD3 = grid.arrange(legend_grob, hist.4, nrow=2, heights=c(0.05,1), top=textGrob(prot, gp=gpar(fontsize=16, fontface="bold")), left='Density')
p_list[[prot]] = CD3

prot='CD19'; pos=c('B'); pos.name='B Cells'
roc.list = lapply(list(Raw=raw.counts,DSB=dsb.counts, DecontPro=decont.counts,scAR=scar.counts), 
                  function(x) assign_prob(counts=x, prot=prot,pos=pos,types=type_vec))
roc.list[['scPDA']] = roc_auc(pi, prot=prot,pos=pos,types=type_vec)
roc_vec = sapply(roc.list, function(x) as.numeric(x$auc))
p1 = overlap.gg(dsb.counts, prot=prot, pos=pos, title='DSB', pos.name=pos.name, types=type_vec, xlab='DSB Value', AUC=roc_vec['DSB']) + NoLegend() + labs(y=NULL)
p2 = overlap.gg(decont.counts, prot=prot, pos=pos, title='DecontPro', pos.name=pos.name, types=type_vec, AUC=roc_vec['DecontPro']) + NoLegend() + labs(y=NULL)
p3 = overlap.gg(scar.counts, prot=prot, pos=pos, title='scAR', pos.name=pos.name, types=type_vec, AUC=roc_vec['scAR']) + NoLegend() + labs(y=NULL)
p4 = overlap.gg(scPDA.counts, prot=prot, pos=pos, title='scPDA', pos.name=pos.name, types=type_vec, AUC=roc_vec['scPDA']) + NoLegend() + labs(y=NULL)
p1_with_leg = overlap.gg(raw.counts, prot=prot, pos=pos, title='Raw', pos.name=pos.name, types=type_vec, AUC=roc_vec['Raw']) + theme(legend.position = "top",  legend.justification = "center", legend.box = "horizontal",legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8)) +guides(fill=guide_legend(title=NULL))
g <- ggplotGrob(p1_with_leg); legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]
hist.4 = grid.arrange(p1, p2, p3, p4, ncol=2)
CD19 = grid.arrange(legend_grob, hist.4, nrow=2, heights=c(0.05,1), top=textGrob(prot, gp=gpar(fontsize=16, fontface="bold")), left='Density')
p_list[[prot]] = CD19

prot='CD4'; pos1=c("CD4 Memory","CD4 Naive"); pos2=c("CD14 Mono", "CD16 Mono", "Dendritic")
neg=c("B", "NK", "CD8 Naive", "CD8 Memory"); pos.name=c('CD4 T cells', 'Monocytes/DC')
p1=overlap3.gg(dsb.counts, prot=prot, types=type_vec, pos1=pos1, pos2 = pos2, neg = neg, title='DSB', pos.name=pos.name, xlab='DSB Value') + NoLegend() + labs(y=NULL)
p2=overlap3.gg(decont.counts, prot=prot, types=type_vec, pos1=pos1, pos2 = pos2, neg = neg, title='Raw', pos.name=pos.name) + NoLegend() + labs(y=NULL)
p3=overlap3.gg(scar.counts, prot=prot, types=type_vec, pos1=pos1, pos2 = pos2, neg = neg, title='scAR', pos.name=pos.name) + NoLegend() + labs(y=NULL)
p4=overlap3.gg(scPDA.counts, prot=prot, types=type_vec, pos1=pos1, pos2 = pos2, neg = neg, title='scPDA', pos.name=pos.name) + NoLegend() + labs(y=NULL)
p1_with_leg = overlap3.gg(raw.counts, prot=prot, types=type_vec, pos1=pos1, pos2 = pos2, neg = neg, title='Raw', pos.name=pos.name) + theme(legend.position = "top",  legend.justification = "center", legend.box = "horizontal",legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8)) +guides(fill=guide_legend(title=NULL))
g <- ggplotGrob(p1_with_leg); legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]
hist.4 = grid.arrange(p1, p2, p3, p4, ncol=2)
CD4 = grid.arrange(legend_grob, hist.4, nrow=2, heights=c(0.05,1), top=textGrob(prot, gp=gpar(fontsize=16, fontface="bold")), left='Density')
p_list[[prot]] = CD4

combined_plot <- grid.arrange( grobs = p_list, ncol = 3)
ggsave("fig_reprod/Figure4_a-c.png", combined_plot, width = 18, height = 6, dpi = 200)
################################ Figure d-f ################################
plots = function(raw, dsb, decont, scar, scPDA, pi, prot, pos, pos.name, types){
  p1=overlap.gg(raw, prot, pos, pos.name, title='Raw', types) + NoLegend()
  p2=overlap.gg(dsb, prot, pos, pos.name, title='DSB', types, xlab='DSB Value') + NoLegend()
  p3=overlap.gg(decont, prot, pos, pos.name, title='DecontPro', types) + NoLegend()  
  p4=overlap.gg(scar, prot, pos, pos.name, title='scAR', types) + NoLegend()
  p5=overlap.gg(scPDA, prot, pos, pos.name, title='scPDA', types) + NoLegend()
  p6=overlap.gg(expm1(pi), prot, pos, pos.name, title='scPDA Background Probability', 
                types, ifLFC=FALSE, xlab='scPDA Background Probability') + NoLegend()
  roc.list = lapply(list(Raw=raw,DSB=dsb, DecontPro=decont,scAR=scar), 
                    function(x) assign_prob(counts=x, prot=prot,pos=pos,types=types))
  roc.list[['scPDA']] = roc_auc(pi, prot=prot,pos=pos,types=types)
  p7=curve.gg(roc.list)
  
  # Fetch the Legend
  p1_with_leg = overlap.gg(raw, prot, pos, pos.name, title='Raw', types) + 
    theme(legend.position = "top",  legend.justification = "center", legend.box = "horizontal",
          legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8))+
    guides(fill=guide_legend(title=NULL))
  g <- ggplotGrob(p1_with_leg); legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]
  
  p.all = grid.arrange(p1,p2,p3,p4,p5,p6,p7, nrow=1)
  p.all = grid.arrange(legend_grob, p.all, heights=c(0.05,1), 
                       top=textGrob(prot, gp=gpar(fontsize=16, fontface="bold")))
  return(p.all)
}
p_list = list()
prot='CD40'; pos=c('B', 'Dendritic', 'CD16 Mono', 'CD14 Mono'); pos.name='DC/B/Mono'
p = plots(raw.counts, dsb.counts, decont.counts, scar.counts, scPDA.counts, pi, prot, pos, pos.name, types=type_vec)
p_list[[prot]] = p

prot='IgD'; pos=c('B'); pos.name='B'
p = plots(raw.counts, dsb.counts, decont.counts, scar.counts, scPDA.counts, pi, prot, pos, pos.name, types=type_vec)
p_list[[prot]] = p

prot='CD192'; pos=c('CD14 Mono'); pos.name='CD14 Mono'
p = plots(raw.counts, dsb.counts, decont.counts, scar.counts, scPDA.counts, pi, prot, pos, pos.name, types=type_vec)
p_list[[prot]] = p

combined_plot <- grid.arrange( grobs = p_list, ncol = 1)
ggsave("fig_reprod/Figure4_d-f.png", combined_plot, width = 20, height = 9, dpi = 200)
################################ All Plots ########################################
plots = function(raw, dsb, decont, scar, scPDA, prot, pos, pos.name, types){
  par(mfrow=c(1,7))
  lfc.raw = overlap(raw, prot, pos, pos.name, title='Raw', types)
  lfc.dsb = overlap(dsb, prot, pos, pos.name, title='DSB', types, xlab='DSB Value')
  lfc.decont = overlap(decont, prot, pos, pos.name, title='DecontPro', types)  
  lfc.scar = overlap(scar, prot, pos, pos.name, title='scAR', types)
  lfc.scPDA = overlap(scPDA, prot, pos, pos.name, title='scPDA', types)
  overlap(pi, prot, pos, pos.name, title='scPDA Background Probability', 
          types, ifprob=TRUE, xlab='scPDA Background Probability')
  roc.list = lapply(list(Raw=raw,DSB=dsb, DecontPro=decont, scAR=scar), 
                    function(x) assign_prob(counts=x, prot=prot,pos=pos,types=types))
  roc.list[['scPDA']] = roc_auc(pi, prot=prot,pos=pos,types=types)
  curve(roc.list, prot)
  auc=sapply(roc.list, function(x) x$auc)
  return(list(Raw=lfc.raw, DSB=lfc.dsb, DecontPro=lfc.decont, scAR=lfc.scar, scPDA=lfc.scPDA, AUC=auc))
}
pdf('fig_reprod/supp/full_list_bg.pdf', width=18, height=3, title='Background Probability for Marker Proteins')
lfc.df = data.frame(Raw=numeric(0), DSB=numeric(0), scAR=numeric(0), scPDA=numeric(0)); auc.df = lfc.df
par(mfrow=c(23,1))
marker.list = list()
marker=list(prot='CD14', pos=c('CD14 Mono', 'Dendritic'), pos.name='DC/CD14 Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD11c', pos=c('CD14 Mono', 'CD16 Mono', 'NK','Dendritic'), pos.name='DC/Mono/NK'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD123', pos=c('CD14 Mono', 'CD16 Mono', 'B','Dendritic'), pos.name='DC/Mono/B'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD141', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD16', pos=c('CD16 Mono', 'Dendritic', 'NK'), pos.name='DC/CD16 Mono/NK'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD172a', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD185', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD19', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD192', pos=c('CD14 Mono'), pos.name='CD14 Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD21', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD24', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD27', pos=c('CD4 Memory', 'CD4 Naive', 'CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='T Cells'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD3', pos=c('CD4 Memory', 'CD4 Naive', 'CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='T Cells'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD40', pos=c('B', 'Dendritic'), pos.name='B/DC'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD56', pos=c('MAIT', 'NK'), pos.name='NK/MAIT'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD86', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic'), pos.name='DC/Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='CD8a', pos=c('CD8 Memory', 'CD8 Naive', 'MAIT'), pos.name='CD8 T/MAIT'); marker.list = append(marker.list, list(marker))
marker=list(prot='FceRI', pos=c('Dendritic'), pos.name='Dendritic'); marker.list = append(marker.list, list(marker))
marker=list(prot='HLA-DR', pos=c('CD14 Mono', 'CD16 Mono', 'Dendritic', 'B'), pos.name='B/DC/Mono'); marker.list = append(marker.list, list(marker))
marker=list(prot='IgD', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='IgM', pos=c('B'), pos.name='B'); marker.list = append(marker.list, list(marker))
marker=list(prot='TCR-Va7.2', pos=c('MAIT'), pos.name='MAIT'); marker.list = append(marker.list, list(marker))
for (marker in marker.list){
  prot = marker$prot; pos=marker$pos; pos.name=marker$pos.name
  result = plots(raw=raw.counts, dsb=dsb.counts, decont=decont.counts, scar=scar.counts, 
                 scPDA=scPDA.counts, prot=prot, pos=pos, pos.name=pos.name, types=type_vec)
  prot.lfc = data.frame(Raw=result$Raw, DSB=result$DSB, DecontPro=result$DecontPro, scAR=result$scAR, scPDA=result$scPDA, row.names=prot); 
  prot.auc = data.frame(Raw=NA, DSB=NA, DecontPro=NA, scAR=NA, scPDA=NA, row.names=prot)
  prot.auc[1,]=as.numeric(result$AUC)
  lfc.df=rbind(lfc.df, prot.lfc); auc.df=rbind(auc.df, prot.auc)
}
dev.off()
################################################################################
########################### Fig g-i ############################################
color_values = c("Raw" = "grey", "DSB"="steelblue", "DecontPro"="orange2", "scAR"="pink3", "scPDA"="red3")
df_long.lfc <- tidyr::gather(lfc.df %>% rownames_to_column(var='protein'), key = "method", value = "expression", -protein)
df_long.lfc$method <- factor(df_long.lfc$method, levels = c("Raw", "DSB", "DecontPro", "scAR", "scPDA"))
p.lfc = ggplot(df_long.lfc, aes(x = protein, y = expression, color = method)) + geom_point(size = 2) + 
  theme_minimal() + title.center + box + coord_flip() + 
  labs(title='LFC Comparison', x='Proteins', y='Log Fold Change') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top", legend.justification = 0.5, 
        legend.title = element_blank(), legend.background = element_blank())+
  scale_color_manual(values = color_values)+
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

df_long.auc <- tidyr::gather(auc.df %>% rownames_to_column(var='protein'), key = "method", value = "expression", -protein)
df_long.auc$method <- factor(df_long.auc$method, levels = c("Raw", "DSB", "DecontPro", "scAR", "scPDA"))
p.auc = ggplot(df_long.auc, aes(x = protein, y = expression, color = method)) + geom_point(size = 2) + 
  theme_minimal() + title.center + box + coord_flip() + 
  labs(title='AUC Comparison', x='Proteins', y='AUC') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top", legend.justification = 0.5, 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), legend.background = element_blank())+
  scale_color_manual(values = color_values)+
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

p_list = list()
p_list[['LFC']] = p.lfc
p_list[['AUC']] = p.auc
combined_plot <- grid.arrange( grobs = p_list, ncol = 2)
ggsave("fig_reprod/Figure4_g-h.png", combined_plot, width = 6, height = 8, dpi = 300)

################################################################################
########################### Fig i ##############################################
set.seed(2)
prot='CD40'; pos=c('B', 'Dendritic', 'CD14 Mono', 'CD16 Mono'); pos.name='DC/B/Mono'
prot_raw=log1p(raw.counts[prot,]); prot_pi=pi[prot,]
gmm = Mclust(prot_raw, G=2, warn = FALSE, verbose = FALSE)
raw.foreground = gmm$z[,2,drop=FALSE] %>% t()
rownames(raw.foreground)=prot

prot_leg = overlap.gg(raw.foreground, prot, pos, pos.name, title='Raw', types=type_vec)
prob_dens = overlap.gg(1-pi, prot, pos, pos.name, title='1 - Background Probability', 
                       types=type_vec, ifLFC=FALSE, xlab='1 - Background Probability')+
  NoLegend() + coord_flip() + cowplot::theme_nothing()
prot_dens = prot_leg + ggtitle('CD40') + theme_void() + NoLegend() + title.center
g <- ggplotGrob(prot_leg)
legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]

prot_raw=log1p(raw.counts[prot,]); prot_pi=pi[prot,]
annot = data.frame(x=rep(c(0,1),2), y=c(0,0,1,1),t=c("A","B","C","D"))
df = data.frame(raw_prob=raw.foreground[1,], scPDA_prob=1-prot_pi, type=tea$celltype)
df$section <- with(df, factor( 
  ifelse(raw_prob < 0.5 & scPDA_prob < 0.5, "SectionA",
         ifelse(raw_prob >= 0.5 & scPDA_prob < 0.5, "SectionB",
                ifelse(raw_prob < 0.5 & scPDA_prob >= 0.5, "SectionC", "SectionD")))))
scatter=ggplot(df)+geom_point(aes(x=raw_prob, y=scPDA_prob, color=section), size=0.2)+
  scale_color_manual(values=c("SectionA"="green4", "SectionB"="red3", "SectionC"="purple3", "SectionD"="steelblue")) +
  cowplot::theme_nothing()+xlab('1 - GMM Background Probability')+ylab('1 - scPDA Background Probability')+
  geom_hline(yintercept = c(0.5), linetype='dotted', linewidth=1)+
  geom_vline(xintercept = c(0.5), linetype='dotted', linewidth=1)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title.x = element_text(face="plain", color="black", size=11),
        axis.title.y = element_text(face="plain", color="black", size=11,angle=90, vjust=0.5))+
  geom_text(data=annot, aes(x=x,y=y,label=t), color='black', size=10, fontface="bold")

scatter = grid.arrange(prot_dens, legend_grob, scatter, prob_dens, widths=c(6,1), heights=c(1,6))

p1 = DimPlot(tea, reduction='CLR_wnn_umap', group.by='celltype', 
             cells.highlight = rownames(df %>% filter(section == 'SectionA')),
             cols.highlight = 'green4', pt.size=0.2, sizes.highlight=0.05, label=TRUE)+
  theme_void() + NoLegend() + theme.border + ggtitle('Region A') + title.center

p2 = DimPlot(tea, reduction='CLR_wnn_umap', group.by='celltype', 
             cells.highlight = rownames(df %>% filter(section == 'SectionB')),
             cols.highlight = 'red3', pt.size=0.2, sizes.highlight=0.05, label=TRUE)+
  theme_void() + NoLegend() + theme.border + ggtitle('Region B') + title.center

p3 = DimPlot(tea, reduction='CLR_wnn_umap', group.by='celltype', 
             cells.highlight = rownames(df %>% filter(section == 'SectionC')),
             cols.highlight = 'purple3', pt.size=0.2, sizes.highlight=0.05, label=TRUE)+
  theme_void() + NoLegend() + theme.border + ggtitle('Region C') + title.center

p4 = DimPlot(tea, reduction='CLR_wnn_umap', group.by='celltype', 
             cells.highlight = rownames(df %>% filter(section == 'SectionD')),
             cols.highlight = 'steelblue', pt.size=0.2, sizes.highlight=0.05, label=TRUE)+
  theme_void() + NoLegend() + theme.border + ggtitle('Region D') + title.center

umaps = grid.arrange(p3, spacer, p4, p1, spacer, p2, widths=c(20,1,20), heights=c(20,20))

all = grid.arrange(scatter, spacer, umaps, widths=c(10,1,10))
ggsave('fig_reprod/Figure4_i.png', plot=all, width=16, height=8, dpi=300)


################################################################################
########################### Combine ############################################
row1 <- image_read("fig_reprod/Figure4_a-c.png")
target_width <- image_info(row1)$width

row2 <- image_read("fig_reprod/Figure4_d-f.png")
row2 <- image_resize(row2, paste0(target_width, "x"))

img3 <- image_read("fig_reprod/Figure4_g-h.png")
img4 <- image_read("fig_reprod/Figure4_i.png")
ratio_d = image_info(img3)$width/(image_info(img3)$width+image_info(img4)$width)
img3 = image_resize(img3, paste0(target_width*ratio_d, "x"))
img4 = image_resize(img4, paste0(target_width*(1-ratio_d), "x"))
row3 <- image_append(c(img3, img4))

final_image <- image_append(c(row1, row2, row3), stack = TRUE)
image_write(final_image, path = "fig_reprod/Figure4.jpg", format="jpeg", quality=60)

file.remove("fig_reprod/Figure4_a-c.png",
            "fig_reprod/Figure4_d-f.png",
            "fig_reprod/Figure4_g-h.png",
            "fig_reprod/Figure4_i.png")