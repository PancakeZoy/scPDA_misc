library(here)
setwd(here())
source('code/tools.R')

dsb_dir = 'results/titr188/titr188_dsb.csv'
gmm_dir = 'results/titr188/titr188_gmm.csv'
decontpro_dir = 'results/titr188/titr188_decontpro.csv'
scar_dir = 'results/titr188/titr188_scar.csv'
scpda_dir = 'results/titr188/titr188_scpda.csv'

theme.border = theme(panel.background = element_rect(colour = "black", fill = NA, linewidth = 1))
title.center = theme(plot.title = element_text(hjust = 0.5, face='plain'))
box = annotation_custom(grob = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 2)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

####################################################################################
# Load Raw
Titr188 = readRDS('data/titr188.rds') %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
Idents(Titr188) = 'type'
raw.counts = Titr188 %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)
# Load scPDA
scPDA.counts = read.csv(scpda_dir, row.names = 1, check.names = FALSE) %>% as.matrix() %>% round()
Titr188.scPDA = SetAssayData(Titr188, slot='counts', new.data=scPDA.counts) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
# Load DSB
dsb.counts =  read.csv(dsb_dir, row.names = 1, check.names = FALSE) %>% as.matrix()
Titr188.DSB = SetAssayData(Titr188, slot='counts', new.data=expm1(dsb.counts)) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
# Load scAR
scar.counts = read.csv(scar_dir, row.names = 1, check.names = FALSE) %>% as.matrix()
rownames(scar.counts) = prot.names; colnames(scar.counts) = cell.names
Titr188.scAR = SetAssayData(Titr188, slot='counts', new.data=scar.counts) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
# Load GMM
gmm.counts = read.csv(gmm_dir, row.names = 1, check.names = FALSE) %>% as.matrix() %>% round()
Titr188.gmm = SetAssayData(Titr188, slot='counts', new.data=gmm.counts) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
# Load DecontPro
decont.counts = read.csv(decontpro_dir, row.names = 1, check.names = FALSE) %>% as.matrix()
Titr188.DecontPro = SetAssayData(Titr188, slot='counts', new.data=decont.counts) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)
#################################### Figure a ####################################
mat.cut <- matrix(c(1, 1, 1, 1, 4, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3), nrow=5, byrow=TRUE)
mat.all <- matrix(1:6, ncol=3)
spacer <- plot_spacer() + theme_void() + theme(plot.background = element_blank())

type_leg <- FeatureScatter(Titr188.scPDA %>% subset(type!='remaining'), feature1='CD4', feature2='CD14', jitter = TRUE) + 
  labs(color = NULL)+ guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.text = element_text(size = 8),
        legend.key.height = unit(0.5, "line"))
type_leg <- ggplotGrob(type_leg)
type_leg <- type_leg$grobs[[which(type_leg$layout$name == "guide-box")]]

cut1 = FeatureScatter(Titr188.scPDA %>% subset(type!='remaining'), feature1='CD4', feature2='CD14',jitter = TRUE)+
  geom_vline(xintercept = 2, linetype="dashed")+
  geom_hline(yintercept = 0.3, linetype="dashed")+
  cowplot::theme_nothing()+xlab('CD4')+ylab('CD14')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(face="plain", color="black", size=11),
        axis.title.y = element_text(face="plain", color="black", size=11,angle=90, vjust=0.5))
data.pda = Titr188.scPDA %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
x_density1 <- ggplot(data.pda, aes(x=CD4)) + geom_density(fill="steelblue", bw = 0.1) + 
  theme_void() + ggtitle('scPDA Denoised') + theme(plot.title = element_text(hjust = 0.5, vjust=-5, face='plain'))
y_density1 <- ggplot(data.pda, aes(x=CD14)) + geom_density(fill="orange", bw = 0.05) + theme_void() +coord_flip()
cut1 = grid.arrange(x_density1, cut1, y_density1, type_leg, layout_matrix=mat.cut)

cut2 = FeatureScatter(Titr188 %>% subset(type!='remaining'), feature1='CD4', feature2='CD14', jitter = TRUE)+NoLegend()+
  geom_vline(xintercept = 2, linetype="dashed")+ annotate("text", x = 3.8, y = 0.28, label = "455 CD4T Cells") +
  geom_hline(yintercept = 0.35, linetype="dashed")+
  cowplot::theme_nothing()+xlab('CD4')+ylab('CD14')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(face="plain", color="black", size=11),
        axis.title.y = element_text(face="plain", color="black", size=11,angle=90, vjust=0.5))
data.raw = Titr188 %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
x_density2 <- ggplot(data.raw, aes(x=CD4)) + geom_density(fill="steelblue", bw = 0.1) + 
  theme_void() +ggtitle('Raw Counts') + theme(plot.title = element_text(hjust = 0.5, vjust=-5, face='plain'))
y_density2 <- ggplot(data.raw, aes(x=CD14)) + geom_density(fill="orange", bw = 0.05) + theme_void() +coord_flip()
cut2 = grid.arrange(x_density2, cut2, y_density2, type_leg, layout_matrix=mat.cut)

LFC = function(raw, denoise, protein, pos, neg=NA){
  if (all(is.na(neg))){
    types = c("NK","CD8T","CD4T","B","CM")
    neg = setdiff(types, pos)
  }
  raw.pos = raw %>% filter(type %in% pos) %>% pull(protein) %>% mean()
  raw.neg = raw %>% filter(type %in% neg) %>% pull(protein) %>% mean()
  raw.lfc = log(raw.pos/raw.neg, 2)
  denoise.pos = denoise %>% filter(type %in% pos) %>% pull(protein) %>% mean()
  denoise.neg = denoise %>% filter(type %in% neg) %>% pull(protein) %>% mean()
  denoise.lfc = log(denoise.pos/denoise.neg, 2)
  lfc = data.frame(protein=protein, scPDA=denoise.lfc, Raw=raw.lfc)
  return(lfc)
}
cd4t.cd4=LFC(raw=data.raw, denoise=data.pda, protein='CD4', pos=c('CD4T'), neg=c('B','NK','CD8T'))
cm.cd4=LFC(raw=data.raw, denoise=data.pda, protein='CD4', pos=c('CM'), neg=c('B','NK','CD8T'))
cd4.ridge.raw = RidgePlot(Titr188 %>% subset(type!='remaining'), features='CD4', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD4 Raw Counts') + title.center +
  annotate("text", x = 3.6, y = 3.3, label = sprintf("LFC: %.2f", cd4t.cd4$Raw)) +
  annotate("text", x = 1, y = 5.2, label = sprintf("LFC: %.2f", cm.cd4$Raw)) + labs(x='CLR Value')
cd4.ridge.scPDA = RidgePlot(Titr188.scPDA %>% subset(type!='remaining'), features='CD4', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD4 scPDA Denoised')+ title.center +
  annotate("text", x = 4, y = 3.1, label = sprintf("LFC: %.2f", cd4t.cd4$scPDA)) + 
  annotate("text", x = 1, y = 5.1, label = sprintf("LFC: %.2f", cm.cd4$scPDA)) + labs(x='CLR Value')

cm.cd14 = LFC(raw=data.raw, denoise=data.pda, protein='CD14', pos=c('CM'))
cd14.ridge.raw = RidgePlot(Titr188 %>% subset(type!='remaining'), features='CD14', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD14 Raw Counts') + title.center +
  annotate("text", x = 0.6, y = 5.2, label = sprintf("LFC: %.2f", cm.cd14$Raw)) + labs(x='CLR Value')
cd14.ridge.scPDA = RidgePlot(Titr188.scPDA %>% subset(type!='remaining'), features='CD14', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD14 scPDA Denoised')+ title.center +
  annotate("text", x = 0.75, y = 5.1, label = sprintf("LFC: %.2f", cm.cd14$scPDA)) + labs(x='CLR Value')

p1 = grid.arrange(cut2, cut1,cd4.ridge.raw,cd4.ridge.scPDA,cd14.ridge.raw, cd14.ridge.scPDA, layout_matrix=mat.all)
##################################### Figure b #######################################################
cut1 = FeatureScatter(Titr188.scPDA %>% subset(type!='remaining'), feature1='CD8', feature2='CD19',jitter = TRUE)+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  cowplot::theme_nothing()+xlab('CD8')+ylab('CD19')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(face="plain", color="black", size=11),
        axis.title.y = element_text(face="plain", color="black", size=11,angle=90, vjust=0.5))
data.pda = Titr188.scPDA %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
x_density1 <- ggplot(data.pda, aes(x=CD8)) + geom_density(fill="steelblue", bw = 0.1) + 
  theme_void() + ggtitle('scPDA Denoised') + theme(plot.title = element_text(hjust = 0.5, vjust=-5, face='plain'))
y_density1 <- ggplot(data.pda, aes(x=CD19)) + geom_density(fill="orange", bw = 0.05) + theme_void() +coord_flip()
cut1 = grid.arrange(x_density1, cut1, y_density1, type_leg, layout_matrix=mat.cut)

cut2 = FeatureScatter(Titr188 %>% subset(type!='remaining'), feature1='CD8', feature2='CD19', jitter = TRUE)+NoLegend()+
  geom_vline(xintercept = 1.5, linetype="dashed")+
  geom_hline(yintercept = 1, linetype="dashed")+
  cowplot::theme_nothing()+xlab('CD8')+ylab('CD19')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.x = element_text(face="plain", color="black", size=11),
        axis.title.y = element_text(face="plain", color="black", size=11,angle=90, vjust=0.5))
data.raw = Titr188 %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
x_density2 <- ggplot(data.raw, aes(x=CD8)) + geom_density(fill="steelblue", bw = 0.1) + 
  theme_void() +ggtitle('Raw Counts') + theme(plot.title = element_text(hjust = 0.5, vjust=-5, face='plain'))
y_density2 <- ggplot(data.raw, aes(x=CD19)) + geom_density(fill="orange", bw = 0.05) + theme_void() +coord_flip()
cut2 = grid.arrange(x_density2, cut2, y_density2, type_leg, layout_matrix=mat.cut)

cd8t.cd8=LFC(raw=data.raw, denoise=data.pda, protein='CD8', pos=c('CD8T'))
cd8.ridge.raw = RidgePlot(Titr188 %>% subset(type!='remaining'), features='CD8', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD8 Raw Counts') + title.center +
  annotate("text", x = 3.7, y = 2.2, label = sprintf("LFC: %.2f", cd8t.cd8$Raw)) + labs(x='CLR Value')
cd8.ridge.scPDA = RidgePlot(Titr188.scPDA %>% subset(type!='remaining'), features='CD8', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD8 scPDA Denoised')+ title.center +
  annotate("text", x = 4, y = 2.5, label = sprintf("LFC: %.2f", cd8t.cd8$scPDA)) + labs(x='CLR Value')

B.cd19 = LFC(raw=data.raw, denoise=data.pda, protein='CD19', pos=c('B'))
cd19.ridge.raw = RidgePlot(Titr188 %>% subset(type!='remaining'), features='CD19', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD19 Raw Counts') + title.center +
  annotate("text", x = 2.2, y = 4.2, label = sprintf("LFC: %.2f", B.cd19$Raw)) + labs(x='CLR Value')
cd19.ridge.scPDA = RidgePlot(Titr188.scPDA %>% subset(type!='remaining'), features='CD19', slot='data') + 
  NoLegend() + theme(axis.title.y = element_blank()) + ggtitle('CD19 scPDA Denoised')+ title.center +
  annotate("text", x = 2.5, y = 4.5, label = sprintf("LFC: %.2f", B.cd19$scPDA)) + labs(x='CLR Value')
p2 = grid.arrange(cut2, cut1,cd8.ridge.raw,cd8.ridge.scPDA,cd19.ridge.raw,cd19.ridge.scPDA, layout_matrix=mat.all)

plots_list = list()
plots_list[['A']] = p1
plots_list[['B']] = p2
combined_plot <- grid.arrange(grobs = plots_list, ncol = 2)
ggsave("fig_reprod/Figure3_a-b.png", combined_plot, width = 24, height = 8, dpi = 300)

##################################################################################
#################################### Figure Log Fold Change ######################
LFC = function(df.list, protein, pos, neg=NA){
  if (all(is.na(neg))){
    types = c("NK","CD8T","CD4T","B","CM")
    neg = setdiff(types, pos)
  }
  lfc.helper = function(df){
    df.pos = df %>% filter(type %in% pos) %>% pull(protein) %>% mean()
    df.neg = df %>% filter(type %in% neg) %>% pull(protein) %>% mean()
    lfc = log(df.pos/df.neg, 2)
    return(lfc)
  }
  lfc = sapply(df.list, lfc.helper) %>% t() %>% as.data.frame(row.names = protein)
  return(lfc)
}
data.raw = Titr188 %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
data.pda = Titr188.scPDA %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
data.scar = Titr188.scAR %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
data.dsb = Titr188.DSB %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
data.gmm = Titr188.gmm %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
data.decont = Titr188.DecontPro %>% GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame() %>% mutate(type=Titr188.scPDA$type)
df.list = list('Raw'=data.raw, 'GMM'=data.gmm, 'scAR'=data.scar, 'DSB'=data.dsb, 'scPDA'=data.pda, 'DecontPro'=data.decont)

lfc.df = data.frame()
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD14', pos=c('CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD16', pos=c('NK')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD19', pos=c('B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD20', pos=c('B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD21', pos=c('B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD3', pos=c('CD4T', 'CD8T')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD32', pos=c('CM', 'B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD328', pos=c('CM', 'NK')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD33', pos=c('CM', 'NK')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD35', pos=c('CM', 'B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD36', pos=c('CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD4', pos=c('CD4T','CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD41', pos=c('CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD5', pos=c('CD4T', 'CD8T'), neg=c('NK','CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD56', pos=c('NK')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD7', pos=c('CD4T', 'CD8T','NK')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CD8', pos=c('CD8T')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='CLEC12A', pos=c('CM')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='HLA.DR', pos=c('CM', 'B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='IgD', pos=c('B')))
lfc.df = lfc.df %>% rbind(LFC(df.list, protein='IgM', pos=c('B')))

df_long.lfc <- tidyr::gather(lfc.df %>% rownames_to_column(var='Protein'), key = "method", value = "LFC", -Protein)
df_long.lfc$method <- factor(df_long.lfc$method, levels = c("Raw", "scAR", 'DSB', 'DecontPro', 'GMM', 'scPDA'))
method.name = levels(df_long.lfc$method)
# df_long.lfc <- df_long.lfc %>% filter(method %in% c('Raw', 'scPDA'))
colors <- setNames(rgb_alpha(c("grey", 'green4', 'purple4', "steelblue", "orange1", "red3"), 0.6), method.name)
p.lfc = ggplot(df_long.lfc, aes(x = Protein, y = LFC, color = method)) + geom_point(size = 3) + 
  theme_minimal() + title.center + box + 
  labs(title='Log Fold Change Comparison', x='Proteins', y='Log Fold Change') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.justification = 0.5, 
        legend.title = element_blank(), legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face='plain'))+
  scale_color_manual(values = colors)+coord_flip()+
  guides(color = guide_legend(nrow = 2))
ggsave('fig_reprod/Figure3_d.png', plot=p.lfc, width=6, height=9, dpi=300)
#################################### Figure UMAP #################################
##################################################################################
remain = Titr188.scPDA %>% subset(type=='remaining') %>% 
  GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame()
new.type = rep('Unclassified', nrow(remain)); names(new.type) = rownames(remain)
for (i in 1:nrow(remain)){
  cell = remain[i,]
  # CD4
  if (cell$CD3>1 & cell$CD4>2 & cell$CD8<2){
    new.type[i] = 'CD4T'
  }
  # CD8
  if (cell$CD3>1 & cell$CD4<2 & cell$CD8>2){
    new.type[i] = 'CD8T'
  }
  # B
  if (cell$CD3<1 & cell$CD19>1){
    new.type[i] = 'B'
  }
  # CM
  if (cell$CD3<1 & cell$CD19<1 & cell$CD14>0.3 & cell$CD16<1.5){
    new.type[i] = 'CM'
  }
  # NK
  if (cell$CD3<1 & cell$CD19<1 & cell$CD14<0.3 & cell$CD56>0.5){
    new.type[i] = 'NK'
  }
}
new.type = as.data.frame(new.type)
new.type$new.type = factor(new.type$new.type, levels=c("NK","CD8T","CD4T","B","CM","Unclassified"))

plt.cells = Titr188.scPDA; margin = theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
temp = plt.cells
temp$type[temp$type=='remaining'] = 'Unclassified'
Idents(temp)='type'
Idents(temp) <- factor(Idents(temp), levels = c('NK', "CD8T", "CD4T", "B", "CM", "Unclassified"))
my.col = hue_pal()(6); my.col[6]='grey'
p1 = DimPlot(temp, reduction = 'rnaumap', pt.size=0.5, cols = my.col) + theme_void() + 
  theme.border + ggtitle('3358 All Cells') + title.center + margin + 
  theme(legend.position = c(1, 1), legend.justification = c(1.1, 1),
        legend.key.height = unit(0.5, "line"))
p2 = DimPlot(plt.cells %>% subset(type == 'remaining'), reduction = 'rnaumap', cols='grey', pt.size=0.5) + theme_void() + 
  NoLegend() + theme.border + ggtitle(sprintf('%d Unclassified Cells', nrow(new.type))) + title.center + margin

remain.subset = plt.cells %>% subset(type == 'remaining') %>% 
  AddMetaData(metadata = new.type) %>% subset(new.type != 'Unclassified')
p3 = DimPlot(remain.subset, reduction = 'rnaumap',group.by = 'new.type', label=TRUE) + theme_void() + 
  theme.border + title.center + margin +
  ggtitle(sprintf('%d Cells Classified After Denoising', sum(new.type$new.type != 'Unclassified'))) +
  theme(legend.position = c(1, 1), legend.justification = c(1.1, 1),
        legend.key.height = unit(0.5, "line"))
layout_mat <- matrix(c(1, 2, 3), ncol = 3, byrow = TRUE)
# Arrange your plots with labels
combUmap = grid.arrange(p1, p2, p3, layout_matrix = layout_mat,
                        bottom = textGrob("UMAP 1"), left = textGrob("UMAP 2", rot = 90))
ggsave("fig_reprod/Figure3_c.png", combUmap, width = 12, height = 3, dpi = 500)
#################################### DotPlotmap #####################################
##################################################################################
mydot = function(dataset, title){
  reshape <- dataset %>% 
    gather(key = "Protein", value = "Expression", -type) %>% 
    group_by(type, Protein) %>%
    summarize(`Median Expression` = log1p(median(Expression)), 
              `Percent Expressed` = mean(Expression > 0) * 100, 
              .groups = 'drop')
  E.min = min(reshape$`Median Expression`); E.max=max(reshape$`Median Expression`)
  sparsity = sum(reshape$`Median Expression` <= 0)/nrow(reshape)
  reshape <- reshape %>% filter(`Percent Expressed` >= 10)
  ggplot(reshape, aes(x=type, y=Protein, size=`Percent Expressed`, color=`Median Expression`)) +
    geom_point() + theme_minimal()+
    guides(size = guide_legend(override.aes = list(color = "gray")), 
           color = guide_colorbar()) +
    labs(title=sprintf('%s (%.2f)', title, sparsity), x=NULL, y=NULL) + 
    scale_color_gradientn(
      colours = c("grey", "red4", "red2"),
      values = rescale(c(0, 0.15, 0.75, 1)),
      breaks = c(E.min, mean(c(E.min, E.max)),E.max),
      labels = c("Low", "Mid","High")) + 
    title.center+box+RotatedAxis()
}

# Raw
raw.dot=raw.counts %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.raw=mydot(raw.dot, title='Raw') + NoLegend()
# DSB
dsb.dot = expm1(dsb.counts) %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.dsb=mydot(dsb.dot, title='DSB')+NoLegend()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# DecontPro
decont.dot = decont.counts %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.decont=mydot(decont.dot, title='DecontPro')+NoLegend()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# scAR
scar.dot = scar.counts %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.scar=mydot(scar.dot, title='scAR')+NoLegend()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# GMM
gmm.dot=gmm.counts %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.gmm=mydot(gmm.dot, title='GMM')+ NoLegend()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
# scPDA
scpda.dot = scPDA.counts %>% t %>% as.data.frame() %>% mutate(type=Titr188$type) %>% filter(type!='remaining')
dotmap.scpda=mydot(scpda.dot, title='scPDA')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

heat_all = dotmap.raw+dotmap.dsb+dotmap.decont+dotmap.scar+dotmap.gmm+dotmap.scpda+plot_layout(ncol = 6)
ggsave("fig_reprod/Figure3_e.png", heat_all, width = 12, height = 9, dpi = 300)

##################### Combine Plots #####################
row1 <- image_read("fig_reprod/Figure3_a-b.png")
row2 <- image_read("fig_reprod/Figure3_c.png")
img3 <- image_read("fig_reprod/Figure3_d.png") %>% image_background(color='white')
img4 <- image_read("fig_reprod/Figure3_e.png")

target_width <- image_info(row1)$width
row2 <- image_resize(row2, paste0(target_width, "x"))
ratio_d = image_info(img3)$width/(image_info(img3)$width+image_info(img4)$width)
img3 = image_resize(img3, paste0(target_width*ratio_d, "x"))
img4 = image_resize(img4, paste0(target_width*(1-ratio_d), "x"))
row3 <- image_append(c(img3, img4))

final_image <- image_append(c(row1, row2, row3), stack = TRUE)
image_write(final_image, path = "fig_reprod/Figure3.jpg", format="jpeg", quality=60)

file.remove("fig_reprod/Figure3_a-b.png",
            "fig_reprod/Figure3_c.png",
            "fig_reprod/Figure3_d.png",
            "fig_reprod/Figure3_e.png")