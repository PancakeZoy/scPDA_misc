require(dplyr)
require(Seurat)
require(ggplot2)
require(patchwork)
require(cowplot)

Titr188 = readRDS('data/titr188.rds')
Titr188@reductions = list(); Titr188@graphs = list(); Idents(Titr188) = 'type'
DefaultAssay(Titr188) = 'RNA'
Titr188 <- NormalizeData(Titr188) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca_rna')
DefaultAssay(Titr188) = 'Protein'

scpda_dir = 'results/titr188/titr188_scPDA.csv'
scPDA.counts = read.csv(scpda_dir, row.names = 1, check.names = FALSE) %>% as.matrix()
Titr188.scPDA = SetAssayData(Titr188, slot='counts', new.data=scPDA.counts)

pipeline = function(obj, res=0.1){
  obj = obj %>% NormalizeData(normalization.method = 'CLR', margin=2) %>% 
    ScaleData() %>% RunPCA(assay="Protein", reduction.name = 'pca_protein') %>% 
    FindMultiModalNeighbors(
      reduction.list = list("pca_rna", "pca_protein"), 
      dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight") %>% 
    RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap") %>% 
    FindClusters(graph.name = "wsnn", algorithm = 3, resolution = res)
  
  obj = obj %>% RunUMAP(reduction='pca_protein', dims=1:20, reduction.name='umap_protein')
}
res=0.4; wnn_label = paste0('wsnn_res.', res)
Titr188 = pipeline(Titr188, res)
Titr188.scPDA = pipeline(Titr188.scPDA, res)


my_theme = theme(
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, size = 1),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank(),
  plot.title = element_text(size = 14)
)

Titr188 = Titr188 %>% subset(type!='remaining')
Titr188.scPDA = Titr188.scPDA %>% subset(type!='remaining')
reduction = 'wnn.umap'
p1 = DimPlot(Titr188, reduction=reduction, group.by="type",label=TRUE)+
  ggtitle('RNA + Raw Protein (Cell Type) ')+my_theme+NoLegend()
p2 = DimPlot(Titr188.scPDA, reduction=reduction, group.by="type",label=TRUE)+
  ggtitle('RNA + Denoised Protein (Cell Type)')+my_theme+NoLegend()
p3 = DimPlot(Titr188, reduction=reduction, group.by=wnn_label,label=TRUE)+
  ggtitle('RNA + Raw Protein (Cluster)')+my_theme+NoLegend()
p4 = DimPlot(Titr188.scPDA, reduction=reduction, group.by=wnn_label,label=TRUE)+
  ggtitle('RNA + Denoised Protein (Cluster)')+my_theme+NoLegend()
fig.a = (p1+p2)/(p3+p4)


prot = 'CD45RA'
p5 = FeaturePlot(Titr188, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Raw)'))+my_theme
p6 = FeaturePlot(Titr188.scPDA, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Denoised)'))+my_theme
prot = 'CD45RO'
p7 = FeaturePlot(Titr188, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Raw)'))+my_theme
p8 = FeaturePlot(Titr188.scPDA, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Denoised)'))+my_theme
prot = 'CD14'
p9 = FeaturePlot(Titr188, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Raw)'))+my_theme
p10 = FeaturePlot(Titr188.scPDA, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Denoised)'))+my_theme
fig.b = (p5|p7|p9)/(p6|p8|p10)


prot = 'CD31'
p11 = FeaturePlot(Titr188.scPDA, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Denoised)'))+my_theme
prot = 'CD38'
p12 = FeaturePlot(Titr188.scPDA, reduction=reduction, features=prot)+
  NoLegend()+ggtitle(paste0(prot, ' (Denoised)'))+my_theme
fig.c = (p11/p12)


# row1
margin1 = 15; margin2 = 20
p1  <- p1  + theme(plot.margin = margin(margin2,  margin1,   margin1,   margin1))
p2  <- p2  + theme(plot.margin = margin(margin2, margin2,   margin1,   margin1))
p11 <- p11 + theme(plot.margin = margin(margin2,  margin1,   margin1,  margin2))
# row2
p3  <- p3  + theme(plot.margin = margin(margin1,  margin1,  margin2,   margin1))
p4  <- p4  + theme(plot.margin = margin(margin1, margin2,  margin2,   margin1))
p12 <- p12 + theme(plot.margin = margin(margin1,  margin1,  margin2,  margin2))
# row3
p5  <- p5  + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))
p7  <- p7  + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))
p9  <- p9  + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))
# row4
p6  <- p6  + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))
p8  <- p8  + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))
p10 <- p10 + theme(plot.margin = margin(margin2, margin1,   margin1,   margin1))

# —— 2) 把 12 个小图按行拼成 4×3 —— 
labels_vec <- c(
  "a","b","k",
  "c","d","l",
  "e","f","g",
  "h","i","j"
)
grid_plot <- plot_grid(
  p1,  p2,  p11,
  p3,  p4,  p12,
  p5,  p7,   p9,
  p6,  p8,  p10,
  ncol           = 3,
  byrow          = TRUE,
  labels         = labels_vec,
  label_size     = 20,
  label_fontface = "bold",
  label_x        = 0,
  label_y        = 1.02
)

ggsave(
  filename = "fig_reprod/Figure5.png",
  plot     = grid_plot,
  width    = 10.5,
  height   = 14,
  dpi      = 200
)
