library(here)
setwd(here())
source('code/tools.R')
plot_protein <- function(prot, pos.name=NULL, pos=NULL, counts_list, types, use_negpro=FALSE, plot_title=NULL) {
  if (is.null(plot_title)) {
    plot_title <- prot
  }
  
  if (!use_negpro) {
    # For standard proteins
    p.raw <- overlap.f1(counts=counts_list$raw.counts, prot=prot, pos=pos, pos.name=pos.name, title='Raw', types=types) +
      labs(x=NULL, y=NULL) +
      theme(legend.position = "top", legend.justification = "center", legend.box = "horizontal",
            legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8)) +
      guides(fill=guide_legend(title=NULL))
    
    p.scPDA <- overlap.f1(counts=counts_list$scPDA.counts, prot=prot, pos=pos, pos.name=pos.name, title='scPDA Denoised', types=types) +
      NoLegend() + labs(x=NULL, y=NULL)
    
    # Extract legend from p.raw
    g <- ggplotGrob(p.raw)
    legend_grob <- g$grobs[[which(g$layout$name == "guide-box")]]
    
    p.raw <- p.raw + NoLegend()
    p.all <- grid.arrange(p.raw, p.scPDA, ncol=2)
    p.all <- grid.arrange(legend_grob, p.all, heights=c(0.05,1), top=textGrob(plot_title, gp=gpar(fontsize=16, fontface="bold")),
                          left='Density', bottom='Counts')
  } else {
    # For negative control proteins
    p.raw <- overlap.negpro(counts_list$raw.counts, prot=prot, title='Raw') +
      labs(x=NULL, y=NULL) + NoLegend() +
      theme(legend.position = "top", legend.justification = "center", legend.box = "horizontal",
            legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 8)) +
      guides(fill=guide_legend(title=NULL))
    
    p.scPDA <- overlap.negpro(counts_list$scPDA.counts, prot=prot, title='scPDA Denoised', col='orange') +
      NoLegend() + labs(x=NULL, y=NULL)
    
    # Since there is no legend, we skip legend extraction
    p.all <- grid.arrange(p.raw, p.scPDA, ncol=2,
                          top=textGrob(plot_title, gp=gpar(fontsize=16, fontface="bold")),
                          left='Density', bottom='Counts')
  }
  
  # Return the combined plot
  return(p.all)
}

#################### Load data ############################
dsb = readRDS('data/dsb.rds')
prot.names = rownames(dsb)
raw.counts = dsb %>% GetAssayData(slot='counts') %>% as.matrix()
cell.counts = dsb %>% subset(Type1 != 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
unstain.counts = dsb %>% subset(Type1 == 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
types = dsb$Type1

pi = h5read('results/dsb/dsb_scPDA.h5', 'pi') %>% as.matrix()
colnames(pi) = colnames(raw.counts)
rownames(pi) = rownames(raw.counts)
scPDA.counts = round((1 - pi) * raw.counts)

counts_list <- list(raw.counts = raw.counts, scPDA.counts = scPDA.counts)
############################################################
plots_list <- list()

prot='CD8'
pos.name= 'CD8 T Cells'
pos= c('CD8+ Memory T','CD8+ Naive T (some DNT)')
plots_list[['CD8']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

prot='CD4'
pos.name= 'CD4 T Cells'
pos= c('CD4+ Memory T','CD4+ Naive T')
plots_list[['CD4']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

prot='CD3'
pos.name= 'T Cells'
pos= c('CD4+ Memory T','CD4+ Naive T','CD8+ Memory T','CD8+ Naive T (some DNT)')
plots_list[['CD3']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

prot='CD19'
pos.name= 'B Cells'
pos= c('B cells')
plots_list[['CD19']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

prot='CD20'
pos.name= 'B Cells'
pos= c('B cells')
plots_list[['CD20']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

prot='CD14'
pos.name= 'Classical Monocytes and mDC'
pos= c('Classical Monocytes and mDC')
plots_list[['CD14']] <- plot_protein(prot=prot, pos.name=pos.name, pos=pos, counts_list=counts_list, types=types)

# Negative control proteins
prot <- 'CD117'
plot_title <- 'CD117'
plots_list[['CD117']] <- plot_protein(prot=prot, counts_list=counts_list, types=types, use_negpro=TRUE, plot_title=plot_title)

prot <- 'CD138'
plot_title <- 'CD138'
plots_list[['CD138']] <- plot_protein(prot=prot, counts_list=counts_list, types=types, use_negpro=TRUE, plot_title=plot_title)

prot <- 'MouseIgG1kappaisotype'
plot_title <- 'Isotype1'
plots_list[['Isotype1']] <- plot_protein(prot=prot, counts_list=counts_list, types=types, use_negpro=TRUE, plot_title=plot_title)

combined_plot <- grid.arrange( grobs = plots_list, ncol = 3)
ggsave("fig_reprod/Figure1_a-i.png", combined_plot, width = 15, height = 9, dpi = 200)