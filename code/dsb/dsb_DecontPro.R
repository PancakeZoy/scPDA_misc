require(here)
setwd(here())
require(dplyr)
require(decontX)
require(Seurat)
set.seed(1)
dsb = readRDS('data/dsb.rds')
raw.counts = dsb %>% GetAssayData(slot='counts') %>% as.matrix()
adt_seurat <- CreateSeuratObject(raw.counts, assay = "Protein")
adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "Protein") %>%
  RunPCA(assay = "Protein", features = rownames(adt_seurat), npcs = 10,
         reduction.name = "pca_prot") %>%
  FindNeighbors(dims = 1:10, assay = "Protein", reduction = "pca_prot") %>%
  FindClusters(resolution = 0.2)
clusters <- as.integer(Idents(adt_seurat))
counts <- as.matrix(raw.counts)
start_time <- Sys.time()
out <- decontPro(counts,clusters)
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)
decontaminated_counts <- out$decontaminated_counts
write.csv(decontaminated_counts, file='results/dsb/dsb_DecontPro.csv', row.names = TRUE)
