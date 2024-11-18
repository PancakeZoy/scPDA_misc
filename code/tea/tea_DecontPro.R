require(here)
setwd(here())
require(dplyr)
require(decontX)
require(Seurat)
set.seed(1)
raw.counts = read.csv('results/tea/tea_raw.csv', row.names = 1)
adt_seurat <- CreateSeuratObject(raw.counts, assay = "Protein")
adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "Protein") %>%
  RunPCA(assay = "Protein", features = rownames(adt_seurat), npcs = 20,
         reduction.name = "pca_prot", nfeatures.print=3) %>%
  FindNeighbors(dims = 1:20, assay = "Protein", reduction = "pca_prot") %>% FindClusters()
clusters <- as.integer(Idents(adt_seurat))
counts <- as.matrix(raw.counts)
start_time <- Sys.time()
out <- decontPro(counts,clusters)
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)
decontaminated_counts <- out$decontaminated_counts
write.csv(decontaminated_counts, file='results/tea/tea_DecontPro.csv', row.names = TRUE)