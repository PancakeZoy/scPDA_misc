require(here)
setwd(here())
source('code/tools.R')
gmm_mu1_dir = 'results/titr188/titr188_GMM_mu1.csv'
meta_dir = 'results/titr188/titr188_meta.csv'
raw_dir = 'results/titr188/titr188_raw.csv'
dsb_dir = 'results/titr188/titr188_DSB.csv'
gmm_dir = 'results/titr188/titr188_GMM.csv'
decontpro_dir = 'results/titr188/titr188_Decontpro.csv'
scar_dir = 'results/titr188/titr188_scAR.csv'
scpda_dir = 'results/titr188/titr188_scPDA.csv'

# Save Data Ready for Comparison
Titr188 = readRDS('data/titr188.rds')
raw.counts = Titr188 %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)

set.seed(1)
bg_noise_list = apply(raw.counts, 1, function(x) {Mclust(x, G = 2:3, warn = FALSE, verbose = FALSE)})
mu1 = sapply(bg_noise_list, function(x) sort(x$parameters$mean)[1] %>% unname) %>% as.data.frame()
colnames(mu1)='noise'
write.csv(raw.counts, file=raw_dir, row.names = TRUE)
write.csv(mu1, file=gmm_mu1_dir, row.names = TRUE)
write.csv(Titr188@meta.data, file=meta_dir, row.names = TRUE)

####################### Run scPDA and scAR #######################
folder.path = here('code','titr188')
sh.path = file.path(folder.path, 'titr188_py.sh')
sh.script <- c("source ~/.zshrc",
               "conda activate denoise",
               paste0("cd ", folder.path),
               "python titr188_scPDA.py",
               "python titr188_scAR.py",
               "conda deactivate")
writeLines(sh.script, sh.path)
system(paste0('chmod +x ', sh.path))
system(sh.path)
file.remove(sh.path)

pi = h5read('results/titr188/titr188_scPDA.h5', 'pi') %>% as.matrix()
colnames(pi)=colnames(raw.counts); rownames(pi)=rownames(raw.counts)
scPDA.counts = raw.counts * (1-pi)
write.csv(scPDA.counts, file=scpda_dir, row.names = TRUE)

scar.counts = h5read('results/titr188/titr188_scAR.h5', 'count') %>% as.matrix()
rownames(scar.counts) = prot.names; colnames(scar.counts) = cell.names
write.csv(scar.counts, file=scar_dir, row.names = TRUE)
####################### DSB #######################
set.seed(1)
dsb.counts =  ModelNegativeADTnorm(cell_protein_matrix = raw.counts,
                                   denoise.counts = TRUE,
                                   use.isotype.control = FALSE,
                                   define.pseudocount = 1)
write.csv(dsb.counts, dsb_dir, row.names = TRUE)

####################### GMM #######################
set.seed(1)
bg_noise = apply(log1p(raw.counts), 1, function(x) {Mclust(x, G = 2:3, warn = FALSE, verbose = FALSE)})
gmm.pi = sapply(bg_noise, function(x) x$z[,1] %>% unname) %>% t
gmm.counts = raw.counts * (1-gmm.pi)
write.csv(gmm.counts, gmm_dir, row.names = TRUE)

####################### DecontPro #######################
require(decontX)
set.seed(1)
adt_seurat <- CreateSeuratObject(raw.counts, assay = "Protein")
adt_seurat <- NormalizeData(adt_seurat, normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "Protein") %>%
  RunPCA(assay = "Protein", features = rownames(adt_seurat), npcs = 20,
         reduction.name = "pca_prot", nfeatures.print=3) %>%
  FindNeighbors(dims = 1:20, assay = "Protein", reduction = "pca_prot") %>% FindClusters()
clusters <- as.integer(Idents(adt_seurat))
counts <- as.matrix(raw.counts)
start_time <- Sys.time()
output <- decontPro(counts,clusters)
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)
decont.counts <- output$decontaminated_counts
write.csv(decont.counts, file=decontpro_dir, row.names = TRUE)