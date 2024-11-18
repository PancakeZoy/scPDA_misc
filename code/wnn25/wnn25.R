require(here)
setwd(here())
source('code/tools.R')
count_dir = 'results/wnn25/wnn25_raw.csv'
gmm_mu1_dir = 'results/wnn25/wnn25_GMM_mu1.csv'
meta_dir = 'results/wnn25/wnn25_meta.csv'
dsb_dir = 'results/wnn25/wnn25_DSB.rds'
gmm_dir = 'results/wnn25/wnn25_GMM.rds'

wnn25 <- readRDS('data/wnn25.rds')
DefaultAssay(wnn25) <- 'ADT'
raw.counts = wnn25 %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)

set.seed(1)
bg_noise_list = apply(raw.counts, 1, function(x) {Mclust(x, G = 1:3, modelNames="V", warn = FALSE, verbose = FALSE)})
mu1 = sapply(bg_noise_list, function(x) sort(x$parameters$mean)[1] %>% unname) %>% as.data.frame()
colnames(mu1)='noise'
write.csv(raw.counts, file=count_dir, row.names = TRUE)
write.csv(mu1, file=gmm_mu1_dir, row.names = TRUE)
write.csv(wnn25@meta.data, file=meta_dir, row.names = TRUE)

####################### Run scPDA and scAR #######################
folder.path = here('code','wnn25')
sh.path = file.path(folder.path, 'wnn25_py.sh')
sh.script <- c("source ~/.zshrc",
               "conda activate denoise",
               paste0("cd ", folder.path),
               "python wnn25_scPDA.py",
               "python wnn25_scAR.py",
               "conda deactivate")
writeLines(sh.script, sh.path)
system(paste0('chmod +x ', sh.path))
system(sh.path)
file.remove(sh.path)

####################### DSB #######################
set.seed(1)
dsb.counts =  ModelNegativeADTnorm(cell_protein_matrix = raw.counts, 
                                   denoise.counts = TRUE, 
                                   use.isotype.control = FALSE,
                                   define.pseudocount = 0.5)
saveRDS(dsb.counts, dsb_dir)

####################### GMM #######################
set.seed(1)
bg_noise = apply(log1p(raw.counts), 1, function(x) {Mclust(x, G = 2, warn = FALSE, verbose = FALSE)})
gmm.pi = sapply(bg_noise, function(x) x$z[,1] %>% unname) %>% t
gmm.counts = round(raw.counts * (1-gmm.pi))
saveRDS(gmm.counts, gmm_dir)

####################### DecontPro #######################
source('code/wnn25/wnn25_DecontPro.R')