require(here)
setwd(here())
source('code/tools.R')
count_dir = 'results/dsb/dsb_raw.csv'
gmm_mu1_dir = 'results/dsb/dsb_GMM_mu1.csv'
abp_dir = 'results/dsb/dsb_EmptyProfile.csv'
dsb_dir = 'results/dsb/dsb_DSB.rds'
gmm_dir = 'results/dsb/dsb_GMM.rds'

# Save Data Ready for Comparison
dsb = readRDS('data/dsb.rds')
raw.counts = dsb %>% GetAssayData(slot='counts') %>% as.matrix()
cell.counts = dsb %>% subset(Type1 != 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
unstain.counts = dsb %>% subset(Type1 == 'Neg Control') %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); cell.names=colnames(cell.counts); unstain.names=colnames(unstain.counts)
true.noise = rowMeans(unstain.counts)

set.seed(1)
bg_noise_list = apply(raw.counts, 1, function(x) {Mclust(x, G = 1:3, modelNames="V", warn = FALSE, verbose = FALSE)})
mu1 = sapply(bg_noise_list, function(x) sort(x$parameters$mean)[1] %>% unname) %>% as.data.frame()
colnames(mu1)='noise'
amb.profile = readRDS('data/dsb_empty.rds') %>% 
  list.cbind %>% as.matrix() %>% rowMeans() %>% proportions()

write.csv(raw.counts, file=count_dir, row.names = TRUE)
write.csv(mu1, file=gmm_mu1_dir, row.names = TRUE)
write.csv(amb.profile, file=abp_dir, row.names = TRUE)
####################### Run scPDA and scAR #######################
folder.path = here('code','dsb')
sh.path = file.path(folder.path, 'DSB_py.sh')
sh.script <- c("source ~/.zshrc",
               "conda activate denoise",
               paste0("cd ", folder.path),
               "python DSB_scPDA.py",
               "python DSB_scAR.py",
               "conda deactivate")
writeLines(sh.script, sh.path)
system(paste0('chmod +x ', sh.path))
system(sh.path)
file.remove(sh.path)

####################### DSB #######################
set.seed(1)
hash.empty = readRDS('data/dsb_empty.rds') %>% list.cbind %>% as.matrix()
rownames(hash.empty)=prot.names
isotypes = prot.names[grep('type', prot.names)]
dsb.counts = DSBNormalizeProtein(cell_protein_matrix = raw.counts,
                                 empty_drop_matrix = hash.empty,
                                 denoise.counts = TRUE, 
                                 use.isotype.control = TRUE, 
                                 isotype.control.name.vec = isotypes,  
                                 define.pseudocount = TRUE,
                                 pseudocount.use = 1)
saveRDS(dsb.counts, dsb_dir)

####################### GMM #######################
set.seed(1)
bg_noise = apply(log1p(raw.counts), 1, function(x) {Mclust(x, G = 2, warn = FALSE, verbose = FALSE)})
gmm.pi = sapply(bg_noise, function(x) x$z[,1] %>% unname) %>% t
gmm.counts = round(raw.counts * (1-gmm.pi))
saveRDS(gmm.counts, gmm_dir)

####################### DecontPro #######################
source('code/dsb/dsb_DecontPro.R')