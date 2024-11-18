require(here)
setwd(here())
source('code/tools.R')
abp_dir = 'results/tea/tea_EmptyProfile.csv'
count_dir = 'results/tea/tea_raw.csv'
gmm_mu1_dir = 'results/tea/tea_GMM_mu1.csv'
meta_dir = 'results/tea/tea_meta.csv'
dsb_dir = 'results/tea/tea_dsb.rds'

tea <- readRDS('data/teaseq.rds')
empty = tea$empty; tea=tea$teaseq %>% subset(celltype != 'Unknown'); Idents(tea) = 'celltype'
raw.counts = tea %>% GetAssayData(slot='counts') %>% as.matrix()
prot.names=rownames(raw.counts); cell.names=colnames(raw.counts)
noise = rowMeans(empty); type_vec=tea$celltype

set.seed(1)
bg_noise_list = apply(raw.counts, 1, function(x) {Mclust(x, G = 1:3, modelNames="V", warn = FALSE, verbose = FALSE)})
mu1 = sapply(bg_noise_list, function(x) sort(x$parameters$mean)[1] %>% unname) %>% as.data.frame()
colnames(mu1)='noise'; amb.profile = as.data.frame(noise %>% proportions())
write.csv(raw.counts, file=count_dir, row.names = TRUE)
write.csv(amb.profile, file=abp_dir, row.names = TRUE)
write.csv(mu1, file=gmm_mu1_dir, row.names = TRUE)
write.csv(tea@meta.data, file=meta_dir, row.names = TRUE)

####################### Run scPDA and scAR #######################
folder.path = here('code','tea')
sh.path = file.path(folder.path, 'tea_py.sh')
sh.script <- c("source ~/.zshrc",
               "conda activate denoise",
               paste0("cd ", folder.path),
               "python tea_scPDA.py",
               "python tea_scAR.py",
               "conda deactivate")
writeLines(sh.script, sh.path)
system(paste0('chmod +x ', sh.path))
system(sh.path)
file.remove(sh.path)

################################# DSB #################################
set.seed(1)
isotypes=rownames(raw.counts)[grep('type', rownames(raw.counts))]
dsb.counts =  DSBNormalizeProtein(cell_protein_matrix = raw.counts, empty_drop_matrix = empty,
                                  denoise.counts = TRUE, use.isotype.control = TRUE,
                                  isotype.control.name.vec = isotypes, define.pseudocount = 1)
dsb.counts = expm1(dsb.counts)
saveRDS(dsb.counts, dsb_dir)

############################### DecontPro #############################
source('code/tea/tea_DecontPro.R')