require(ADTnorm)
require(dplyr)

back = function(x) (sinh(x)-1)*5
counts = readRDS('data/DSB_counts.rds')
colnames(counts)[grepl("sotype", colnames(counts))] = paste0("Isotype", 1:4)
colnames(counts) = gsub(" ", "", colnames(counts))
meta = read.csv('data/DSB_meta.csv', row.names = 1)
features = data.frame(sample=rep('DSB', nrow(counts)))
features$type = meta$Type1
all_markers = counts %>% colnames

######################## Neg Control Cells ########################
neg_idx = meta$Type1=='Neg Control'
error_prots = c('CD123', 'CD71', "CD16", "CD278", "CD314","HLA-ABC","HLA-DR")
valid_prots = sort(setdiff(all_markers, error_prots))
cell_x_adt_norm <- ADTnorm(
  cell_x_adt = counts, 
  cell_x_feature = features, 
  save_outpath = '~/data', 
  study_name = 'DSB',
  marker_to_process = valid_prots,
  bimodal_marker = NULL, 
  trimodal_marker = NULL, 
  positive_peak = NULL,
  brewer_palettes = "Dark2",
  save_fig = FALSE,
  target_landmark_location = c(arcsinh_transform(0), 3))
before = counts[neg_idx,] %>% select(all_of(valid_prots))
after = cell_x_adt_norm[neg_idx,] %>% select(all_of(valid_prots)) %>% back()
nrc.unstain = 1-colMeans(after)/colMeans(before)
mean(nrc.unstain); median(nrc.unstain)
saveRDS(nrc.unstain, 'data/dsb_unstain_nrc.rds')


######################## Non-expressive Proteins ########################
nonexpr_markers = c('CD117','CD137','CD138','CD206','CD223','CD273','CD80', 
                    "Isotype1","Isotype2","Isotype3","Isotype4")
cell_x_adt_norm <- ADTnorm(
  cell_x_adt = counts, 
  cell_x_feature = features, 
  save_outpath = '~/data', 
  study_name = 'DSB',
  marker_to_process = nonexpr_markers,
  bimodal_marker = NULL, 
  trimodal_marker = NULL, 
  positive_peak = NULL,
  brewer_palettes = "Dark2",
  save_fig = FALSE,
  target_landmark_location = c(arcsinh_transform(0), 3))
before = counts %>% select(all_of(nonexpr_markers))
after = cell_x_adt_norm %>% select(all_of(nonexpr_markers)) %>% back()
nrc.nonexpr = 1-colMeans(after)/colMeans(before)
names(nrc.nonexpr) = nonexpr_markers
mean(nrc.nonexpr); median(nrc.nonexpr)
saveRDS(nrc.nonexpr, 'data/dsb_nonexpr_nrc.rds')

