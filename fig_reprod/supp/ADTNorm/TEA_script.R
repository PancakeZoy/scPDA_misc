require(ADTnorm)
require(dplyr)

back = function(x) (sinh(x)-1)*5
counts = readRDS('data/TEA_counts.rds') %>% t()
features = data.frame(sample=rep('TEA', nrow(counts)))
all_markers = colnames(counts)

error_prots = c("HLA-DR","IgG1-K-Isotype-Control","TCR-Va24-Ja18","TCR-Va7.2","TCR-a/b","TCR-g/d")
valid_prots = sort(setdiff(all_markers, error_prots))
cell_x_adt_norm <- ADTnorm(
  cell_x_adt = counts, 
  cell_x_feature = features, 
  save_outpath = '~/data', 
  study_name = 'TEA',
  marker_to_process = valid_prots,
  bimodal_marker = NULL, 
  trimodal_marker = NULL, 
  positive_peak = NULL,
  brewer_palettes = "Dark2",
  save_fig = FALSE,
  target_landmark_location = c(arcsinh_transform(0), 3))
cell_x_adt_norm = cell_x_adt_norm %>% select(all_of(valid_prots))
after_counts = back(cell_x_adt_norm)

saveRDS(after_counts, 'data/TEA_denoised.rds')
