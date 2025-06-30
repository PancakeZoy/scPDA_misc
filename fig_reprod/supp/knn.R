require(Seurat)
require(dplyr)
require(patchwork)

####################################################################################
# Load scPDA
Titr188 = readRDS('Data/Titr188/Titr188.rds')
Idents(Titr188) = 'type'
scpda_dir = 'Data/Titr188/Titr188_scpda.csv'
scPDA.counts = read.csv(scpda_dir, row.names = 1, check.names = FALSE) %>% as.matrix() %>% round()
Titr188.scPDA = SetAssayData(Titr188, slot='counts', new.data=scPDA.counts) %>% 
  NormalizeData(normalization.method = 'CLR', margin=2)

#################################### Figure UMAP #################################
##################################################################################
remain = Titr188.scPDA %>% subset(type=='remaining') %>% 
  GetAssayData(slot='data') %>% as.matrix %>% t %>% as.data.frame()
new.type = rep('Unclassified', nrow(remain)); names(new.type) = rownames(remain)
for (i in 1:nrow(remain)){
  cell = remain[i,]
  # CD4
  if (cell$CD3>1 & cell$CD4>2 & cell$CD8<2){
    new.type[i] = 'CD4T'
  }
  # CD8
  if (cell$CD3>1 & cell$CD4<2 & cell$CD8>2){
    new.type[i] = 'CD8T'
  }
  # B
  if (cell$CD3<1 & cell$CD19>1){
    new.type[i] = 'B'
  }
  # CM
  if (cell$CD3<1 & cell$CD19<1 & cell$CD14>0.3 & cell$CD16<1.5){
    new.type[i] = 'CM'
  }
  # NK
  if (cell$CD3<1 & cell$CD19<1 & cell$CD14<0.3 & cell$CD56>0.5){
    new.type[i] = 'NK'
  }
}

df = data.frame(Titr188.scPDA@reductions[["rnaumap"]]@cell.embeddings) %>% 
  mutate(old_type=Titr188.scPDA$type, new_type=new.type[colnames(Titr188)])

#### Tune a KNN on classified cells ####
require(caret)
set.seed(42)
tune_df = df %>% filter(old_type != 'remaining')
tune_df$old_type <- as.factor(tune_df$old_type)
X <- tune_df[, c("UMAP_1", "UMAP_2")]
Y <- tune_df$old_type

ctrl <- trainControl(method = "cv", number = 10)
knn_fit <- train(
  x = X,
  y = Y,
  method = "knn",
  trControl = ctrl,
  tuneGrid = expand.grid(k = 1:15)
)
print(knn_fit)

pred_df = df %>% filter(old_type == 'remaining') %>% 
  filter(new_type != 'Unclassified') # only focus on those originally unclassifed cells but now can be classified because of scPDA
pred_X = pred_df[, c("UMAP_1", "UMAP_2")]
Y = pred_df$new_type
pred <- predict(knn_fit, newdata = pred_X) %>% as.character()
acc = mean(pred == Y)

conf_mat <- confusionMatrix(as.factor(pred), as.factor(Y))
per_class_acc <- conf_mat$byClass[, "Sensitivity"]
per_class_acc

cm_df <- as.data.frame(conf_mat$table)
colnames(cm_df) <- c("Predicted","Reference","Count")
p_conf = ggplot(cm_df, aes(x=Reference, y=Predicted, fill=Count)) +
  geom_tile(color="white") +
  geom_text(aes(label=Count), size=4) +
  scale_fill_gradient(low="white", high="steelblue") +
  labs(x="scPDA-gating label", y="KNN-predicted label",
       title="Confusion Matrix for 314 Rescued Cells") +
  theme_minimal() + NoLegend() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle=45, hjust=1))

p_cv <- ggplot(knn_fit, aes(x = K, y = Accuracy)) +
  geom_line(color = "steelblue") + 
  geom_point(color = "steelblue") +
  labs(x = "K (number of neighbors)",
       y = "Accuracy",
       title = 'KNN Accuracy vs. Number of Neighbors') +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.ticks.x.top   = element_blank(),
    axis.text.x.top    = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right  = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  )

knn_all = (p_cv + plot_spacer() + p_conf) + 
  plot_layout(widths = c(1, 0.3, 1)) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = "bold", size = 20, hjust = 0.5, vjust = 0.5))
ggsave('fig_reprod/supp/knn.png', knn_all, width = 10, height = 4, dpi = 200)