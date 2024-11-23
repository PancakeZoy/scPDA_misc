require(cowplot)
require(data.table)
require(dplyr)
require(dsb)
library(forcats)
require(ggplot2)
require(ggrepel)
library(ggridges)
require(grid)
require(gridExtra)
require(magick)
require(mclust)
require(patchwork)
require(pROC)
require(qpdf)
require(reticulate)
require(rhdf5)
require(rlist)
require(scales)
require(Seurat)
require(SeuratData)
require(tibble)
library(pROC)
library(tidyr)

'%ni%' = Negate('%in%')
log10p = function(x){log10(1+x)}

ridge = function(obj, prot, slot, log, groupby, line=NULL, color){
  p = RidgePlot(obj, features=prot, group.by = groupby, slot=slot, log=log)+NoLegend()
  if (!is.null(line)){
    v = sapply(line, function(x) x[prot])
    p = p + geom_vline(xintercept = v, col=color, size=1)
  }
  return(list(p))
}

getambient = function(matrix, G=1:2, quan=0.9){
  bg_noise_list = apply(matrix, 1, function(x) {Mclust(x, G = G, warn = FALSE, verbose = FALSE)})
  bg_noise = c()
  for (prot in names(bg_noise_list)){
    mu = bg_noise_list[[prot]]$parameters$mean
    k = length(mu)
    if (k==1){
      q = quantile(matrix[prot,], quan)
      mu = mean(matrix[prot,]<=q)
    }else{
      mu = mu[1]
    }
    bg_noise=c(bg_noise, mu)
  }
  names(bg_noise)=names(bg_noise_list)
  ambient_profile = bg_noise %>% as.data.frame()
  colnames(ambient_profile) = 'Ambient'
  return(ambient_profile)
}

rgb_alpha <- function(color_name, alpha=1) {
  cols <- col2rgb(color_name)
  r <- cols[1,] / 255
  g <- cols[2,] / 255
  b <- cols[3,] / 255
  return(rgb(r, g, b, alpha=alpha))
}

assign_prob = function(counts, prot, pos, types){
  label = if_else(types %in% pos, 1, 0)
  gmm = Mclust(log1p(counts[prot,]), G=2, warn = FALSE, verbose = FALSE)
  roc <- roc(label, gmm$z[,1]); auc = sprintf('%.2f', auc(roc))
  return(list(roc=roc, auc=auc))
}

roc_auc = function(prob, prot, pos, types){
  label = if_else(types %in% pos, 1, 0)
  roc <- roc(label, prob[prot,]); auc = sprintf('%.2f', auc(roc))
  return(list(roc=roc, auc=auc))
}

curve = function(prob.list, prot){
  n = length(prob.list); colors=hue_pal()(n); 
  methods=names(prob.list); auc=sapply(prob.list, function(x) x$auc)
  plot(prob.list[[1]]$roc, col=colors[1], main=sprintf('%s ROC Curve', prot))
  for (i in 2:n){
    lines(prob.list[[i]]$roc, col=colors[i])
  }
  legend('bottomright', col=colors, lty=1, bty="n", cex=0.8, lwd=2,
         legend=paste0(methods, ' (AUC:', auc,')'))
}

overlap = function(counts, prot, pos, pos.name, title, types, neg=NA, xlab='log(counts + 1)',
                   ifprob=FALSE, col.pos='navyblue', col.neg='orange', neg.name='Others'){
  col.pos = rgb_alpha(col.pos, 0.4); col.neg=rgb_alpha(col.neg, 0.4)
  if (all(is.na(neg))){
    alltypes = unique(types)
    neg = setdiff(alltypes, pos)
  }
  counts=log1p(counts)
  pos_cells = types %in% pos; neg_cells = types %in% neg
  pos_counts = counts[prot, pos_cells]; neg_counts = counts[prot, neg_cells]
  lfc = log(mean(pos_counts)/mean(neg_counts),2)
  main = if (ifprob) sprintf('%s %s', prot, title) else sprintf('%s %s (%.2f)', prot, title, lfc)
  hist(neg_counts, probability = TRUE, col=col.neg, breaks=50, 
       main=main, xlab=xlab)
  hist(pos_counts, probability = TRUE, add=TRUE, col=col.pos, breaks=50)
  legend('topright', legend=c(neg.name, pos.name), fill=c(col.neg, col.pos), bty="n", cex=0.8)
  return(lfc)
}

overlap.gg = function(counts, prot, pos, pos.name, title, types, neg=NA, xlab='log(counts + 1)',
                      ifLFC=TRUE, AUC=NULL, col.pos='navyblue', col.neg='orange', neg.name='Others'){
  col.pos = rgb_alpha(col.pos, 0.4); col.neg=rgb_alpha(col.neg, 0.4)
  if (all(is.na(neg))){
    alltypes = unique(types)
    neg = setdiff(alltypes, pos)
  }
  counts=log1p(counts)
  pos_cells = types %in% pos; neg_cells = types %in% neg
  pos_counts = counts[prot, pos_cells]; neg_counts = counts[prot, neg_cells]
  lfc = log(mean(pos_counts)/mean(neg_counts),2)
  
  descriptors = c()
  if (ifLFC) descriptors <- c(descriptors, sprintf('LFC:%.2f', lfc))
  if (!is.null(AUC)) descriptors <- c(descriptors, sprintf('AUC:%.2f', AUC))
  main = if (length(descriptors)==0) title else paste0(title, ' (', paste(descriptors, collapse='/'), ')')
  
  data <- data.frame(counts = c(neg_counts, pos_counts),
                     group = c(rep(neg.name, length(neg_counts)), 
                               rep(pos.name, length(pos_counts))))
  p1_w_legend = ggplot(data, aes(x=counts, fill=group)) +
    geom_histogram(aes(y=after_stat(density)), position="identity", bins=100, 
                   alpha=0.6, color="black", linewidth=0.1) +
    scale_fill_manual(values=c(col.neg, col.pos)) + theme_minimal() +
    scale_y_continuous(expand = c(0,0)) + labs(x=xlab, title = main) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1), 
          legend.box.just="right", 
          legend.margin=margin(0,0,0,0),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line = element_line(colour = "black")) + 
    guides(fill=guide_legend(title=NULL))
  return(p1_w_legend)
}

overlap3 = function(counts, prot, pos1, pos2, neg, title, col.neg='navyblue', types,
                    col.pos1='orange', col.pos2='red3', pos.name, neg.name='Others'){
  col.neg=rgb_alpha(col.neg, 0.4); col.pos1 = rgb_alpha(col.pos1, 0.4); col.pos2 = rgb_alpha(col.pos2, 0.4)
  pos_cells1 = types %in% pos1; pos_cells2 = types %in% pos2
  pos_counts1 = counts[prot, pos_cells1]; pos_counts2 = counts[prot, pos_cells2]
  neg_cells = types %in% neg; neg_counts = counts[prot, neg_cells]
  hist(log1p(neg_counts), probability = TRUE, col=col.neg, breaks=50, 
       main=paste(prot, title), xlab='log(counts + 1)')
  hist(log1p(pos_counts1), probability = TRUE, add=TRUE, col=col.pos1, breaks=50)
  hist(log1p(pos_counts2), probability = TRUE, add=TRUE, col=col.pos2, breaks=50)
  legend('topright', legend=c(neg.name, pos.name), fill=c(col.neg, col.pos1, col.pos2), bty="n", cex=0.8)
}

overlap3.gg = function(counts, prot, pos1, pos2, neg, pos.name, title, types, xlab='log(counts + 1)',
                       col.neg='navyblue', col.pos1='orange', col.pos2='red3', neg.name='Others'){
  col.neg=rgb_alpha(col.neg, 0.4); col.pos1 = rgb_alpha(col.pos1, 0.4); col.pos2 = rgb_alpha(col.pos2, 0.4)
  pos_cells1 = types %in% pos1; pos_cells2 = types %in% pos2; counts=log1p(counts)
  pos_counts1 = counts[prot, pos_cells1]; pos_counts2 = counts[prot, pos_cells2]
  neg_cells = types %in% neg; neg_counts = counts[prot, neg_cells]
  lfc.1 = log(mean(pos_counts1)/mean(neg_counts),2) %>% round(2)
  lfc.2 = log(mean(pos_counts2)/mean(neg_counts),2) %>% round(2)
  main = bquote(.(title) ~ (LFC[CD4]: .(lfc.1) / LFC[Mono] ~ .(lfc.2)))
  data <- data.frame(counts = c(neg_counts, pos_counts1, pos_counts2),
                     group = factor(c(rep(neg.name, length(neg_counts)), 
                                      rep('CD4 T Cells', length(pos_counts1)),
                                      rep('Monocytes', length(pos_counts2))), 
                                    levels=c('Others', 'Monocytes', 'CD4 T Cells')))
  p1_w_legend = ggplot(data, aes(x=counts, fill=group)) +
    geom_histogram(aes(y=after_stat(density)), position="identity", bins=100, 
                   alpha=0.6, color="black", linewidth=0.1) +
    scale_fill_manual(values=c(col.neg, col.pos1, col.pos2)) + theme_minimal() +
    scale_y_continuous(expand = c(0,0)) + labs(x=xlab, title = main) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1), 
          legend.box.just="right", 
          legend.margin=margin(0,0,0,0),
          plot.title = element_text(hjust = 0.5, size=10),
          axis.line = element_line(colour = "black")) + 
    guides(fill=guide_legend(title=NULL))
  return(p1_w_legend)
}


curve.gg = function(prob.list){
  n = length(prob.list); auc=sapply(prob.list, function(x) x$auc)
  colors=hue_pal()(n); methods=names(prob.list)
  for (i in 1:n){
    methods[i] = sprintf('%s (AUC: %s)', methods[i], auc[methods[i]])
  }
  sens = lapply(prob.list, function(x) data.frame(sens=x$roc$sensitivities))
  spec = lapply(prob.list, function(x) data.frame(spec=x$roc$specificities))
  roc_data <- do.call(rbind, lapply(1:n, function(i) {
    data.frame(method = methods[i],sens = sens[[i]]$sens,spec = spec[[i]]$spec)
  }))
  roc_data$method = factor(roc_data$method, levels=methods)
  ggplot(roc_data, aes(x = sens, y = spec, color = method)) +
    geom_line() + scale_color_manual(values = colors) + 
    geom_abline(slope=1, intercept=1, color='grey')+
    labs(title = "ROC Curves", x = "Sensitivity", y = "Specificity") +
    theme_minimal() + scale_x_reverse() + theme.border + title.center + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(1, 0), 
          legend.justification = c(1, 0),
          legend.box.just = "right") + 
    guides(color=guide_legend(title=NULL))
}


overlap.f1 = function(counts, prot, pos, pos.name, title, types, neg='Neg Control', xlab='Counts',
                      col.pos='navyblue', col.neg='orange', neg.name='Unstained Cells'){
  colors = c(rgb_alpha(col.pos, 0.4), rgb_alpha(col.neg, 0.4))
  names(colors) = c(pos.name, neg.name)
  log10_1p_trans <- trans_new(
    name = "log10_1p",
    transform = function(x) log10(1 + x),
    inverse = function(x) 10^x - 1
  )
  breaks <- c(0, 10, 100, 1000)
  pos_cells = types %in% pos; neg_cells = types %in% neg
  pos_counts = counts[prot, pos_cells]; neg_counts = counts[prot, neg_cells]
  data <- data.frame(counts = c(neg_counts, pos_counts),
                     group = c(rep(neg.name, length(neg_counts)), 
                               rep(pos.name, length(pos_counts))))
  data$group <- factor(data$group, levels = c(neg.name, pos.name))
  p1_w_legend = ggplot(data, aes(x=counts, fill=group)) +
    geom_histogram(aes(y=after_stat(density)), position="identity", bins=100, 
                   alpha=0.6, color="black", linewidth=0.1) +
    scale_x_continuous(trans = log10_1p_trans, breaks = breaks, labels = breaks)+
    scale_fill_manual(values=colors) + theme_minimal() +
    scale_y_continuous(expand = c(0,0)) + labs(x=xlab, title = title) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1), 
          legend.box.just="right", 
          legend.margin=margin(0,0,0,0),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black")) + 
    guides(fill=guide_legend(title=NULL))
  return(p1_w_legend)
}


overlap.negpro = function(counts, prot, title, xlab='Counts', col='navyblue'){
  col = rgb_alpha(col, 0.4)
  log10_1p_trans <- trans_new(
    name = "log10_1p",
    transform = function(x) log10(1 + x),
    inverse = function(x) 10^x - 1
  )
  breaks <- c(0, 10, 100, 1000)
  data <- data.frame(counts=counts[prot, ])
  p1_w_legend = ggplot(data, aes(x=counts)) +
    geom_histogram(aes(y=after_stat(density)), position="identity", bins=100, fill=col,
                   alpha=0.6, color="black", linewidth=0.1) + theme_minimal() +
    scale_x_continuous(trans = log10_1p_trans, breaks = breaks, labels = breaks)+ 
    scale_y_continuous(expand = c(0,0)) + labs(x=xlab, title = title) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1), 
          legend.box.just="right", 
          legend.margin=margin(0,0,0,0),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black")) + 
    guides(fill=guide_legend(title=NULL))
  return(p1_w_legend)
}