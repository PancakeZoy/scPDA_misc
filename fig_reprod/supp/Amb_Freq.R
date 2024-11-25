require(rlist)
require(dplyr)
require(tibble)
require(stringr)
require(Seurat)
require(latex2exp)
require(data.table)
require(magrittr)
require(here)
'%ni%' = Negate('%in%')
setwd(here())

######################## 10K ########################
cells = readRDS('data/supp/10k_pos_prot.rds')
neg_prot = readRDS('data/supp/10k_neg_prot1.rds')
cells.profile = rowMeans(cells) %>% proportions()
empty.profile = rowMeans(neg_prot) %>% proportions()
names(cells.profile) = gsub("_TotalSeqB", "", names(cells.profile))
names(empty.profile) = gsub("_TotalSeqB", "", names(empty.profile))

jpeg('fig_reprod/supp/Prop_10K.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(cells.profile, empty.profile, pch=16,
     xlab='Protein Proportion Estimated from Cell Population',
     ylab='Protein Proportion Estimated from Empty Droplets',)
abline(0,1)
ann_prot = sort(abs(empty.profile-cells.profile)) %>% tail(5) %>% names()
text(x=cells.profile[ann_prot], y=empty.profile[ann_prot], labels=ann_prot, 
     pos=c(2,1,1,4,1), cex=0.8)
dev.off()

######################## 5K ########################
cells = readRDS('data/supp/5k_pos_prot.rds')
neg_prot = readRDS('data/supp/5k_neg_prot2.rds')
cells.profile = rowMeans(cells) %>% proportions()
empty.profile = rowMeans(neg_prot) %>% proportions()
names(cells.profile) = gsub("_TotalSeqB", "", names(cells.profile))
names(empty.profile) = gsub("_TotalSeqB", "", names(empty.profile))

jpeg('fig_reprod/supp/Prop_5K.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(cells.profile, empty.profile, pch=16,
     xlab='Protein Proportion Estimated from Cell Population',
     ylab='Protein Proportion Estimated from Empty Droplets',)
abline(0,1)
ann_prot = sort(abs(empty.profile-cells.profile)) %>% tail(5) %>% names()
text(x=cells.profile[ann_prot], y=empty.profile[ann_prot], labels=ann_prot, 
     pos=c(2,1,1,2,1), cex=0.8)
dev.off()

######################## 5 Prime ########################
cells = readRDS('data/supp/5`_pos_prot.rds')
neg_prot = readRDS('data/supp/5`_neg_prot2.rds')
cells.profile = rowMeans(cells) %>% proportions()
empty.profile = rowMeans(neg_prot) %>% proportions()
names(cells.profile) = gsub("_TotalSeqC", "", names(cells.profile))
names(empty.profile) = gsub("_TotalSeqC", "", names(empty.profile))

jpeg('fig_reprod/supp/Prop_5prime.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(cells.profile, empty.profile, pch=16,
     xlab='Protein Proportion Estimated from Cell Population',
     ylab='Protein Proportion Estimated from Empty Droplets',)
abline(0,1)
ann_prot = sort(abs(empty.profile-cells.profile)) %>% tail(5) %>% names()
text(x=cells.profile[ann_prot], y=empty.profile[ann_prot], labels=ann_prot, 
     pos=c(2,1,1,2,1), cex=0.8)
dev.off()
