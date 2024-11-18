require(rlist)
require(dplyr)
require(tibble)
require(stringr)
require(latex2exp)
require(data.table)
require(magrittr)
require(here)
'%ni%' = Negate('%in%')
setwd(here())

cells = readRDS('data/supp/dsb_cells.rds')
lib = readRDS("data/supp/dsb_lib.rds")[[1]]
rna_size = lib$nUMI_mRNA; lib=lib[,-1]; lib=lib %>% column_to_rownames('bc')
hash = readRDS("data/supp/dsb_hash.rds")[[1]] %>% as.matrix()
lib.profile = colMeans(log1p(lib)); hash.profile = rowMeans(log1p(hash))

coef = lm(hash.profile~lib.profile)$coefficients
test = t.test(lib.profile, hash.profile, paired = TRUE)

jpeg('fig_reprod/supp/Hash_vs_Lib.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(lib.profile, hash.profile, pch=16, xlab='Empty Droplets Defined by Library Size', 
     ylab='Empty Droplets Defined by Hashing', main=NULL)
abline(0,1)
mtext(sprintf('Hashing = %.2f + %.2f * Library', coef[1], coef[2]), side=3, line=-1)
mtext(TeX("$R^2: $ 0.00"), side=1, line=-1)
dev.off()


######################## Thresholds ########################
neg_prot1 = readRDS('data/supp/10k_neg_prot1.rds') %>% log1p
neg_prot2 = readRDS('data/supp/10k_neg_prot2.rds') %>% log1p
neg_prot3 = readRDS('data/supp/10k_neg_prot3.rds') %>% log1p
neg_prot4 = readRDS('data/supp/10k_neg_prot4.rds') %>% log1p

mean.lib1 = rowMeans(neg_prot1)
mean.lib2 = rowMeans(neg_prot2)
mean.lib3 = rowMeans(neg_prot3)
mean.lib4 = rowMeans(neg_prot4)

jpeg('fig_reprod/supp/Tr1_vs_Tr2.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(mean.lib1, mean.lib2, pch=16, 
     xlim=c(0,max(mean.lib1)*1.2),
     ylim = c(0, max(mean.lib2)*1.2), 
     xlab='Empty Droplets Defined by Threshold 1',
     ylab = 'Empty Droplets Defined by Threshold 2')
abline(0,1)
coef = lm(mean.lib2~mean.lib1)$coefficients
test = t.test(mean.lib1, mean.lib2, paired = TRUE)$p.value
mtext(sprintf('Threshold 2 = %.2f + %.2f * Threshold 1', coef[1], coef[2]), side=3, line=-1)
mtext(TeX("$R^2: $ 0.00"), side=1, line=-1)
dev.off()

jpeg('fig_reprod/supp/Tr2_vs_Tr3.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(mean.lib2, mean.lib3, pch=16, 
     xlim=c(0,max(mean.lib2)*1.5),
     ylim = c(0, max(mean.lib3)*2), 
     xlab='Empty Droplets Defined by Threshold 2',
     ylab = 'Empty Droplets Defined by Threshold 3')
abline(0,1)
coef = lm(mean.lib3~mean.lib2)$coefficients
test = t.test(mean.lib3, mean.lib2, paired = TRUE)$p.value
mtext(sprintf('Threshold 3 = %.2f + %.2f * Threshold 2', coef[1], coef[2]), side=3, line=-1)
mtext(TeX("$R^2: $ 0.00"), side=1, line=-1)
dev.off()

jpeg('fig_reprod/supp/Tr3_vs_Tr4.jpeg', width=800, height=800, res=150)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(mean.lib3, mean.lib4, pch=16, 
     xlim=c(0,max(mean.lib3)*1.5),
     ylim = c(0, max(mean.lib4))*1.2, 
     xlab='Empty Droplets Defined by Threshold 3',
     ylab = 'Empty Droplets Defined by Threshold 4')
abline(0,1)
coef = lm(mean.lib4~mean.lib3)$coefficients
test = t.test(mean.lib3, mean.lib4, paired = TRUE)$p.value
mtext(sprintf('Threshold 4 = %.2f + %.2f * Threshold 3', coef[1], coef[2]), side=3, line=-1)
mtext(TeX("$R^2: $ 0.00"), side=1, line=-1)
dev.off()