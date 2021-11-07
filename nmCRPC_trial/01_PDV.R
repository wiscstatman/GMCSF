library(tidyverse)
library(fdrtool)
library(RColorBrewer)
library(Rtsne)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

# function to plot PCA loadings
PCload.func <- function(lab, U, D, x, y, pca.vec, title, cex_txt=1){
  Z <- U %*% diag(D)
  plot(Z[,x], Z[,y], type="n", main = title, las = 1,
       xlab = paste0("PC",x," (", round(pca.vec[x]*100, 1) ,"% variance explained)"),
       ylab = paste0("PC",y," (", round(pca.vec[y]*100, 1) ,"% variance explained)") 
  )
  text(Z[,x],Z[,y], lab, cex=cex_txt)
}


# function to count peptides at different FDR thresholds
count.func <- function(pval.vec, thresh.vec){
  counter <- NULL
  for (i in thresh.vec ){
    counter <- c(counter, length(pval.vec[pval.vec <= i]))
  }
  countab <- rbind(c("FDR threshold", thresh.vec),
                   c("Peptide counts", counter))
  return(countab)
}

####################################################################################### 
#                                      Setup Dataset                                  #
####################################################################################### 

dat <- read_tsv("Processed_aggregate_data.txt")
colnames(dat)
colnames(dat)[18:197] <- paste0( substr(colnames(dat)[18:197],1,7), substr(colnames(dat)[18:197],16,18) )
colnames(dat) #check

sample_dat <- read_csv("PDV_sample.csv")
sum(as.numeric(sample_dat$Array_ID %in% colnames(dat)[18:197])) == nrow(sample_dat) # check

# arrange samples according time and patients
sample_dat <- sample_dat %>%
  arrange(desc(Time), Patient_ID)

# check patients' post and pre records aligned properly
sum(as.numeric( sample_dat$Patient_ID[sample_dat$Time == "post"] == 
                  sample_dat$Patient_ID[sample_dat$Time == "pre"] )) == length(unique(sample_dat$Patient_ID))

# subset columns of dat
dat2 <- subset(dat, select = c("SEQ_ID", "PROBE_ID", sample_dat$Array_ID))

# check array ID are aligned between sample_dat and dat2
sum(as.numeric( colnames(dat2)[3:ncol(dat2)] == sample_dat$Array_ID )) == nrow(sample_dat)

# take log2 then take post - pre
dat3 = as.matrix(log2( dat2[,3:ncol(dat2)] ))
dat3 = dat3[,(ncol(dat3)/2 + 1):ncol(dat3)] - dat3[,1:(ncol(dat3)/2)] 
colnames(dat3) <- sample_dat$Patient_ID[1:ncol(dat3)]

# save probe ID and protein ID as well
dat3_probeID <- dat2$PROBE_ID
dat3_seqID <- dat2$SEQ_ID

rm(dat, dat2); gc()


####################################################################################### 
#                           Check Fluorescence Normalization                          #
####################################################################################### 

long_dat <- log2(dat2[,3:ncol(dat2)]) 
colnames(long_dat) <- paste( sample_dat$Patient_ID, sample_dat$Time, sep = "_" )
long_dat <- long_dat %>% 
  pivot_longer(cols = everything(), names_to = "patient_time", values_to = "log2_fluorescence") 

# check
nrow(long_dat) == nrow(dat3) * nrow(sample_dat)
head(long_dat)

# set fill color
long_dat$time <- factor(rep(sample_dat$Time, nrow(dat2)))

# sort order of patients in boxplot
long_dat$patient_time <- factor(long_dat$patient_time, 
                                levels = paste( sample_dat$Patient_ID, sample_dat$Time, sep = "_" ))

# set boxplot color
boxplot_pal <- brewer.pal(9,"Set1")[c(7,5)]
names(boxplot_pal) <- levels(long_dat$time )

tiff("fluorescence_normalization.tiff", res=600, compression = "lzw", width=15, height=10, units="in", pointsize = 10)
ggplot(long_dat, aes(x = patient_time, y = log2_fluorescence, fill = time)) +
  geom_boxplot(outlier.shape = ".") +
  scale_fill_manual(name = "time", values = boxplot_pal[levels(long_dat$time)]) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for Patients (pre- and post- treatment) ", 
       x = "Patient_Time", y = "Log2 Fluorescence Levels") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
dev.off()

rm(dat2, long_dat, boxplot_pal, dat); gc()

####################################################################################### 
#                               Preliminary Analysis -- Clustering                    #
####################################################################################### 

clust.dat <- t(dat3)
hc_average <- hclust(dist(clust.dat), "average")
dend <- as.dendrogram(hc_average)

# coloring
labels_colors(dend) <- cols[order.dendrogram(dend)]
par(cex=.4)
plot(dend, xlab = "Patients", ylab = "Dissimilarity")

####################################################################################### 
#                               Preliminary Analysis -- PCA                           #
####################################################################################### 

# svd
sv.dat <- sweep(t(dat3), 2, colMeans(t(dat3)), "-") # centering
sv <- svd(sv.dat)
V <- sv$v
D <- sv$d
U <- sv$u

# variance explained
pca.var <- D^2/sum(D^2) 
pca.cumvar <- cumsum(pca.var)

# plot PCA
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.6,0.4,0),cex=.5,
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func(colnames(dat3), U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
PCload.func(colnames(dat3), U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)

####################################################################################### 
#                              Preliminary Analysis -- t-SNE                          #
####################################################################################### 


# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf
tsnedat <- unname(t(dat3)) # dim(X) = N samples by D dimensions 
initdim <- 12
perplex <- 3 
lab <- row.names(t(dat3))

set.seed(10)
tsne_raw <- Rtsne(tsnedat, dims = 2L, initial_dims = initdim, perplexity = perplex,
                  theta = 0, check_duplicates = F, PCA = T)

# t-SNE plot
# plot(tsne_raw$Y, ylab = "", xlab = "", main = "tSNE plot", col = cols, pch = shapes)
# legend('bottomright', pch = shp, col = pal, cex = 0.5, c("A",  "B") )
plot(tsne_raw$Y,t = "n", ylab = "", xlab = "")
text(tsne_raw$Y, labels=lab)

####################################################################################### 
#                                      test if mu != 0                                 #
#######################################################################################

# initialize
pval <- rep(NA, nrow(dat3))
est <- rep(NA, nrow(dat3))
StdErr <- rep(NA, nrow(dat3))
Wilcox_pval <- rep(NA, nrow(dat3))

for(i in 1:nrow(dat3)){
  fit1 <- t.test(as.numeric(dat3[i, ]),alternative = "two.sided")
  pval[i] <- fit1$p.value
  est[i] <- fit1$estimate
  StdErr[i] <- fit1$stderr
  Wilcox_pval[i] <- wilcox.test(as.numeric(dat3[i, ]), alternative = "two.sided",mu=0)$p.value
  
  if(i %% 1000 == 0 | i == nrow(dat3)){
    print(i)
  }
}

par(mfrow=c(1,2))
hist(pval, xlab = "t-test p-val", freq = F, col = "grey", breaks= 50, 
     main = "t-test p-value histogram")
hist(Wilcox_pval, xlab = "Wilcox p-val", freq = F, col = "grey", 
     main = "Wilcoxon-test p-value histogram")
par(mfrow=c(1,1))

BHpval_ttest <- p.adjust(pval, method = "BH") # t-test BH-adj p-val
BHpval_wilcox <- p.adjust(Wilcox_pval, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_ttest, seq(.15, .25,by=.01))
count.func(BHpval_wilcox, seq(.15, .25,by=.01))

####################################################################################### 
#                          Compare against JITC List                                  #
####################################################################################### 

# JITC list
JITC_list <- read_csv("JITC_PAP_Longitudinal.csv")
JITC_ind <- which(dat3_probeID %in% JITC_list$PROBE_ID)

sum(as.numeric( est[JITC_ind] >= 1 ))/nrow(JITC_list) 
# 99% of JITC-identified peptides showed at least one fold-change increase 

sum(as.numeric( (est[JITC_ind] >= 1)&(pval[JITC_ind] <= .05) ))/nrow(JITC_list) 
# 81.44% of JITC-identified peptides showed at least one fold-change increase with raw p-values <= .05

sum(as.numeric( est[JITC_ind] < 0 )) 
# all JITC-identified peptides show increase in peptide activities 

signif_ind <- (est[JITC_ind] >= 1)&(pval[JITC_ind] <= .05)

png("volcano_plot.png", res=600, width = 6, height = 6, units = "in")
plot(x = est, y = -log10(pval), pch = ".",
     xlab = "fold changes", ylab = "-log10(t-test p-values)", main = "volcano plot" )
lines(x = est[JITC_ind], y = (-log10(pval))[JITC_ind], 
      type = "p", pch = 20, cex = .3, col = "red")
# lines(x = est[JITC_ind][signif_ind], y = (-log10(pval))[JITC_ind][signif_ind], 
#       type = "p", pch = 20, cex = .3, col = "blue")
abline(v = 1, col = "blue", lty = 2)
abline(h = -log10(.05), col = "blue", lty = 2)
dev.off()