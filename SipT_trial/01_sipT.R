library(tidyverse)
library(RColorBrewer)
library(fdrtool)
library(gridExtra)
library(dendextend) # dendrogram
library(Rtsne)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

# specified color and shape scheme 
pal <- c("navy", "darkorange1")
names(pal) <- c("A",  "B")
shp <- c(8, 16)
names(shp) <- names(pal)

# function to plot PCA loadings
PCload.func <- function(cols, shapes, U, D, x, y, pca.vec, title){
  Z <- U %*% diag(D)
  plot(Z[,x], Z[,y], col = cols, pch = shapes, main = title, las = 1,
       xlab = paste0("PC",x," (", round(pca.vec[x]*100, 1) ,"% variance explained)"),
       ylab = paste0("PC",y," (", round(pca.vec[y]*100, 1) ,"% variance explained)") 
  )
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

# load t-test & Wilcoxon rank-sum test results
load("sipT_results.RData")
load("sipT_muB.RData")
load("sipT_all.RData")

####################################################################################### 
#                                      Setup Dataset                                  #
####################################################################################### 

# read array data
dat <- read_csv("MCV_data.csv")

# check
colnames(dat)
length(unique(dat$PROBE_ID)) == nrow(dat)
length(unique(dat$PROBE_SEQUENCE)) == nrow(dat)

# read sample data
samp <- read_csv("MCV_samples.csv")

# check
length(unique(samp$Array_ID)) == nrow(samp)
samp %>%  group_by(Group, Time) %>% summarise(n=n())
sum(as.numeric( samp$Array_ID %in% colnames(dat) )) == nrow(samp)

# rearrange rows in samp
samp <- samp %>%
  arrange(desc(Time), Group, desc(Period), Patient_ID)

# check
sum(as.numeric( samp$Patient_ID[1:(nrow(samp)/2)] == samp$Patient_ID[(nrow(samp)/2+1): nrow(samp)] )) == nrow(samp)/2

# rearrange columns in dat
arr_iii <- match( samp$Array_ID, colnames(dat[,4:ncol(dat)]) )
dat <- dat[, c(1:3, arr_iii+3)]

# check
sum(as.numeric( colnames(dat[, 4:ncol(dat)]) == samp$Array_ID )) == nrow(samp)
rm(arr_iii); gc()

# take log2 then find take (post - pre)
dat2 <- log2( dat[,4:ncol(dat)] )
dat2 <- dat2[,(ncol(dat2)/2 + 1):ncol(dat2)] - dat2[,1:(ncol(dat2)/2)] 
colnames(dat2) <- samp$Patient_ID[1:ncol(dat2)]

# check
ncol(dat2) == length(unique(colnames(dat2)))
samp$Period[1:ncol(dat2)]
samp$Group[1:ncol(dat2)]

dat2 <- cbind(dat$PROBE_ID, dat$SEQ_ID, dat2)

# colors and shapes for the visualization techniques
cols = pal[ match(samp$Group[samp$Time == "Post"], names(pal)) ]
shapes = shp[ match(samp$Group[samp$Time == "Post"], names(shp)) ]

####################################################################################### 
#                           Check Fluorescence Normalization                          #
####################################################################################### 

long_dat <- log2(dat[,4:ncol(dat)]) 
colnames(long_dat) <- paste( samp$Patient_ID, samp$Time, sep = "_" )
long_dat <- long_dat %>% 
  pivot_longer(cols = everything(), names_to = "id_time", values_to = "log2_fluorescence") 

# check
nrow(long_dat) == nrow(dat2) * nrow(samp)
head(long_dat)

# set fill color
long_dat$group_time <- factor(rep(paste(samp$Group,samp$Time,sep="_"), nrow(dat2)))

# sort order of patients in boxplot
long_dat$id_time <- factor(long_dat$id_time, levels = paste( samp$Patient_ID, samp$Time, sep = "_" ))

# set boxplot color
boxplot_pal <- brewer.pal(9,"Set1")[c(2,6,7,5)]
names(boxplot_pal) <- levels(long_dat$group_time )

tiff("fluorescence_normalization.tiff", res=600, compression = "lzw", width=15, height=10, units="in", pointsize = 10)
ggplot(long_dat, aes(x = id_time, y = log2_fluorescence, fill = group_time)) +
  geom_boxplot(outlier.shape = ".") +
  scale_fill_manual(name = "group_time", values = boxplot_pal[levels(long_dat$group_time)]) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for Patients (pre- and post- treatment) ", 
       x = "Patient_Time", y = "Log2 Fluorescence Levels") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
dev.off()

####################################################################################### 
#                               Preliminary Analysis -- Clustering                    #
####################################################################################### 

clust.dat <- t(dat2[,3:ncol(dat2)])
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
sv.dat <- sweep(t(dat2[,3:ncol(dat2)]), 2, colMeans(t(dat2[,3:ncol(dat2)])), "-") # centering
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
PCload.func(cols, shapes, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('bottomright', pch = shp, col = pal, c("A",  "B") )
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)
legend('bottomright', pch = shp, col = pal, c("A",  "B") )

# verify the outlier points
pca_Z <- U %*% diag(D)
colnames(dat2[,3:ncol(dat2)])[which(pca_Z[,2]<=-350)] # ID010
colnames(dat2[,3:ncol(dat2)])[which(pca_Z[,1]>150)] # ID004, 005, 006, 007, 008


####################################################################################### 
#                              Preliminary Analysis -- t-SNE                          #
####################################################################################### 

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf
tsnedat <- unname(t(dat2[,3:ncol(dat2)])) # dim(X) = N samples by D dimensions 
initdim <- 18 
perplex <- 5 
lab <- row.names(t(dat2[,3:ncol(dat2)]))

set.seed(10)
tsne_raw <- Rtsne(tsnedat, dims = 2L, initial_dims = initdim, perplexity = perplex,
                  theta = 0, check_duplicates = F, PCA = T)

# t-SNE plot
# plot(tsne_raw$Y, ylab = "", xlab = "", main = "tSNE plot", col = cols, pch = shapes)
# legend('bottomright', pch = shp, col = pal, cex = 0.5, c("A",  "B") )
plot(tsne_raw$Y,t = "n", ylab = "", xlab = "")
text(tsne_raw$Y, labels=lab, col=cols)
legend('bottomright', col = pal, cex = 0.5, c("A",  "B") )

# check which 5 patients clustered away
# colnames(dat2[,3:ncol(dat2)])[which(tsne_raw$Y[,1]>100)] # ID004, 005, 006, 007, 008

####################################################################################### 
#                      2-sample t-test & Wilcoxon Rank-Sum Test                       #
#                             to Assess Treatment Effect                              #
####################################################################################### 

rm(dat, long_dat, pal); gc()

# getting variables ready
gr <- samp$Group[1:(nrow(samp)/2)]
pd <- samp$Period[1:(nrow(samp)/2)]
gr_drop <- gr[pd != "3-month"]
dat2_drop <- dat2[,3:20]
dat2_drop <- dat2_drop[,(pd != "3-month")]

# check
table(gr); table(gr_drop)

#---------------------------------------------------------------------------------------
# t-test

# initialize
pval <- rep(NA, nrow(dat2))
est <- matrix(NA, nrow = nrow(dat2), ncol = 2)
StdErr <- rep(NA, nrow(dat2))

pval_drop <- rep(NA, nrow(dat2))
est_drop <- matrix(NA, nrow = nrow(dat2), ncol = 2)
StdErr_drop <- rep(NA, nrow(dat2))

for(i in 1:nrow(dat2)){
  fit1 <- t.test( as.numeric(dat2[i, 3:20]) ~ gr )
  pval[i] <- fit1$p.value
  est[i,] <- fit1$estimate
  StdErr[i] <- fit1$stderr 
  
  fit_drop <- t.test( as.numeric(dat2_drop[i,]) ~ gr_drop )
  pval_drop[i] <- fit_drop$p.value
  est_drop[i,] <- fit_drop$estimate
  StdErr_drop[i] <- fit_drop$stderr 

  if(i %% 1000 == 0 | i == nrow(dat2)){
    print(i)
  }
}

#---------------------------------------------------------------------------------------
# Wilcoxon rank-sum test

# initialize
Wilcox_pval <- rep(NA, nrow(dat2))
Wilcox_pval_drop <- rep(NA, nrow(dat2))

for(i in 1:nrow(dat2)){
  Wilcox_pval[i] <- wilcox.test( as.numeric(dat2[i, 3:20]) ~ gr )$p.value
  Wilcox_pval_drop[i] <- wilcox.test( as.numeric(dat2_drop[i,]) ~ gr_drop )$p.value
  if(i %% 1000 == 0 | i == nrow(dat2)){
    print(i)
  }
}

save(pval, est, StdErr,
     pval_drop, est_drop, StdErr_drop,
     Wilcox_pval, Wilcox_pval_drop,
     file = "sipT_results.RData")

####################################################################################### 
#              verify if 3-month period significantly different from 6-month          #
####################################################################################### 

# compare t-test p-val
tiff("compare_ttestpval.tiff", res=600, compression = "lzw", width=6, height=6, units="in")
plot(x = pval, y = pval_drop, type = "p", pch = ".", 
     xlab = "t-test pval (all patients)", ylab = "t-test pval (exclude 3-month)",
     main = "compare t-test p-values")
abline(a=0,b=1,lty = 2, col = "blue")
dev.off()

# compare Wilcox p-val
tiff("compare_wilcoxpval.tiff", res=600, compression = "lzw", width=6, height=6, units="in")
plot(x = Wilcox_pval, y = Wilcox_pval_drop, type = "p", pch = ".", 
     xlab = "Wilcox pval (all patients)", ylab = "Wilcox pval (exclude 3-month)",
     main = "compare Wilcoxon  test p-values")
abline(a=0,b=1,lty = 2, col = "blue")
dev.off()

# compare within-treatment-group diff = post - pre 
# between 6-month (average) and 3-month (individual)
tiff("compare_diff.tiff", res=600, compression = "lzw", width=10, height=6, units="in")
par(mfrow=c(1,2))
plot(x = est_drop[,1], y = dat2[,3:20][,(gr=="A") & (pd=="3-month")], 
     type = "p", pch = ".", cex.lab = .7, cex.main = .8,
     xlab = "average change (post - pre) in log2(fluorescence) among 6-month patients", 
     ylab = "change (post - pre) in log2(fluorescence) for 3-month patient",
     main = "compare change (post - pre) in log2(fluorescence) \nin group A (sipT only)")
abline(a=0,b=1,lty = 2, col = "blue")
plot(x = est_drop[,2], y = dat2[,3:20][,(gr=="B") & (pd=="3-month")], 
     type = "p", pch = ".", cex.lab = .7,cex.main = .8,
     xlab = "average change in log2(fluorescence) among 6-month patients", 
     ylab = "change in log2(fluorescence) for 3-month patient",
     main = "compare change (post - pre) in log2(fluorescence) \nin group B (sipT + vaccine)")
abline(a=0,b=1,lty = 2, col = "blue")
dev.off()

####################################################################################### 
#                       Maybe drop 3-month patients and proceed?                      #
####################################################################################### 

# get p-val histograms
tiff("pval_exclude3month.tiff", res=600, compression = "lzw", width=10, height=6, units="in")
par(mfrow=c(1,2))
hist(pval_drop, xlab = "t-test p-val", freq = F, 
     main = "t-test p-value histogram \n(Exclude 3-month patients)")
hist(Wilcox_pval_drop, xlab = "Wilcox p-val", freq = F, 
     main = "Wilcoxon-test p-value histogram \n(Exclude 3-month patients)")
dev.off()

# get BH-adjusted p-val
BHpval_drop <- p.adjust(pval_drop, method = "BH") # t-test BH-adj p-val
Wilcox_BHpval_drop <- p.adjust(Wilcox_pval_drop, method = "BH") # wilcoxon BH-adj p-val

count.func(BHpval_drop, seq(.1, .5,by=.1))
count.func(Wilcox_BHpval_drop, seq(.1, .5,by=.1))

# no peptide is significant after BH adjustment (for both t-test and Wilcoxon test)...

####################################################################################### 
#                         If we include all patients anyway...                        #
####################################################################################### 

# get p-val histograms
tiff("pval_AllPatients.tiff", res=600, compression = "lzw", width=10, height=6, units="in")
par(mfrow=c(1,2))
hist(pval, xlab = "t-test p-val", freq = F, 
     main = "t-test p-value histogram \n(include all patients)")
hist(Wilcox_pval, xlab = "Wilcox p-val", freq = F, 
     main = "Wilcoxon-test p-value histogram \n(include all patients)")
dev.off()

# get BH-adjusted p-val
BHpval <- p.adjust(pval, method = "BH") # t-test BH-adj p-val
Wilcox_BHpval <- p.adjust(Wilcox_pval, method = "BH") # wilcoxon BH-adj p-val

count.func(BHpval, seq(.1, .5,by=.1))
count.func(Wilcox_BHpval, seq(.1, .5,by=.1))

# no peptide is significant after BH adjustment (for both t-test and Wilcoxon test)...

####################################################################################### 
#                          Compare against JITC List                                  #
####################################################################################### 
# get data frame ready

# JITC list
JITC_list <- read_csv("JITC_PAP_Longitudinal.csv")

# get group and overall mean, stdev & studentized fold changes
ind_A <- which( samp$Group[samp$Time=="Post"] == "A" )
ind_B <- which( samp$Group[samp$Time=="Post"] == "B" )

# overall mean, stdev & studentized fold changes
foldchange_overall <- t(apply( dat2[,3:ncol(dat2)], 1, function(x){
  mean_fold <- mean(x)
  stdev_fold <- sd(x)
  studentized_fold <- mean_fold/stdev_fold  
  return( c( mean_fold, stdev_fold, studentized_fold  ) )
} ))

# group A mean, stdev & studentized fold changes
foldchange_A <- t(apply( dat2[,3:ncol(dat2)][,ind_A], 1, function(x){
  mean_fold <- mean(x)
  stdev_fold <- sd(x)
  studentized_fold <- mean_fold/stdev_fold  
  return( c( mean_fold, stdev_fold, studentized_fold  ) )
} ))

# group B mean, stdev & studentized fold changes
foldchange_B <- t(apply( dat2[,3:ncol(dat2)][,ind_B], 1, function(x){
  mean_fold <- mean(x)
  stdev_fold <- sd(x)
  studentized_fold <- mean_fold/stdev_fold  
  return( c( mean_fold, stdev_fold, studentized_fold  ) )
} ))

# combine into single data frame
compare_df <- data.frame(
  dat$PROBE_ID, 
  round(foldchange_A,4), 
  round(foldchange_B,4), 
  round(foldchange_overall,4) 
)
colnames(compare_df) <- c(
  "PROBE_ID",
  paste0(c("mean_","stdev_","studentized_"),"A"),
  paste0(c("mean_","stdev_","studentized_"),"B"),
  paste0(c("mean_","stdev_","studentized_"),"overall")
)
compare_df <- compare_df %>% arrange(desc(studentized_overall))

compare_df$JITC_signif <- as.numeric(compare_df$PROBE_ID %in% JITC_list$PROBE_ID)

write.table(compare_df, "fold_change_all_peptides.csv", sep = ",", col.names = T, row.names = F)



#---------------------------------------------------------------------------------
# compare 

compare_df <- read_csv("fold_change_all_peptides.csv")

compare_plot.func <- function(colmn, xlab){
  df <- data.frame(
    val = colmn,
    JITC_signif = as.factor(compare_df$JITC_signif)
  )
  
  grid.arrange(
    ggplot(df, aes(x = val, fill = JITC_signif)) +
      geom_histogram(aes(y=..density..), bins = 100, position = "identity", alpha = .4) +
      labs(x = xlab, title = paste0("Density histograms of \n", xlab)) + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "bottom"),
    ggplot(df, aes(x = val, fill = JITC_signif)) +
      geom_histogram(bins = 100, position = "identity", alpha = .4) +
      labs(x = xlab, title = paste0("Frequency histograms of \n", xlab)) + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "bottom"),
    ncol = 2
  )
}

compare_plot.func(compare_df$mean_A, "mean fold changes (group A)")
compare_plot.func(compare_df$studentized_A, "studentized fold changes (group A)")

compare_plot.func(compare_df$mean_B, "mean fold changes (group B)")
compare_plot.func(compare_df$studentized_B, "studentized fold changes (group B)")

compare_plot.func(compare_df$mean_overall, "mean fold changes (all patients)")
compare_plot.func(compare_df$studentized_overall, "studentized fold changes (all patients)")


####################################################################################### 
#                         Check p-value histogram & FDR control                       #
#                         by subsetting only those in JITC list                       #
#######################################################################################

JITC_signif <- (dat2$'dat$PROBE_ID' %in% JITC_list$PROBE_ID)

# include all patients
par(mfrow=c(1,2))
hist(pval[JITC_signif], xlab = "t-test p-val", freq = F, 
     main = "t-test p-value histogram \n(include all patients for 5680 JITC peptides)")
hist(Wilcox_pval[JITC_signif], xlab = "Wilcox p-val", freq = F, 
     main = "Wilcoxon-test p-value histogram \n(include all patients for 5680 JITC peptides)")
par(mfrow=c(1,1))

BHpval_JITC <- p.adjust(pval[JITC_signif], method = "BH") # t-test BH-adj p-val
Wilcox_BHpval_JITC <- p.adjust(Wilcox_pval[JITC_signif], method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_JITC, seq(.1, .5,by=.1))
count.func(Wilcox_BHpval_JITC, seq(.1, .5,by=.1))

# exclude 3-month patients
par(mfrow=c(1,2))
hist(pval_drop[JITC_signif], xlab = "t-test p-val", freq = F, 
     main = "t-test p-value histogram \n(Exclude 3-month patients for 5680 JITC peptides)")
hist(Wilcox_pval_drop[JITC_signif], xlab = "Wilcox p-val", freq = F, 
     main = "Wilcoxon-test p-value histogram \n(Exclude 3-month patients for 5680 JITC peptides)")
par(mfrow=c(1,1))

BHpval_drop_JITC <- p.adjust(pval_drop[JITC_signif], method = "BH") # t-test BH-adj p-val
Wilcox_BHpval_drop_JITC <- p.adjust(Wilcox_pval_drop[JITC_signif], method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_drop_JITC, seq(.1, .5,by=.1))
count.func(Wilcox_BHpval_drop_JITC, seq(.1, .5,by=.1))
