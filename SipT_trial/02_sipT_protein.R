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

# sum fluorescence according to protein
dat2 <- dat[,3:ncol(dat)]
dat2 <- aggregate(.~SEQ_ID,dat2,sum)
as.numeric(dat2[dat2$SEQ_ID=="10_RPS28_6234",2:ncol(dat2)]) == 
  as.numeric(colSums(dat[dat$SEQ_ID=="10_RPS28_6234",4:ncol(dat)]))

# check
sum(as.numeric( colnames(dat2[2:ncol(dat2)]) == colnames(dat[4:ncol(dat)]) )) == nrow(samp)
nrow(dat2) == length(unique(dat$SEQ_ID))

# take log2 then find take (post - pre)
dat3 = as.matrix(log2( dat2[,2:ncol(dat2)] ))
dat3 = dat3[,(ncol(dat3)/2 + 1):ncol(dat3)] - dat3[,1:(ncol(dat3)/2)] 
colnames(dat3) <- samp$Patient_ID[1:ncol(dat3)]

# check
ncol(dat3) == length(unique(colnames(dat3)))
samp$Period[1:ncol(dat3)]
samp$Group[1:ncol(dat3)]

# colors and shapes for the visualization techniques
cols = pal[ match(samp$Group[samp$Time == "Post"], names(pal)) ]
shapes = shp[ match(samp$Group[samp$Time == "Post"], names(shp)) ]

####################################################################################### 
#                      2-sample t-test & Wilcoxon Rank-Sum Test                       #
#                             to Assess Treatment Effect                              #
####################################################################################### 

dat3_SEQ_ID = dat2$SEQ_ID
rm(dat, pal, dat2); gc()

# getting variables ready
gr <- samp$Group[1:(nrow(samp)/2)]
pd <- samp$Period[1:(nrow(samp)/2)]
gr_drop <- gr[pd != "3-month"]
dat3_drop <- dat3[,(pd != "3-month")]

# check
table(gr); table(gr_drop)

#---------------------------------------------------------------------------------------
# t-test

# initialize
pval <- rep(NA, nrow(dat3))
est <- matrix(NA, nrow = nrow(dat3), ncol = 2)
StdErr <- rep(NA, nrow(dat3))

pval_drop <- rep(NA, nrow(dat3))
est_drop <- matrix(NA, nrow = nrow(dat3), ncol = 2)
StdErr_drop <- rep(NA, nrow(dat3))

for(i in 1:nrow(dat3)){
  fit1 <- t.test( as.numeric(dat3[i, ]) ~ gr )
  pval[i] <- fit1$p.value
  est[i,] <- fit1$estimate
  StdErr[i] <- fit1$stderr 
  
  fit_drop <- t.test( as.numeric(dat3_drop[i,]) ~ gr_drop )
  pval_drop[i] <- fit_drop$p.value
  est_drop[i,] <- fit_drop$estimate
  StdErr_drop[i] <- fit_drop$stderr 
  
  if(i %% 1000 == 0 | i == nrow(dat3)){
    print(i)
  }
}

#---------------------------------------------------------------------------------------
# Wilcoxon rank-sum test

# initialize
Wilcox_pval <- rep(NA, nrow(dat3))
Wilcox_pval_drop <- rep(NA, nrow(dat3))

for(i in 1:nrow(dat3)){
  Wilcox_pval[i] <- wilcox.test( as.numeric(dat3[i, ]) ~ gr )$p.value
  Wilcox_pval_drop[i] <- wilcox.test( as.numeric(dat3_drop[i,]) ~ gr_drop )$p.value
  if(i %% 1000 == 0 | i == nrow(dat3)){
    print(i)
  }
}

#---------------------------------------------------------------------------------------

# exclude 3-month patient
BHpval_drop <- p.adjust(pval_drop, method = "BH") # t-test BH-adj p-val
Wilcox_BHpval_drop <- p.adjust(Wilcox_pval_drop, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_drop, seq(.05, .3,by=.05))
count.func(Wilcox_BHpval_drop, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(pval_drop, xlab = "t-test p-val", freq = F,  breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(Exclude 3-month patients)")
hist(Wilcox_pval_drop, xlab = "Wilcox p-val", freq = F,  breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(Exclude 3-month patients)")


# get BH-adjusted p-val
BHpval <- p.adjust(pval, method = "BH") # t-test BH-adj p-val
Wilcox_BHpval <- p.adjust(Wilcox_pval, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval, seq(.05, .3,by=.05))
count.func(Wilcox_BHpval, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(pval, xlab = "t-test p-val", freq = F, breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(include all patients)")
hist(Wilcox_pval, xlab = "Wilcox p-val", freq = F, breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(include all patients)")
par(mfrow=c(1,1))

# no peptide is significant after BH adjustment (for both t-test and Wilcoxon test)...

####################################################################################### 
#                          Compare against JITC List                                  #
####################################################################################### 
# get data frame ready

# JITC list
JITC_list <- read_csv("JITC_PAP_Longitudinal.csv")
JITC_signif_protein <- unique(JITC_list$SEQ_ID) 

## 1051 proteins signif out of 1611 --> kinda pointless?

# plot these signif proteins on a volcano plot
signif_ind <- which(dat3_SEQ_ID %in% JITC_signif_protein)

plot(x = est[,1] - est[,2], y = -log10(pval), pch = 20, cex = .3,
     xlab = "fold changes mu_A - mu-B",
     ylab = "-log10(t-test p-values)", 
     main = paste0("volcano plot (", nrow(dat3), " protein) \n comparing group A vs group B") )
lines(x = (est[,1] - est[,2])[signif_ind], y = (-log10(pval))[signif_ind], 
      type = "p", pch = 20, cex = .3, col = "red")


####################################################################################### 
#                      test if mu_A != 0 among group A patients                        #
####################################################################################### 

dat3_grpA <- dat3[,(gr=="A")]
dat3_drop_grpA <- dat3_drop[, (gr_drop == "A")] 

ttest_muA <- rep(NA, nrow(dat3_grpA))
est_muA <- rep(NA, nrow(dat3_grpA))
stderr_muA <- rep(NA, nrow(dat3_grpA))
signedrank_muA <- rep(NA, nrow(dat3_grpA))
ttest_muA_drop <- rep(NA, nrow(dat3_grpA))
est_muA_drop <- rep(NA, nrow(dat3_grpA))
signedrank_muA_drop <- rep(NA, nrow(dat3_grpA))

for(i in 1:nrow(dat3_grpA)){
  fit1 <- t.test(as.numeric(dat3_grpA[i, ]),alternative = "two.sided")
  ttest_muA[i] <- fit1$p.value
  est_muA[i] <- fit1$estimate
  stderr_muA[i] <- fit1$stderr
  signedrank_muA[i] <- wilcox.test(as.numeric(dat3_grpA[i, ]),
                                   alternative = "two.sided",mu=0)$p.value
  
  fit_drop <- t.test(as.numeric(dat3_drop_grpA[i, ]),alternative = "two.sided")
  ttest_muA_drop[i] <- fit_drop$p.value
  est_muA_drop[i] <- fit_drop$estimate
  signedrank_muA_drop[i] <- wilcox.test(as.numeric(dat3_drop_grpA[i, ]),
                                        alternative = "two.sided",mu=0)$p.value
  
  if(i %% 1000 == 0 | i == nrow(dat3_grpA)){
    print(i)
  }
}

# Exclude 3-month patient
BHpval_ttest_muA_drop <- p.adjust(ttest_muA_drop, method = "BH") # t-test BH-adj p-val
BHpval_signedrank_muA_drop <- p.adjust(signedrank_muA_drop, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_ttest_muA_drop, seq(.05, .3,by=.05))
count.func(BHpval_signedrank_muA_drop, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(ttest_muA_drop, xlab = "t-test p-val", freq = F, breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(Exclude 3-month patients)")
hist(signedrank_muA_drop, xlab = "Wilcox p-val", freq = F, breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(Exclude 3-month patients)")

# Include all patients
BHpval_ttest_muA <- p.adjust(ttest_muA, method = "BH") # t-test BH-adj p-val
BHpval_signedrank_muA <- p.adjust(signedrank_muA, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_ttest_muA, seq(.05, .3,by=.05))
count.func(BHpval_signedrank_muA, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(ttest_muA, xlab = "t-test p-val", freq = F, breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(Include all patients)")
hist(signedrank_muA, xlab = "Wilcox p-val", freq = F, breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(Include all patients)")

####################################################################################### 
#                      test if mu != 0 for ALL patients grouped together               #
#######################################################################################

ttest_all <- rep(NA, nrow(dat3))
est_all <- rep(NA, nrow(dat3))
stderr_all <- rep(NA, nrow(dat3))
signedrank_all <- rep(NA, nrow(dat3))
ttest_all_drop <- rep(NA, nrow(dat3_drop))
est_all_drop <- rep(NA, nrow(dat3_drop))
signedrank_all_drop <- rep(NA, nrow(dat3_drop))

for(i in 1:nrow(dat3)){
  y <- as.numeric(dat3[i, ])
  fit1 <- t.test(y,alternative = "two.sided")
  ttest_all[i] <- fit1$p.value
  est_all[i] <- fit1$estimate
  stderr_all[i] <- fit1$stderr
  signedrank_all[i] <- wilcox.test(y, alternative = "two.sided",mu=0)$p.value
  
  fit_drop <- t.test(as.numeric(dat3_drop[i, ]),alternative = "two.sided")
  ttest_all_drop[i] <- fit_drop$p.value
  est_all_drop[i] <- fit_drop$estimate
  signedrank_all_drop[i] <- wilcox.test(as.numeric(dat3_drop[i, ]), alternative = "two.sided",mu=0)$p.value
  
  if(i %% 1000 == 0 | i == nrow(dat3)){
    print(i)
  }
}

# Exclude 3-month patient
BHpval_ttest_all_drop <- p.adjust(ttest_all_drop, method = "BH") # t-test BH-adj p-val
BHpval_signedrank_all_drop <- p.adjust(signedrank_all_drop, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_ttest_all_drop, seq(.05, .3,by=.05))
count.func(BHpval_signedrank_all_drop, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(ttest_all_drop, xlab = "t-test p-val", freq = F, breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(Exclude 3-month patients)")
hist(signedrank_all_drop, xlab = "Wilcox p-val", freq = F, breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(Exclude 3-month patients)")


# Include all patients
BHpval_ttest_all <- p.adjust(ttest_all, method = "BH") # t-test BH-adj p-val
BHpval_signedrank_all <- p.adjust(signedrank_all, method = "BH") # wilcoxon BH-adj p-val
count.func(BHpval_ttest_all, seq(.05, .3,by=.05))
count.func(BHpval_signedrank_all, seq(.05, .3,by=.05))
par(mfrow=c(1,2))
hist(ttest_all, xlab = "t-test p-val", freq = F, breaks= 50, col = "gray",
     main = "t-test p-value histogram \n(Include all patients)")
hist(signedrank_all, xlab = "Wilcox p-val", freq = F, breaks= 50, col = "gray",
     main = "Wilcoxon-test p-value histogram \n(Include all patients)")


which(BHpval_ttest_all <= .05) == which(BHpval_ttest_all_drop <= .2)
dat3_SEQ_ID[which(BHpval_ttest_all <= .05)]
dat3_SEQ_ID[which(BHpval_ttest_all <= .05)] %in% JITC_signif_protein
dat3_SEQ_ID[which(BHpval_signedrank_all <= .05)]
dat3_SEQ_ID[which(BHpval_signedrank_all <= .05)] %in% JITC_signif_protein

save(pval, est, StdErr, pval_drop, est_drop, StdErr_drop, Wilcox_pval, Wilcox_pval_drop,
     ttest_muA, est_muA, stderr_muA, signedrank_muA, ttest_muA_drop, est_muA_drop, signedrank_muA_drop,
     ttest_all, est_all, stderr_all, signedrank_all, ttest_all_drop, est_all_drop, signedrank_all_drop,
     dat3_SEQ_ID,
     file = "sipT_protein_results.RData")
