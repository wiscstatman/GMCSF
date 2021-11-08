library(Rtsne)
library(fdrtool)
library(dendextend) # dendrogram
library(ggplot2)
library(ashr)
library(tidyverse)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

# specified color and shape scheme 
pal <- c("navy", "darkorange1")
names(pal) <- c("GM",  "pTVG-HP")
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
  countab <- rbind(c("Threshold", thresh.vec),
                   c("Peptide counts", counter))
  return(countab)
}

# typical step in analysis
analysis_func <- function(pval, xlab){
  # control FDR
  anova_BH <- p.adjust(pval, method = "BH")
  anova_qval <- fdrtool(pval, statistic = "pvalue", verbose = F, plot  = F)$qval
  anova_qval_eta0 <- unname(fdrtool(pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
  
  # plot histogram of p-values
  all_peptide_hist <- hist(pval, breaks = 70, freq = F, xlab = xlab, las = 1,  
                           main = paste0("Histogram of ", xlab))
  # polygon_ind <- which(all_peptide_hist$density >= anova_qval_eta0)
  # for (i in polygon_ind){
  #   polygon( x = c(all_peptide_hist$breaks[i], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i+1], all_peptide_hist$breaks[i]),
  #            y = c(all_peptide_hist$density[i], all_peptide_hist$density[i], anova_qval_eta0, anova_qval_eta0),
  #            col = "red")
  # }
  # text(x=0.65,y=4, labels = paste( "estimated proportion of \nnon-null peptides =",
  #                                  round( 100*(1 - anova_qval_eta0),2 ),"%" ))
  
  return(list(BH = anova_BH, qval = anova_qval, qval_eta0 = anova_qval_eta0))
}

# post-process allez table for kable output
get_alleztable.func <- function(allez_go_input){
  allez.tab <- allezTable(allez_go_input, symbol = T, nominal.alpha = nom.alpha, n.upp = max_gene_in_set, in.set = T)
  allez.tab$set.size <- paste(allez.tab$in.set, allez.tab$set.size, sep = "/")
  allez.tab <- allez.tab %>% dplyr::select(-c(in.set, genes)) %>%
    mutate(in.genes = str_replace_all(in.genes, ";", "; "))
  return(allez.tab)
}


####################################################################################### 
#                                      Setup Dataset                                  #
####################################################################################### 

dat <- read_tsv("Aggregate_processed_data/Processed_aggregate_data.txt")
colnames(dat)

# unify column names
colnames(dat)[30:33] <- paste0( substr(colnames(dat)[30:33],1,10), substr(colnames(dat)[30:33],18,29) )
colnames(dat)[94:133] <- paste0( substr(colnames(dat)[94:133],1,10), substr(colnames(dat)[94:133],18,29) )
colnames(dat) # check
length(unique(colnames(dat)[18:207])) # check

# remove redundant substring of colnames
table(substr(colnames(dat)[18:207],8,15)) # check
colnames(dat)[18:207] <- paste0( substr(colnames(dat)[18:207],1,7), substr(colnames(dat)[18:207],16,22) )
colnames(dat) # check
length(unique(colnames(dat)[18:207])) # check

#--------------------------------------------------------------------------------------
# we have duplicates for every peptide 
# take (double) log2 transformation
# collapse rows by mean (aka median since only 2 duplicates)

# check unique probe_id 
length(unique(dat$PROBE_ID)) # check
table(as.data.frame(table(dat$PROBE_ID))$Freq) # check

# take (double) log2 of dat
dat_log2 <- dat
dat_log2 <- dat_log2 %>% 
  mutate_at(vars(contains(".dat")), log2)
  # mutate_at(vars(contains(".dat")), function(x){log2(log2(x))})

# arrange rows by PROBE_ID then ROW_NUM then COL_NUM
dat_log2 <- dat_log2 %>% 
  arrange(PROBE_ID, ROW_NUM, COL_NUM) 

# check
dat_log2[1:8, c("PROBE_ID", "ROW_NUM", "COL_NUM")]
dat_log2 %>% 
  select(contains(".dat")) %>%
  as.matrix() %>% as.numeric() %>%
  hist() # summary() 

# extract odd rows and even rows
dat_odd_rows <- dat_log2[seq(1, nrow(dat_log2)-1, by = 2), ]
dat_even_rows <- dat_log2[seq(2, nrow(dat_log2), by = 2), ]

# check both have aligned rows and columns
dat_odd_rows[1:2, 15:25]
dat_even_rows[1:2, 15:25]
length(unique(dat$PROBE_ID)) == length(unique(dat_odd_rows$PROBE_ID))
length(unique(dat$PROBE_ID)) == length(unique(dat_even_rows$PROBE_ID))
sum(as.numeric(dat_odd_rows$PROBE_ID == dat_even_rows$PROBE_ID)) == nrow(dat)/2
sum(as.numeric(colnames(dat_odd_rows) == colnames(dat_even_rows))) == ncol(dat)

# mean 
collapsed_dat <- (dat_odd_rows[, 18:207] + dat_even_rows[, 18:207])/2
## hard-coded column numbers
## do NOT use dplyr's "contains()" because different version of package may mess up sequence of columns

# get unique probe ID
collapsed_dat <- cbind(dat_odd_rows[,"PROBE_ID"], collapsed_dat)
colnames(collapsed_dat) # check

# remove unnecessary data frames to free up memory
rm(dat_odd_rows, dat_even_rows, dat_log2); gc()

#--------------------------------------------------------------------------------------
# desired columns of ref_data_frame: patient ID, before/after, array ID
# desired rows of ref_data_frame: arranged according to before/after, then patient ID
# arrange the columns of collapsed_dat to the rows of this desired ref_data_frame

sample_dat <- read_csv("sample_data.csv")

# exclude PAP099, PAP124, PAP071, PAP096, PAP092, and PAP118
to_exclude <- sample_dat[sample_dat$"Sample_ID" %in% 
                           c("PAP099", "PAP124", "PAP071", "PAP096", "PAP092", "PAP118"),]
collapsed_dat <- collapsed_dat %>%
  select(- matches(to_exclude$Colnames))
sample_dat <- sample_dat %>%
  filter(Patient_ID != "ignore")
rm(to_exclude); gc()

# check
sum(as.numeric( colnames(collapsed_dat %>% select(contains(".dat"))) %in% sample_dat$Colnames )) == nrow(sample_dat)
sum(as.numeric( sample_dat$Colnames %in% colnames(collapsed_dat %>% select(contains(".dat"))) )) == nrow(sample_dat)
length(unique( colnames(collapsed_dat %>% select(contains(".dat"))) )) == nrow(sample_dat)
length(unique( sample_dat$Colnames )) == nrow(sample_dat)

# arrange according to before/after, then patient ID
sample_dat <- sample_dat %>%
  arrange(Before_After, Patient_ID)

# check
sample_dat 
sum(as.numeric( sample_dat$Patient_ID[sample_dat$Before_After == "After"] == 
                  sample_dat$Patient_ID[sample_dat$Before_After == "Before"])) == nrow(sample_dat)/2
sample_dat$Colnames == colnames( collapsed_dat[, ! (names(collapsed_dat) %in% "PROBE_ID")] ) 
## note that colnames not matched

# arrange columns of collapsed_dat 
iii <- match( sample_dat$Colnames, colnames( collapsed_dat[, ! (names(collapsed_dat) %in% "PROBE_ID")] ) )
collapsed_dat <- collapsed_dat[, c(1, iii+1)]
sample_dat$Colnames == colnames( collapsed_dat[, 2: ncol(collapsed_dat)] ) 
## note that colnames are matched now

#--------------------------------------------------------------------------------------
# compute: After - Before 

after_ind <- (1 : (nrow(sample_dat)/2))
before_ind <- ((nrow(sample_dat)/2 + 1) : nrow(sample_dat))
sample_dat$Before_After[after_ind] # check
sample_dat$Before_After[before_ind] # check

after_dat <- collapsed_dat[, after_ind + 1]
before_dat <- collapsed_dat[, before_ind + 1]
diff_mat <- as.matrix(after_dat - before_dat)
colnames(diff_mat) <- sample_dat$Patient_ID[sample_dat$Before_After == "After"]
row.names(diff_mat)<- collapsed_dat$PROBE_ID
head(diff_mat) # check

# remove unnecessary objects
rm(after_ind, before_ind, after_dat, before_dat, iii); gc()

# colors and shapes for the visualization techniques
cols = pal[ match(sample_dat$Treatment[sample_dat$Before_After == "After"], names(pal)) ]
shapes = shp[ match(sample_dat$Treatment[sample_dat$Before_After == "After"], names(shp)) ]

####################################################################################### 
#                               Preliminary Analysis -- Clustering                    #
####################################################################################### 

clust.dat <- t(diff_mat)

# hc_complete <- hclust(dist(clust.dat), "complete")
hc_average <- hclust(dist(clust.dat), "average")

# dend <- as.dendrogram(hc_complete)
dend <- as.dendrogram(hc_average)

# coloring
labels_colors(dend) <- cols[order.dendrogram(dend)]

# tiff("hclust_CompleteLinkage_log2log2.tiff", res=600, compression = "lzw", width=8, height=3.5, units="in", pointsize = 5)
tiff("hclust_AverageLinkage_log2log2.tiff", res=600, compression = "lzw", width=8, height=3.5, units="in", pointsize = 5)
plot(dend, xlab = "Patients", ylab = "Dissimilarity")
dev.off()

####################################################################################### 
#                               Preliminary Analysis -- PCA                           #
####################################################################################### 

# svd
sv.dat <- sweep(t(diff_mat), 2, colMeans(t(diff_mat)), "-") # centering
sv <- svd(sv.dat)
V <- sv$v
D <- sv$d
U <- sv$u

# variance explained
pca.var <- D^2/sum(D^2) 
pca.cumvar <- cumsum(pca.var)

# plot PCA
tiff("PCA_log2log2.tiff", res=600, compression = "lzw", width=6, height=3.5, units="in")
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.6,0.4,0),cex=.5,
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func(cols, shapes, U, D, 1, 2, pca.var, title = "PC2 vs PC1") # PC loadings (PC2 vs PC1)
legend('topright', pch = shp, col = pal, c("GM",  "pTVG-HP") )
PCload.func(cols, shapes, U, D, 3, 2, pca.var, title = "PC2 vs PC3") # PC loadings (PC2 vs PC3)
dev.off()

# zoom-in 
PCload.func2 <- function(cols, shapes, U, D, x, y, pca.vec, title, xlim, ylim){
  Z <- U %*% diag(D)
  plot(Z[,x], Z[,y], col = cols, pch = shapes, main = title, las = 1, xlim = xlim, ylim = ylim,
       xlab = paste0("PC",x," (", round(pca.vec[x]*100, 1) ,"% variance explained)"),
       ylab = paste0("PC",y," (", round(pca.vec[y]*100, 1) ,"% variance explained)") 
  )
}
tiff("PCA_zoomin_log2log2.tiff", res=600, compression = "lzw", width=6, height=3.5, units="in")
par(mfrow = c(1,2), pty = "s", mar = c(2.2,2.3,1.5,0.45), mgp = c(1.4,0.4,0), cex=.5,
    cex.axis = 0.84, cex.lab = 0.84, cex.main = 0.84, tcl = -0.4)
PCload.func2(cols, shapes, U, D, 1, 2, pca.var, title = "PC1 vs PC2", xlim = c(-35,30), ylim = c(-75,65)) 
legend('topright', pch = shp, col = pal, c("GM",  "pTVG-HP") )
PCload.func2(cols, shapes, U, D, 3, 2, pca.var, title = "PC3 vs PC2", xlim = c(-10,20), ylim = c(-75,65)) 
dev.off()

####################################################################################### 
#                              Preliminary Analysis -- t-SNE                          #
####################################################################################### 

# how to specify parameter
# refer: https://lvdmaaten.github.io/tsne/User_guide.pdf
tsnedat <- unname(t(diff_mat)) # dim(X) = N samples by D dimensions 
initdim <- 90 
perplex <- 30 
lab <- row.names(t(diff_mat))

set.seed(10)
tsne_raw <- Rtsne(tsnedat, dims = 2L, initial_dims = initdim, perplexity = perplex,
                  theta = 0, check_duplicates = F, PCA = T)

# t-SNE plot
plot(tsne_raw$Y, ylab = "", xlab = "", main = "tSNE plot", col = cols, pch = shapes)
legend('topleft', pch = shp, col = pal, cex = 0.5, c("GM",  "pTVG-HP") )
dev.off()

####################################################################################### 
#                   T-test & Wilcoxon Rank-Sum Test & ASH                             #
####################################################################################### 

sum(as.numeric( sample_dat$Patient_ID[sample_dat$Before_After == "After"] == colnames(diff_mat) )) == ncol(diff_mat) #check

# get treatment label
treat <- as.factor( sample_dat$Treatment[sample_dat$Before_After == "After"] )
table(treat) # check

# initiate
t_test_pval <- rep(NA, nrow(diff_mat))
wilcox_pval <- rep(NA, nrow(diff_mat))
names(t_test_pval) <- row.names(diff_mat)
names(wilcox_pval) <- row.names(diff_mat)
ash_betahat <- rep(NA, nrow(diff_mat))
ash_sebetahat <- rep(NA, nrow(diff_mat))

# compute t-tests and wilcoxon rank-sum tests p-values
for(i in 1:nrow(diff_mat)){
  t_test_fit <- t.test(diff_mat[i,] ~ treat)
  t_test_pval[i] <- t_test_fit$p.value
  wilcox_pval[i] <- wilcox.test(diff_mat[i,] ~ treat)$p.value
  ash_betahat[i] <- unname(diff(t_test_fit$estimate))
  ash_sebetahat[i] <- t.test(diff_mat[i,] ~ treat, var.equal = T)$stderr

  if( i %% 100 == 0 | i == nrow(diff_mat)){
    print(i)
  }
}

ash_lfsr <- ash( betahat = ash_betahat, sebetahat = ash_sebetahat, df = 90)
ash_lfsr <- ash_lfsr$result$lfsr
  
# save result
# save(t_test_pval, wilcox_pval, ash_lfsr, file = "t_test_wilcox_test_pvals_log2.RData")
# save(t_test_pval, wilcox_pval, ash_lfsr, file = "t_test_wilcox_test_pvals_log2log2.RData")

# load either one
# load("t_test_wilcox_test_pvals_log2.RData")
# load("t_test_wilcox_test_pvals_log2log2.RData")

t_test_analysis <- analysis_func(t_test_pval, "t-test p-values")
wilcox_analysis <- analysis_func(wilcox_pval, "wilcoxon rank-sum test p-values")
hist(ash_lfsr, breaks = 50, xlab = "local false sign rates (LFSR)", 
     main = "Histogram of local false sign rates (LFSR) \nfrom Adaptive Shrinkage (ASH) method")

# peptide counts at various FDR thresholds
count.func(t_test_analysis$BH, seq(0.05, 0.5, by = 0.05))
count.func(wilcox_analysis$BH, seq(0.05, 0.5, by = 0.05))
count.func(ash_lfsr, seq(0.05, 0.5, by = 0.05))

# compare wilcoxon pval vs t-test pval
plot(x = t_test_pval, y = wilcox_pval, pch = ".", xlab = "t-test p-values", 
     ylab = "wilcoxon rank-sum test p-values")
# zoom-in
plot(x = t_test_pval, y = wilcox_pval, xlab = "t-test p-values", pch = ".",
     ylab = "wilcoxon rank-sum test p-values", xlim = c(0,.05), ylim = c(0,.05))


####################################################################################### 
#                        Get fold-change for all  peptides                            #
####################################################################################### 

load("t_test_wilcox_test_pvals_log2.RData")
# 
# t_test_wilcox <- (t_test_pval <= .05 | wilcox_pval <= .05)
# 
# eff <- t( apply( diff_mat[t_test_wilcox,], 1, function(x){
#   t.test(x ~ treat)$estimate
# } ) )

eff <- matrix(NA, nrow = nrow(diff_mat), ncol = 2)
for(i in 1:nrow(diff_mat)){
  eff[i,] <- t.test(diff_mat[i,] ~ treat)$estimate
}
eff <- round(eff, 4)

# eff <- data.frame(PROBE_ID = row.names(eff), eff)
# row.names(eff) <- NULL
eff <- data.frame(PROBE_ID = row.names(diff_mat), eff)

colnames(eff)[2:3] <- c("mean_in_group_GM", "mean_in_group_pTVG-HP")
eff$Difference <- eff$`mean_in_group_pTVG-HP` - eff$mean_in_group_GM 
eff <- eff %>% arrange(desc(Difference))

to_merge <- dat[dat$PROBE_ID %in% eff$PROBE_ID, c("PROBE_ID", "SEQ_ID")]
to_merge <- distinct(to_merge)
seq_iii <- match(eff$PROBE_ID, to_merge$PROBE_ID)
sum(as.numeric( eff$PROBE_ID == to_merge$PROBE_ID[seq_iii] )) == nrow(eff) # check

tab <- eff
tab$SEQ_ID <- to_merge$SEQ_ID[seq_iii]
tab <- tab %>% 
  select(PROBE_ID, SEQ_ID, everything())

tab$T_test_pval <- round(t_test_pval,4)
tab$Wilcox_pval <- round(wilcox_pval,4)

tab$overall_mean <- round( (47*tab$`mean_in_group_pTVG-HP`+45*tab$mean_in_group_GM)/92 , 4)

tab <- tab %>% arrange(desc(overall_mean))

# JITC list
JITC_list <- read_csv("JITC_PAP_Longitudinal.csv")
tab$JITC_signif <- as.numeric(tab$PROBE_ID %in% JITC_list$PROBE_ID)

# rearrange diff_mat rows according to tab
tab_iii <- match(tab$PROBE_ID, row.names(diff_mat))
diff_mat2 <- diff_mat[tab_iii,]
sum(as.numeric(row.names(diff_mat2) == tab$PROBE_ID)) 

# get indices of the 2 groups
ind_pTVG_HP <- which(sample_dat$Treatment[sample_dat$Before_After=="After"] == "pTVG-HP")
ind_GM <- which(sample_dat$Treatment[sample_dat$Before_After=="After"] == "GM")

stdev_pTVG_HP <- apply( diff_mat2[,ind_pTVG_HP], 1, function(x){sd(x)})
stdev_GM <- apply( diff_mat2[,ind_GM], 1, function(x){sd(x)})

tab <- cbind(tab, stdev_pTVG_HP, stdev_GM)
tab$stdev_pTVG_HP <- round(tab$stdev_pTVG_HP, 4)
tab$stdev_GM <- round(tab$stdev_GM, 4)

tab$studentized_pTVG_HP <- round( tab$`mean_in_group_pTVG-HP`/ tab$stdev_pTVG_HP, 4)
tab$studentized_GM <- round( tab$mean_in_group_GM / tab$stdev_GM, 4)

tab$overall_stdev <- round( apply( diff_mat2, 1, function(x){sd(x)}), 4)
tab$overall_studentized <- round( tab$overall_mean / tab$overall_stdev, 4)

tab <- tab %>% arrange(desc(overall_studentized))

# write.table(tab, "top_fold_change.csv", sep = ",", col.names = T, row.names = F)
write.table(tab, "fold_change_all_peptides.csv", sep = ",", col.names = T, row.names = F)


####################################################################################### 
#                          Compare against JITC List                                  #
####################################################################################### 


compare_hist_df <- data.frame(
  mean_fold_changes = tab$overall_mean,
  JITC_signif = as.factor(tab$JITC_signif)
)

ggplot(compare_hist_df, aes(x = mean_fold_changes, fill = JITC_signif)) +
  geom_histogram(aes(y=..density..), bins = 100, position = "identity", alpha = .4) +
  labs(x = "mean fold changes", 
       title = paste0("Density Histograms of Mean Fold Changes")) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")


compare_hist_df2 <- data.frame(
  mean_fold_changes = tab$`mean_in_group_pTVG-HP`,
  JITC_signif = as.factor(tab$JITC_signif)
)


# tab_JITC <- tab[ tab$PROBE_ID %in% JITC_list$PROBE_ID ,]
# 
# summary(tab$Difference)
# 
# # get histogram of difference of fold changes
# 
# hist(tab_JITC$Difference, main = "Histogram of Differences of Fold Changes
#      Between group_pTVG-HP and group_GM
#      for the 5680 Signif Peptides in JITC List", breaks = 70,
#      xlab = "Differences of Fold Changes")
# 
# # plot t-test pval vs wilcox pval, colored by differences of fold changes
# ggplot(tab_JITC, aes(x = T_test_pval, y = Wilcox_pval, 
#                      color = pmin(pmax(Difference,-.5),.7) )) + 
#   geom_point() +
#   scale_color_gradient2(midpoint = median(tab_JITC$Difference), low="blue", mid = "white", 
#                         high="red", name = "Difference of Fold Changes") +
#   labs(title = "T-test p-values vs Wilcoxon p-values colored by Difference of Fold Changes 
#        for the 5680 Significant Peptides in JITC List") +
#   theme(plot.title = element_text(hjust = 0.5))
