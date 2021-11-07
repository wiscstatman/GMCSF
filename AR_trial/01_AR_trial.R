library(lmerTest)
library(lme4)
library(locfdr)
library(tidyverse)
library(RColorBrewer)
library(fdrtool)

####################################################################################### 
#                           Some Common Variables and Functions                       #
####################################################################################### 

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

load("01_ARtrial_results.RData")

####################################################################################### 
#                                      Setup Dataset                                  #
####################################################################################### 

dat <- read_tsv("Aggregate_processed_data/Processed_aggregate_data.txt")
colnames(dat)

length(unique(colnames(dat %>% select(contains(".dat"))))) #check
colnames(dat)[18:197] <- paste0( substr(colnames(dat)[18:197],1,7), substr(colnames(dat)[18:197],16,18) )
colnames(dat) #check

summary(as.numeric(as.matrix(log2(dat[,18:197])))) # check
dat_log2 <- dat
dat_log2[,18:197] <- log2(dat[,18:197])

#--------------------------------------------------------------------------------------

sample_dat <- read_csv("PA_00166_Experiment_Layout.csv")
sum(as.numeric(sample_dat$Array_ID %in% colnames(dat_log2)[18:197])) # check
sum(as.numeric(colnames(dat_log2)[18:197] %in% sample_dat$Array_ID )) # check

# arrange by time, then treatment type, then patient_ID
sample_dat <- sample_dat %>% 
  arrange(factor(Time, levels = c("Pre", "3_month", "6_month")), 
          factor(Treatment, levels = c("VX1_no", "VX1_GM", "VX2_no", "VX2_GM")), 
          Patient_ID)

# remove not-relevant columns in dat_log2?
dat_log2 <- cbind( dat_log2 %>% select(PROBE_DESIGN_ID:Y),
                   subset(dat_log2[, 18:197], 
                          select = colnames(dat_log2)[18:197] %in% sample_dat$Array_ID ))
colnames(dat_log2) # check

# arrange columns of dat_log2 
iii <- match( sample_dat$Array_ID, colnames(dat_log2[,18:ncol(dat_log2)]) )
dat_log2 <- dat_log2[, c(1:17, iii+17)]

# check
colnames(dat_log2) 
sum(as.numeric( sample_dat$Array_ID == colnames( dat_log2[, 18: ncol(dat_log2)] )  )) == nrow(sample_dat)
hist(as.numeric(as.matrix(dat_log2[,18:ncol(dat_log2)])))

# patient_time tag
sample_dat <- sample_dat %>%
  mutate(Time2 = as_factor(Time),
         Time2 = fct_recode(Time2,
           "Pre" = "Pre",
            "3mo" = "3_month",
            "6mo" = "6_month" 
         ))
sample_dat$Patient_Time <- paste(sample_dat$Patient_ID, sample_dat$Time2, sep = "_")
colnames(dat_log2)[18: ncol(dat_log2)] <- sample_dat$Patient_Time 
colnames(dat_log2)  # check

# separate treatment type and schedule
sample_dat$Treat <- factor(sub(".*_","",sample_dat$Treatment), levels = c("no", "GM"))
sample_dat$Schedule <- factor(substr(sample_dat$Treatment, 1, 3), levels = c("VX1", "VX2") )

# set time as numeric variable
time_vec <- c(0,3,6)
names(time_vec) <- levels(sample_dat$Time2)
sample_dat$Time3 <- unname(time_vec[sample_dat$Time2])
rm(time_vec); gc()
  
#--------------------------------------------------------------------------------------
# set color

sample_dat <- sample_dat %>%
  mutate(Treatment = factor(Treatment, levels = c("VX1_no", "VX1_GM", "VX2_no", "VX2_GM")),
         Time = factor(Time, levels = c("Pre", "3_month", "6_month")))

pal <- brewer.pal(9,"Set1")[c(6,5,4,7)]
names(pal) <- levels(sample_dat$Treatment)


####################################################################################### 
#                           Check Fluorescence Normalization                          #
####################################################################################### 

dat_log2_long <- dat_log2[, 18:ncol(dat_log2)] %>% 
  pivot_longer(cols = everything(), names_to = "id_time", values_to = "fluorescence") 

# check
nrow(dat_log2_long) == nrow(dat_log2) * nrow(sample_dat)
head(dat_log2_long)

# set fill color
dat_log2_long$treat <- rep(sample_dat$Treatment, 177604)

# sort order of patients in boxplot
dat_log2_long$id_time <- factor(dat_log2_long$id_time, levels = sample_dat$Patient_Time)

tiff("fluorescence_normalization.tiff", res=600, compression = "lzw", width=15, height=10, units="in", pointsize = 10)
ggplot(dat_log2_long, aes(x = id_time, y = fluorescence, fill = treat)) +
  geom_boxplot(outlier.shape = ".") +
  scale_fill_manual(name = "treatment", values = pal[levels(dat_log2_long$treat)]) +
  labs(title = "Boxplots of Peptide Fluorescence Levels for Patients at 3 time points", 
       x = "Patient_Time", y = "Log2 Fluorescence Levels") +
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))
dev.off()

####################################################################################### 
#                    Linear Mixed Model to Assess Time Effect (REML = TRUE)           #
####################################################################################### 

ncol_dat <- ncol(dat_log2)

Test_Time.func <- function(y, treat_type){
  resp <- y[sample_dat$Treat == treat_type]
  fit1 <- lmer(resp ~ Time3 + (1 + Time3 | Patient_ID), REML = T,
               data = sample_dat[sample_dat$Treat == treat_type,])
  fit0 <- lmer(resp ~ 1 + (1 + Time3 | Patient_ID), REML = T,
               data = sample_dat[sample_dat$Treat == treat_type,])
  resid_fit0 <- unname(round(resid(fit0),4))
  effect_tstat <- coef(summary(fit1))['Time3',c('Estimate', 't value')]
  KR_df_pval <- contest(fit1, c(0,1), ddf = "Kenward-Roger")[c('DenDF', 'Pr(>F)')]
  Satterthwaite_df_pval <- contest(fit1, c(0,1))[c('DenDF', 'Pr(>F)')]
  zval_1sided_KR <- qnorm(pt( as.numeric(effect_tstat['t value']),
                              df = as.numeric(KR_df_pval['DenDF']) ,  lower.tail = T ))
  zval_1sided_Satterthwaite <- qnorm(pt( as.numeric(effect_tstat['t value']),
                                         df = as.numeric(Satterthwaite_df_pval['DenDF']) ,  lower.tail = T ))
  result <- c(
    as.numeric(effect_tstat),
    as.numeric(KR_df_pval),
    as.numeric(Satterthwaite_df_pval),
    zval_1sided_KR,
    zval_1sided_Satterthwaite
  )
  return( list(
    resid_fit0 = resid_fit0,
    result = result
  ) )
}

# initiate
GM_resid_fit0 <- matrix(NA, nrow = nrow(dat_log2), ncol = nrow(sample_dat%>%filter(Treat == "GM")))
noGM_resid_fit0 <- matrix(NA, nrow = nrow(dat_log2), ncol = nrow(sample_dat%>%filter(Treat == "no")))
GM_result <- matrix(NA, nrow = nrow(dat_log2), ncol = 8)
noGM_result <- matrix(NA, nrow = nrow(dat_log2), ncol = 8)
colnames(GM_resid_fit0) <- colnames( subset(dat_log2[, 18:ncol_dat], 
                                            select = sample_dat$Treat == "GM" ) )
colnames(noGM_resid_fit0) <- colnames( subset(dat_log2[, 18:ncol_dat], 
                                              select = sample_dat$Treat == "no" ) )
colnames(GM_result) <- paste0("GM_", c(
  "time_effect",
  "time_tstat",
  "KR_df",
  "KR_Ftest_pval",
  "Satterthwaite_df",
  "Satterthwaite_Ftest_pval",
  "zval_1sided_KR",
  "zval_1sided_Satterthwaite"
))
colnames(noGM_result) <- paste0("noGM_", c(
  "time_effect",
  "time_tstat",
  "KR_df",
  "KR_Ftest_pval",
  "Satterthwaite_df",
  "Satterthwaite_Ftest_pval",
  "zval_1sided_KR",
  "zval_1sided_Satterthwaite"
))

for(i in 1:nrow(dat_log2)){
  y <- as.numeric(dat_log2[i, 18:ncol_dat])
  GM_test <- Test_Time.func(y, "GM")
  noGM_test <- Test_Time.func(y, "no")
  
  GM_resid_fit0[i,] <- GM_test$resid_fit0
  GM_result[i,] <- GM_test$result
  
  noGM_resid_fit0[i,] <- noGM_test$resid_fit0
  noGM_result[i,] <- noGM_test$result
  
  if(i %% 100 == 0){
    print(i)
  }
}

# save(GM_resid_fit0, noGM_resid_fit0, GM_result, noGM_result,
#      file = "01_ARtrial_results.RData")

#---------------------------------------------------------------------------------------------
# F-test p-values based on KR adjustments more conservative than Satterthwaite

GM_Satterth_Ftest_pval <- GM_result[,"GM_Satterthwaite_Ftest_pval"]
GM_KR_Ftest_pval <- GM_result[,"GM_KR_Ftest_pval"]
noGM_Satterth_Ftest_pval <- noGM_result[,"noGM_Satterthwaite_Ftest_pval"]
noGM_KR_Ftest_pval <- noGM_result[,"noGM_KR_Ftest_pval"]

GM_Ftest_KR_BH <- p.adjust(GM_KR_Ftest_pval, method="BH")
GM_Ftest_Satterthwaite_BH <- p.adjust(GM_Satterth_Ftest_pval, method="BH")
noGM_Ftest_KR_BH <- p.adjust(noGM_KR_Ftest_pval,method="BH")
noGM_Ftest_Satterthwaite_BH <- p.adjust(noGM_Satterth_Ftest_pval,method="BH")

# fdrtool's eta0
KR_eta0 <- unname(fdrtool(GM_KR_Ftest_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])
Satter_eta0 <- unname(fdrtool(GM_Satterth_Ftest_pval, statistic = "pvalue", verbose = F, plot  = F)$param[,"eta0"])

par(mfrow=c(1,2))
plot(GM_Satterth_Ftest_pval[GM_Satterth_Ftest_pval <= .2 & GM_KR_Ftest_pval <= .2], 
     GM_KR_Ftest_pval[GM_Satterth_Ftest_pval <= .2 & GM_KR_Ftest_pval <= .2], 
     pch = ".", xlim = c(0,.2), ylim = c(0,.2), las = 1, 
     main = "Time Fixed Effect p-values \nfor GM patients",
     xlab = "Satterthwaite F-test p-values", ylab = "Kenward-Roger (KR) F-test p-values")
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

plot(noGM_Satterth_Ftest_pval[noGM_Satterth_Ftest_pval <= .2 & noGM_KR_Ftest_pval <= .2], 
     noGM_KR_Ftest_pval[noGM_Satterth_Ftest_pval <= .2 & noGM_KR_Ftest_pval <= .2], 
     pch = ".", xlim = c(0,.2), ylim = c(0,.2), las = 1, 
     main = "Time Fixed Effect p-values \nfor noGM patients",
     xlab = "Satterthwaite F-test p-values", ylab = "Kenward-Roger (KR) F-test p-values")
abline(a=0, b=1, col = "red", lty=2, lwd = 2)

dev.off()


count.func(GM_Ftest_KR_BH, seq(.2,.6,by=.05))
count.func(GM_Ftest_Satterthwaite_BH, seq(.2,.6,by=.05))
count.func(noGM_Ftest_KR_BH, seq(.2,.6,by=.05))
count.func(noGM_Ftest_Satterthwaite_BH, seq(.2,.6,by=.05))

#---------------------------------------------------------------------------------------------
#p-value density histograms

KR_Ftest_pval_df <- data.frame(
  treatment = c(
    rep("GM", 177604*2),
    rep("noGM", 177604*2)
  ), 
  method = rep( rep( c("Satterthwaite", "Kenward-Roger"), each = 177604 ), 2) ,
  p_values = c(
    GM_Satterth_Ftest_pval,
    GM_KR_Ftest_pval,
    noGM_Satterth_Ftest_pval,
    noGM_KR_Ftest_pval
  )
)

# ggplot(KR_Ftest_pval_df, aes(x = p_values, fill = method)) +
#   geom_histogram(aes(y=..density..), bins = 100, position = "identity", alpha = .4) +
#   facet_grid(. ~ treatment) +
#   labs(x = "F-test p-values", title = paste0("Density Histograms of F-test p-values")) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position = "bottom")

ggplot(KR_Ftest_pval_df[KR_Ftest_pval_df$treatment == "GM",], aes(x = p_values, fill = method)) +
  geom_histogram(aes(y=..density..), bins = 70, position = "identity", alpha = .4) +
  labs(x = "F-test p-values", title = "Density histograms of F-test p-values for GM patients") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") + 
  annotate("text", x = .5, y=1.8, label = paste0("Estimated proportion of non-null peptides = ", 
                                                  round(1 - KR_eta0,4)*100 , 
                                                 "% based on KR F-test p-values" )) 
  # + annotate("text", x = .5, y=1.5, label = paste0("Estimated proportion of non-null peptides = ", 
  #                                                round(1 - Satter_eta0,4)*100 , 
  #                                                "% based on Satterthwaite F-test p-values" ))

ggplot(KR_Ftest_pval_df[KR_Ftest_pval_df$treatment == "noGM",], aes(x = p_values, fill = method)) +
  geom_histogram(aes(y=..density..), bins = 70, position = "identity", alpha = .4) +
  labs(x = "F-test p-values", title = "Density histograms of F-test p-values for noGM patients") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

#---------------------------------------------------------------------------------------------
# try locfdr

GM_locfdr_KR <- locfdr(GM_result[,"GM_zval_1sided_KR"],
                       nulltype = 0, 
                        df = 9, # to fit estimated f(z) ,
                        # mlests = c(-1.5, 1.4),
                       # mlests = c(-1, 2), # no signal
                        main = "GM's locFDR based on LMER t-stat 1-sided pval with KR-adjusted df")

GM_locfdr_Satterthwaite <- locfdr(GM_result[,"GM_zval_1sided_Satterthwaite"], 
                       nulltype = 0, 
                       df = 11,  # to fit estimated f(z) 
                       # mlests = c(-1.5, 1.4), # initial value for mean & sd of estimating f(0)
                       # mlests = c(-.5, 2), # no signal
                       main = "GM's locFDR based on LMER t-stat 1-sided pval with Satterthwaite-adjusted df")

# let's get what we want
locFDR_GM_KR <- GM_locfdr_KR$fdr
locFDR_GM_Satterthwaite <- GM_locfdr_Satterthwaite$fdr
count.func(locFDR_GM_KR, seq(.01,.1,by=.01))
count.func(locFDR_GM_Satterthwaite, seq(.01,.1,by=.01))

# just to compare
noGM_locfdr_KR <- locfdr(noGM_result[,"noGM_zval_1sided_KR"], 
                         nulltype = 0, 
                         df = 9, # to fit estimated f(z) ,
                       # mlests = c(-1.5, 1.4),
                       # mlests = c(-1, 2), # no signal
                       main = "noGM's locFDR based on LMER t-stat 1-sided pval with KR-adjusted df")

noGM_locfdr_Satterthwaite <- locfdr(noGM_result[,"noGM_zval_1sided_Satterthwaite"], 
                                    nulltype = 0, 
                                    df = 11,  # to fit estimated f(z) 
                                  # mlests = c(-1.5, 1.4), # initial value for mean & sd of estimating f(0)
                                  # mlests = c(-.5, 2), # no signal
                                  main = "noGM's locFDR based on LMER t-stat 1-sided pval with Satterthwaite-adjusted df")

#---------------------------------------------------------------------------------------------
# check volcano plot

plot_eff <- GM_result[,"GM_time_effect"]
plot_pval <- GM_result[,"GM_KR_Ftest_pval"]

# check
raw_cutoff <- which(plot_pval <= .05 & plot_eff >= .333)
length(raw_cutoff) # check

# compare with JITC
JITC_list <- read_csv("JITC_PAP_Longitudinal.csv")
sum(as.numeric( dat_log2$PROBE_ID[raw_cutoff] %in% JITC_list$PROBE_ID ))
plot_JITC_ind <- (dat_log2$PROBE_ID %in% JITC_list$PROBE_ID)

plot(x = plot_eff, y = -log10(plot_pval), 
     pch = ".", las = 1,
     # ylim = c(0,9), xlim = c(-.4, .9),
     xlab = "coefficient of time fixed effect", 
     ylab = "-log10(KR F-test p-values)",
     main = "Volcano Plot for Vaccine + GM-CSF group")
lines(x = plot_eff[ plot_pval <= .05 & plot_eff >= .333 ], 
      y = -log10(plot_pval[ plot_pval <= .05 & plot_eff >= .333 ]),
      type = "p", pch = ".", col = "red")
abline(v = 0.333, col = "blue", lty = 2)
# abline(v = -0.333, col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)

# extra plot no GM
plot_eff2 <- noGM_result[,"noGM_time_effect"]
plot_pval2 <- noGM_result[,"noGM_KR_Ftest_pval"]
plot(x = plot_eff2, y = -log10(plot_pval2), 
     pch = ".", las = 1,
     xlim = c( min(plot_eff), max(plot_eff) ),
     ylim = c( min(-log10(plot_pval)), max(-log10(plot_pval)) ),
     xlab = "coefficient of time fixed effect", 
     ylab = "-log10(KR F-test p-values)",
     main = "Volcano Plot for Vaccine-only group (no GM)")
lines(x = plot_eff2[ plot_pval2 <= .05 & plot_eff2 >= .333 ], 
      y = -log10(plot_pval2[ plot_pval2 <= .05 & plot_eff2 >= .333 ]),
      type = "p", pch = ".", col = "red")
abline(v = 0.333, col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)
plot_eff2[(plot_eff2 >= .333) & (plot_pval2 <= .05)]

# focusing on JITC 5680 peptides (with GM)
plot(x = plot_eff, y = -log10(plot_pval), 
     pch = ".", las = 1,
     # ylim = c(0,9), xlim = c(-.4, .9),
     xlab = "coefficient of time fixed effect", 
     ylab = "-log10(KR F-test p-values)",
     main = "Volcano Plot for Vaccine + GM-CSF group")
lines(x = plot_eff[ plot_JITC_ind ], 
      y = -log10(plot_pval[ plot_JITC_ind ]),
      type = "p", pch = ".", col = "red")
abline(v = 0.333, col = "blue", lty = 2)
# abline(v = -0.333, col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)

# focusing on JITC 5680 peptides (without GM)
plot_eff2 <- noGM_result[,"noGM_time_effect"]
plot_pval2 <- noGM_result[,"noGM_KR_Ftest_pval"]
plot(x = plot_eff2, y = -log10(plot_pval2), 
     pch = ".", las = 1,
     xlim = c( min(plot_eff), max(plot_eff) ),
     ylim = c( min(-log10(plot_pval)), max(-log10(plot_pval)) ),
     xlab = "coefficient of time fixed effect", 
     ylab = "-log10(KR F-test p-values)",
     main = "Volcano Plot for Vaccine-only group (no GM)")
lines(x = plot_eff2[ plot_JITC_ind ], 
      y = -log10(plot_pval2[ plot_JITC_ind ]),
      type = "p", pch = ".", col = "red")
abline(v = 0.333, col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)
