##SECONDARY ANALYSIS AFTER IMPUTED LABORATORY DATASETS OBTAINED (WRITTEN BY JIANG LI)
#CLICK "Code" then "Jump to" each section. 

# required libraries ------------------------------------------------------
library(VIM)
library(FactoMineR)
library(missMDA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(mice) #R mice package for impuation
library(factoextra) #R factoextra package for PCA
library(splitstackshape)
library(miceadds)
library(randomForest)
library(Rmisc)
library(lattice)
library(GGally)
library(visdat)


# missing pattern ---------------------------------------------------------

head(labData_gs_transformed_clean)

PC <- read.csv("path_to_PC_file.csv", header = T, stringsAsFactors = F)  #PCs derived from PCA analysis

labData_gs_transformed_clean_pc <- merge(labData_gs_transformed_clean, PC, by = "ID")
head(labData_gs_transformed_clean_pc)
tail(labData_gs_transformed_clean_pc)

# Scatter plot matrix using the base R plot function
tiff("test_scatterplot.tiff", units="in", width=8, height=8, res=300)
pairs(labData_gs_transformed_clean_pc[labData_gs_transformed_clean_pc$TIME==2, 2:10])
dev.off()

labData_gs_transformed_clean_pc_select <- labData_gs_transformed_clean_pc[, colnames(labData_gs_transformed_clean_pc) %in% c("ID", "TIME", "X13457.7", "X17856.6", "X2160.0", "X2823.3", "X788.0", "X770.8")]

# Scatter plot matrix with lattice  
tiff("test_scatterplot_splom.tiff", units="in", width=8, height=8, res=300)
splom( ~ labData_gs_transformed_clean_pc_select[labData_gs_transformed_clean_pc_select$TIME==2, 2:7])
dev.off()

# Scatter plot matrix colored by groups
tiff("test_scatterplot_splom.tiff", units="in", width=8, height=8, res=300)
splom( ~ labData_gs_transformed_clean_pc_select[,2:7],  pch = 16, col = as.factor(labData_gs_transformed_clean_pc_select$TIME))
dev.off()


# Produce a matrix of plots for the first four variables 
labData.gg <- ggpairs(data = labData_gs_transformed_clean_pc_select, columns = 2:7)
labData.gg


# Color the points by TIME (before the event or after the event) 
tiff("test_scatterplot_splom1.tiff", units="in", width=10, height=8, res=300)
p <- ggpairs(data = labData_gs_transformed_clean_pc_select, mapping = aes(color = TIME, alpha = 0.7), columns = c("X13457.7", "X17856.6", "X2160.0", "X2823.3", "X788.0", "X770.8"), legend = 2) 
p + theme(legend.position = "top")

dev.off()

head(labData_gs_transformed_clean) # original file without holdouts

#GNSIS before
labData_gs_transformed_clean_before <- labData_gs_transformed_clean[labData_gs_transformed_clean$TIME == 2, ] #2 represents "before"

PC <- read.csv("ptIcd9_PC20_20percent_3years.csv", header = T, stringsAsFactors = F)

head(PC)
labData_gs_transformed_clean_before_PC <- merge(labData_gs_transformed_clean_before, PC, by = "ID")
head(labData_gs_transformed_clean_before_PC)

library(naniar)
gg_miss_var(labData_gs_transformed_clean_before_PC[, 1:47])

#GNSIS after
labData_gs_transformed_clean_after <- labData_gs_transformed_clean[labData_gs_transformed_clean$TIME == 1, ] #1 represents "after"

labData_gs_transformed_clean_after_PC <- merge(labData_gs_transformed_clean_after, PC, by = "ID")

gg_miss_var(labData_gs_transformed_clean_after_PC[, 1:47])

all_missing = as.data.frame(apply(labData_gs_transformed_clean_after_PC[, 1:47], 2, function(x){sum(is.na(x))/length(x)}))
all_missing$colnames <- rownames(all_missing)
colnames(all_missing) <- c("percent_missingness", "labs")
all_missing <- all_missing[order(all_missing$percent_missingness),]
labData_gs_transformed_clean_before_PC <- labData_gs_transformed_clean_before_PC[, colnames(labData_gs_transformed_clean_before_PC) %in% all_missing$labs]
labData_gs_transformed_clean_before_PC <- labData_gs_transformed_clean_before_PC %>%
  dplyr::select(all_missing$labs)

labData_gs_transformed_clean_before_PC <- labData_gs_transformed_clean_before_PC[orderCompletely(labData_gs_transformed_clean_before_PC), ]
head(labData_gs_transformed_clean_before_PC)

vis_dat(labData_gs_transformed_clean_before_PC[, 1:47])

labData_gs_transformed_clean_after_PC <- labData_gs_transformed_clean_after_PC[, colnames(labData_gs_transformed_clean_after_PC) %in% all_missing$labs]
labData_gs_transformed_clean_after_PC <- labData_gs_transformed_clean_after_PC %>%
  dplyr::select(all_missing$labs)

labData_gs_transformed_clean_after_PC <- labData_gs_transformed_clean_after_PC[orderCompletely(labData_gs_transformed_clean_after_PC), ]
head(labData_gs_transformed_clean_after_PC)

vis_dat(labData_gs_transformed_clean_after_PC[, 1:47])

png("Rplot_GNSIS_labData_gs_transformed_clean_after_missing.png",width=15.25,height=4.25,units="in",res=1200)
par(mar=c(1.4,1.4,1.4,2.1))
vis_miss(labData_gs_transformed_clean_after_PC[, 1:47])
dev.off()

png("Rplot_GNSIS_labData_gs_transformed_clean_before_missing1.png",width=15.25,height=4.25,units="in",res=1200)
par(mar=c(1.4,1.4,1.4,2.1))
vis_miss(labData_gs_transformed_clean_before_PC[, 1:47])
dev.off()

#marginplot
marginplot(labData_gs_transformed_clean_before_PC[,c("X13457.7","X742.7")])
marginplot(labData_gs_transformed_clean_after_PC[,c("X13457.7","X742.7")])

marginplot(labData_gs_transformed_clean_before_PC[,c("X13457.7","X2823.3")])
marginplot(labData_gs_transformed_clean_after_PC[,c("X13457.7","X2823.3")])

marginplot(labData_gs_transformed_clean_before_PC[,c("X13457.7","X6301.6")])
marginplot(labData_gs_transformed_clean_after_PC[,c("X13457.7","X6301.6")])

marginplot(labData_gs_transformed_clean_before_PC[,c("X13457.7","X788.0")])
marginplot(labData_gs_transformed_clean_after_PC[,c("X13457.7","X788.0")])


# missing pattern for repeated measure ------------------------------------

str(df.data_before)
sum(is.na(df.data_before$LAB_VALUE_TXT)) 
sum(is.na(df.data_after$LAB_VALUE_TXT))  

df.data_before_valid <- df.data_before %>% 
  filter(LAB_VALUE_TXT != 'NA') %>% 
  distinct(PT_ID, LAB_LOINC_CD, diff, .keep_all = T) 
sum(is.na(df.data_before_valid$LAB_VALUE_TXT)) #0

df.data_before_valid_3year_count <- df.data_before_valid %>%
  filter(diff > -1095) %>%
  group_by(PT_ID, LAB_LOINC_CD) %>%
  summarise(count=n()) %>%
  ungroup()

df.data_before_valid_3year_count <- data.frame(df.data_before_valid_3year_count)
df.data_before_valid_3year_count$TIME <- 2

df.data_after_valid_3year_count <- df.data_after_valid %>%
  filter(diff < 1095) %>%
  group_by(PT_ID, LAB_LOINC_CD) %>%
  summarise(count=n()) %>%
  ungroup()

df.data_after_valid_3year_count <- data.frame(df.data_after_valid_3year_count)
df.data_after_valid_3year_count$TIME <- 1

df.data_before_after_valid_3year <- rbind(df.data_before_valid[df.data_before_valid$diff > -1095, ], df.data_after_valid[df.data_after_valid$diff < 1095, ])
df.data_before_after_valid_3year_count <- df.data_before_after_valid_3year %>%
  group_by(PT_ID, LAB_LOINC_CD) %>%
  summarise(count=n()) %>%
  ungroup()

df.data_before_after_valid_3year_count <- data.frame(df.data_before_after_valid_3year_count)
df.data_before_after_valid_3year_count$TIME <- 3

# visualize the distribution before and after imputation  -----------------

#density plots
plot(imp9, c("X511.x", "X548.x", "X601.x", "X814.x", "X10774.x"))

tiff("2lpan_monotone_imp9.tiff", units="in", width=15, height=10, res=600)
densityplot(imp9)
dev.off()

tiff("2lpan_monotone_timeslope_imp10.tiff", units="in", width=15, height=10, res=600)
densityplot(imp10)
dev.off()

tiff("2lpan_monotone_5pc_20percent_imp11.tiff", units="in", width=15, height=10, res=600)
densityplot(imp11)
dev.off()

sum(is.na(complete(imp9))) 
#[1] 0
densityplot(imp5, ~ X17856.6.x, ylim = c(0, 0.015), xlim = c(300, 800))
densityplot(imp5, ~ X511.x, ylim = c(0, 0.30), xlim = c(0, 50))
densityplot(imp5, ~ X624.x, ylim = c(0, 0.06), xlim = c(0, 100))

##First create a long format version of all your datasets:

d.long <- mice::complete(imp9,"long",include = T)
head(d.long)
unique(d.long$.imp)
##Next perform your grouping as normal using base R

d.long.A <- d.long[which(d.long$TIME == 2),]
d.long.B <- d.long[which(d.long$TIME == 1),]
#Then change these back to mids objects, so you can perform mice operations

imp.A <- as.mids(d.long.A)
imp.B <- as.mids(d.long.B)

tiff("densityplot_imp_A_GNSIS_2LPAN_MONOTONE.tiff", units="in", width=15, height=10, res=600)
densityplot(imp.A)
dev.off()

tiff("densityplot_imp_B_GNSIS_2LPAN_MONOTONE.tiff", units="in", width=15, height=10, res=600)
densityplot(imp.B)
dev.off()

d.long <- mice::complete(imp21,"long",include = T)
head(d.long)
unique(d.long$.imp)
d.long.A <- d.long[which(d.long$TIME == 2),]
d.long.B <- d.long[which(d.long$TIME == 1),]
imp.A <- as.mids(d.long.A)
imp.B <- as.mids(d.long.B)

tiff("densityplot_imp_A_HF_2LPAN_MONOTONE.tiff", units="in", width=15, height=10, res=600)
densityplot(imp.A)
dev.off()

tiff("densityplot_imp_B_HF_2LPAN_MONOTONE.tiff", units="in", width=15, height=10, res=600)
densityplot(imp.B)
dev.off()

tiff("densityplot_imp_HF_2LPAN_MONOTONE.tiff", units="in", width=15, height=10, res=600)
densityplot(imp21)
dev.off()


# flux plot ---------------------------------------------------------------

head(labData_gs_before_after_transformed_clean)
flux(labData_gs_before_after_transformed_clean, local = names(labData_gs_before_after_transformed_clean))
#remove some variables if any.
labData_gs_before_after_transformed_clean_remove <- labData_gs_before_after_transformed_clean[, !colnames(labData_gs_before_after_transformed_clean) %in% c("X49541.6", "2965.2")]
#before
q <- flux(labData_gs_before_after_transformed_clean_remove[labData_gs_before_after_transformed_clean_remove$XTIME ==2, 1:(ncol(labData_gs_before_after_transformed_clean_remove)-2)], local = names(labData_gs_before_after_transformed_clean_remove[labData_gs_before_after_transformed_clean_remove$XTIME ==2, 1:(ncol(labData_gs_before_after_transformed_clean_remove)-2)])
)
#after
q <- flux(labData_gs_before_after_transformed_clean_remove[labData_gs_before_after_transformed_clean_remove$XTIME ==1, 1:(ncol(labData_gs_before_after_transformed_clean_remove)-2)], local = names(labData_gs_before_after_transformed_clean_remove[labData_gs_before_after_transformed_clean_remove$XTIME ==2, 1:(ncol(labData_gs_before_after_transformed_clean_remove)-2)])
)

tiff("test12.tiff", units="in", width=8, height=6, res=600)
ggplot(q, aes(x=influx, y=outflux)) +
  geom_point() + # Show dots
  geom_label(
    label=rownames(q), 
    nudge_x = 0.02, nudge_y = 0.02, 
    check_overlap = T,
    fill = alpha(c("white"),0.8)
  ) +
  theme_bw()
dev.off()

# load original file, holdout file, and imputed file ----------------------

load("labData_gs_before_after_transformed_clean_holdout.RData")
head(labData_gs_before_after_transformed_clean_holdout) #18074*49
head(labData_gs_before_after_transformed_clean)
write.table(labData_gs_before_after_transformed_clean, file = "labData_gs_before_after_transformed_clean.txt", col.names = T, row.names = T, quote = F)
write.table(labData_gs_before_after_transformed_clean_holdout, file = "labData_gs_before_after_transformed_clean_holdout.txt", col.names = T, row.names = T, quote = F)

labData_gs_before_after_transformed_clean$ID <- rownames(labData_gs_before_after_transformed_clean)
transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
  select(ID, X13457.7, XTIME) %>%
  filter(XTIME == 2)

sum(is.na(transformed_clean_before$X13457.7))/9037  #0.5513998

head(labData_gs_before_after_transformed_clean_holdout)
dim(labData_gs_before_after_transformed_clean_holdout)
#[1] 18074    49
transformed_clean_holdout_before <- labData_gs_before_after_transformed_clean_holdout %>%
  select(ID, X13457.7, TIME) %>%
  filter(TIME == 2)
sum(is.na(transformed_clean_holdout_before$X13457.7))/9037  #0.5569326

(0.5569326-0.5513998)*9037 = 50
merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")

#short format
imputedData
#long format
imputedDataRepeat 

# functions for error metric ----------------------------------------------

calculateMAPE = function(estimate, parameter, method = "raw") {
  all_vals = c(estimate, parameter)
  df = data.frame(estimate, parameter)
  # 0.01 is added to make sure nonzero value
  mape = sum((abs(df$parameter - df$estimate)) / (abs(parameter) + 0.01)) / sum(!is.na(df$estimate)) * 100  
  return(mape)
}

#nRMSE using raw
calculateRMSE = function(estimate, parameter, method = "raw") {
  all_vals = c(estimate, parameter)
  df = data.frame(estimate, parameter)
  rmse = sqrt( sum((df$estimate - df$parameter)^2, na.rm = TRUE) / sum(!is.na(df$estimate)))
  if (method == "minmax") {
    rmse = rmse / (max(all_vals, na.rm = T) - min(all_vals, na.rm = T))
  } else if (method == "sd") {
    rmse = rmse / sd(all_vals, na.rm = T)
  } else if (method == "raw") {
    rmse = rmse
  }
  return(rmse)
}

#nRMSE using sd
calculateRMSE = function(estimate, parameter, method = "sd") {
  all_vals = c(estimate, parameter)
  df = data.frame(estimate, parameter)
  rmse = sqrt( sum((df$estimate - df$parameter)^2, na.rm = TRUE) / sum(!is.na(df$estimate)))
  if (method == "minmax") {
    rmse = rmse / (max(all_vals, na.rm = T) - min(all_vals, na.rm = T))
  } else if (method == "sd") {
    rmse = rmse / sd(all_vals, na.rm = T)
  } else if (method == "raw") {
    rmse = rmse
  }
  return(rmse)
}

#calculate CR (95% or 99% CI) and AW

calculateCR <- function(ID, estimate, parameter){
  all_vals = c(ID, estimate, parameter)
  df=data.frame(ID, estimate, parameter)
  CI <- df %>%
    group_by(ID) %>%
    dplyr::summarise(avg = ci(estimate)[1], 
                     uci = ci(estimate)[3], 
                     lci = ci(estimate)[2])
  summary_final <- merge(df, CI, by = "ID") %>%
    mutate(coverage = ifelse(estimate <= uci & estimate >= lci, 0, 1)) %>%
    mutate(width = uci-lci)
  summary_final1 <- summary_final %>%
    group_by(ID) %>%
    dplyr::summarise(CR = mean(coverage), AW = mean(width))
  result[1] <- mean(summary_final1$CR)
  result[2] <- CI(summary_final1$CR)[1]
  result[3] <- CI(summary_final1$CR)[3]
  result[4] <- mean(summary_final1$AW)
  result[5] <- CI(summary_final1$AW)[1]
  result[6] <- CI(summary_final1$AW)[3]
  return(result)
}

calculateCR99 <- function(ID, estimate, parameter){
  all_vals = c(ID, estimate, parameter)
  df=data.frame(ID, estimate, parameter)
  CI <- df %>%
    group_by(ID) %>%
    dplyr::summarise(avg = mean(estimate), 
                     uci = mean(estimate) + 1.645*(sd(estimate)/sqrt(50)), 
                     lci = mean(estimate) - 1.645*(sd(estimate)/sqrt(50)))
  summary_final <- merge(df, CI, by = "ID") %>%
    mutate(coverage = ifelse(estimate <= uci & estimate >= lci, 0, 1)) %>%
    mutate(width = uci-lci)
  summary_final1 <- summary_final %>%
    group_by(ID) %>%
    dplyr::summarise(CR = mean(coverage), AW = mean(width))
  result[1] <- mean(summary_final1$CR)
  result[2] <- CI(summary_final1$CR)[1]
  result[3] <- CI(summary_final1$CR)[3]
  result[4] <- mean(summary_final1$AW)
  result[5] <- CI(summary_final1$AW)[1]
  result[6] <- CI(summary_final1$AW)[3]
  return(result)
}

#conduct CR & AW for imputed complete sets in long format
CR_AW <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()

Result_CR_AW = function(original_file, holdout_file, imputedData_pmm) {
  labs <- unique(imputedData_pmm$labs) #for GNSIS monotone or nonmonotone (long format)
  imputedData_pmm <- imputedData_pmm %>%
    filter(XTIME == 2) %>% #only select data for "before" index date. change to XTIME when working on holdout_case
    dplyr::select(ID, labs, score, repeatnum) %>%
    spread(labs, score) %>%
    dplyr::select(sort(names(.)))
  for (i in labs){
      print(i)
      transformed_clean_before <- original_file %>%
        dplyr::select(ID, i, XTIME) %>% 
        filter(XTIME == 2)   
      transformed_clean_holdout_before <- holdout_file %>%
        dplyr::select(ID, i, XTIME) %>%
        filter(XTIME == 2)
      merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
      merge <- cbind(merge_test,imputedData_pmm[, c(i, "ID")])
      colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", i, "Index")
      merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
      results[i] = data.frame(calculateCR(merge[,1], merge[,6], merge[,2]))  #calculate 99%CI
  }
  CR_AW <- as.data.frame(do.call(rbind,results)) #dat$lab <- i
  colnames(CR_AW) <- c("mean_coverage_rate", "UCI_coverage_rate", "LCI_coverage_rate", "mean_average_width", "UCI_average_width", "LCI_average_width")
  return(CR_AW)
}

#conduct CR & AW for imputed complete sets in short format
CR_AW <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()

Result_CR_AW = function(original_file, holdout_file, imputedData_pmm) {
  labs <- colnames(imputedData_pmm[, 1:(ncol(imputedData_pmm) -3)]) #for GNSIS monotone or nonmonotone (short format)
  for (i in labs){
    print(i)
    transformed_clean_before <- original_file %>%
      dplyr::select(ID, i, TIME) %>%   
      filter(TIME == 2)   
    transformed_clean_holdout_before <- holdout_file %>%
      dplyr::select(ID, i, TIME) %>%
      filter(TIME == 2)
    merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
    merge <- cbind(merge_test,imputedData_pmm[, c(i, "ID")])
    colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", i, "Index")
    merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
    results[i] = data.frame(calculateCR(merge[,1], merge[,6], merge[,2])) #calculate 99%CI change into calculatedCR99
  }
  CR_AW <- as.data.frame(do.call(rbind,results))
  colnames(CR_AW) <- c("mean_coverage_rate", "UCI_coverage_rate", "LCI_coverage_rate", "mean_average_width", "UCI_average_width", "LCI_average_width")
  return(CR_AW)
}

#holdout value
Result_CR_AW_pmm <- Result_CR_AW(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedData_pmm_monotone)
#holdout case
Result_CR_AW_pmm <- Result_CR_AW(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout_case, imputedData_pmm_monotone)

head(Result_CR_AW_pmm)
Result_CR_AW_pmm$lab <- rownames(Result_CR_AW_pmm)
Result_CR_AW_pmm$model <- "2LPAN_FCS"
Result_CR_AW_pmm_2lpan_monotone <- Result_CR_AW_pmm

Result_2lpan <- rbind(Result_CR_AW_pmm_2lpan_monotone, Result_CR_AW_pmm_2lpan_nonmonotone)
Result_pmm <- rbind(Result_CR_AW_pmm_pmm_monotone, Result_CR_AW_pmm_pmm_nonmonotone)
Result_ALL <- rbind(Result_pmm, Result_2lpan)
Result_UNIVARIATE <- rbind(Result_CR_AW_pmm_2lpan_univariate_0pc, Result_CR_AW_pmm_2lpan_univariate_5pc)
Result <- rbind(Result_ALL, Result_UNIVARIATE)
head(Result_ALL)

tiff("test12.tiff", units="in", width=15, height=6, res=600)
ggplot(Result_ALL, aes(x = lab, y = mean_coverage_rate, ymin = LCI_coverage_rate, ymax = UCI_coverage_rate)) + 
  geom_pointrange(aes(col = factor(model)), position=position_dodge(width=0.8)) + 
  geom_point(aes(x = lab, y = mean_coverage_rate, size = 0.5, col = factor(model)), position=position_dodge(width=0.8)) +
  ylab("Mean Coverage Rate") + geom_hline(aes(yintercept = 0.90)) + scale_color_discrete(name = "model") + xlab("Laboratory Variables") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

tiff("test13.tiff", units="in", width=15, height=6, res=600)
ggplot(Result, aes(x = lab, y = mean_average_width, ymin = LCI_average_width, ymax = UCI_average_width)) + 
  geom_pointrange(aes(col = factor(model)), position=position_dodge(width=0.8)) + 
  scale_y_log10() +
  geom_point(aes(x = lab, y = mean_average_width, size = 0.5, col = factor(model)), position=position_dodge(width=0.8)) +
  ylab("Mean Average Width") + scale_color_discrete(name = "model") + xlab("Laboratory Variables") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


# calculate nRMSE for imputed files in long and short format --------------

#the following code is to calculate nRMSE for 2l.pan longformat imputed files
RMSE <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()

Result_RMSE = function(original_file, holdout_file, imputedData_pmm) {
  labs <- unique(imputedData_pmm$labs)
  for (j in 1:50){
    print(j)
    imputedData <- imputedData_pmm %>%
      filter(repeatnum == j) %>%
      filter(XTIME == 2) %>% #only select data for "before" index date. change to XTIME when working on holdout_case
      dplyr::select(ID, labs, score) %>%  
      spread(labs, score) %>%
      dplyr::select(sort(names(.)))
    for (i in labs){
      print(i)
      transformed_clean_before <- original_file %>%
        dplyr::select(ID, i, XTIME) %>%  
        filter(XTIME == 2) 
      transformed_clean_holdout_before <- holdout_file %>%
        dplyr::select(ID, i, XTIME) %>%  
        filter(XTIME == 2)  
      merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
      merge <- cbind(merge_test,imputedData[, c(i, "ID")])
      colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", i, "Index")
      merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
      results[i] = calculateRMSE(merge[,2], merge[,6]) #should change 2 to 6 and 6 to 2
    }
    dat <- data.frame(results)
    dat$repeatnum <- j
    datalist[[j]] <- dat 
    print(datalist)
  }
  RMSE <- as.data.frame(do.call(rbind,datalist))
  return(RMSE)
}

#holdout value
Result_RMSE_pmm <- Result_RMSE(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedDataRepeat)

dim(Result_RMSE_pmm) #for 50 repeats 2300*2
Result_RMSE_pmm$model <- "2LPAN_UNIVARIATE"
Result_RMSE_pmm$PCA <- "5pcs"
Result_RMSE_pmm$CM <- "20_percent"
Result_RMSE_pmm$lab <- rownames(Result_RMSE_pmm)

#the following function is to calculate nRMSE for monotone imputed data

RMSE <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()

Result_RMSE = function(original_file, holdout_file, imputedData_pmm) {
  labs <- colnames(imputedData_pmm[, 1:(ncol(imputedData_pmm) -3)]) #for GNSIS monotone or nonmonotone
  for (j in 1:50){
    print(j)
    imputedData_pmm$ID <- gsub("PT", "", imputedData_pmm$ID)
    imputedData <- imputedData_pmm[imputedData_pmm$repeatnum == j, ]
    for (i in labs){
      print(i)
      transformed_clean_before <- original_file %>%
        dplyr::select(ID, i, TIME) %>%   #change TIME to XTIME for HF
        filter(TIME == 2)   #change TIME to XTIME for HF
      transformed_clean_holdout_before <- holdout_file %>%
        dplyr::select(ID, i, TIME) %>%
        filter(TIME == 2)
      merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
      merge <- cbind(merge_test,imputedData[, c(i, "ID")])
      #print(head(merge))
      colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", i, "Index")
      #merge <- cbind(merge_test,imputedData$i)
      merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
      #print(head(merge))
      results[i] = calculateRMSE(merge[,2], merge[,6])
    }
    dat <- data.frame(results)
    dat$repeatnum <- j
    datalist[[j]] <- dat
    print(datalist)
  }
  RMSE <- as.data.frame(do.call(rbind,datalist))
  return(RMSE)
}

Result_RMSE_pmm <- Result_RMSE(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedData_pmm_monotone)
dim(Result_RMSE_pmm)
Result_RMSE_pmm$model <- "2LPAN_FCS"
Result_RMSE_pmm$PCA <- "0pcs"
Result_RMSE_pmm$CM <- "0_percent"
Result_RMSE_pmm$lab <- rownames(Result_RMSE_pmm)
head(Result_RMSE_pmm)


# get summary statistics of nRMSE -----------------------------------------
#short format
RMSE <- data.frame()
results <- NULL
labs <- list()
imputed_result <- data.frame()

Result_RMSE_pmm_addlabs <- function(Result_RMSE_pmm, imputedData_pmm_monotone){
  labs <- colnames(imputedData_pmm_monotone[, 3:(ncol(imputedData_pmm_monotone) -1)])
  for (j in 1:length(unique(Result_RMSE_pmm$repeatnum))){
    print(j)
    imputed_result <- Result_RMSE_pmm[Result_RMSE_pmm$repeatnum == j,]
    imputed_result$lab <- labs
    results[[j]] <- imputed_result
  }
  RMSE <- as.data.frame(do.call(rbind,results))
  return(RMSE)
}

Result_RMSE_split <- Result_RMSE_pmm_addlabs(Result_RMSE_pmm, imputedData_pmm_monotone)
head(Result_RMSE_split)
data <- Result_RMSE_split %>%
  group_by(lab) %>%
  summarise(rmse_median = median(results), rmse_mean = mean(results), rmse_sd = sd(results)) #tibble
data <- data.frame(data)


# long format
RMSE <- data.frame()
results <- NULL
labs <- list()
imputed_result <- data.frame()

Result_RMSE_pmm_addlabs <- function(Result_RMSE_pmm, imputedDataRepeat){
  labs <- unique(imputedDataRepeat$labs)
  for (j in 1:length(unique(Result_RMSE_pmm$repeatnum))){
    print(j)
    imputed_result <- Result_RMSE_pmm[Result_RMSE_pmm$repeatnum == j,]
    imputed_result$lab <- labs
    results[[j]] <- imputed_result
  }
  RMSE <- as.data.frame(do.call(rbind,results))
  return(RMSE)
}


Result_RMSE_split <- Result_RMSE_pmm_addlabs(Result_RMSE_pmm, imputedDataRepeat)
head(Result_RMSE_split)
data <- Result_RMSE_split %>%
  group_by(lab) %>%
  summarise(rmse_median = median(results), rmse_mean = mean(results), rmse_sd = sd(results)) #tibble
data <- data.frame(data)


# t-test of nRMSE from two 'Result_RMSE' files ----------------------------
#write a function to conduct summary of t-test of two Result_RMSE files
#PCA
df <- NULL
summary <- NULL
df_final <- NULL
labs <- list()
summary_final <- NULL
dflist <- list()
summarylist <- list()
#n <- length(unique(Result_RMSE_merge_split$lab_1))
Result_RMSE_ttest = function(Result_RMSE_merge_split) {
  labs <- unique(Result_RMSE_merge_split$lab_1)
  for (i in labs){
    print(i)
    df <- Result_RMSE_merge_split %>% 
      filter(lab_1 == i) %>% 
      do(tidy(t.test(results ~ PCA, data = ., 
                     alt = "two.sided", 
                     paired = FALSE, 
                     conf.level = 0.99)))
    df <- data.frame(df)
    df$lab <- i
    dflist[[i]] <- df
    print(dflist)
  }
  df_final <- as.data.frame(do.call(rbind,dflist))
  return(df_final)
}

#model
df <- NULL
summary <- NULL
df_final <- NULL
labs <- list()
summary_final <- NULL
dflist <- list()
summarylist <- list()
Result_RMSE_ttest = function(Result_RMSE_merge_split) {
  labs <- unique(Result_RMSE_merge_split$lab_1)
  for (i in labs){
    print(i)
    df <- Result_RMSE_merge_split %>% 
      filter(lab_1 == i) %>% 
      do(tidy(t.test(results ~ model, data = ., 
                     alt = "two.sided", 
                     paired = FALSE, 
                     conf.level = 0.99)))
    df <- data.frame(df)
    df$lab <- i
    dflist[[i]] <- df
    print(dflist)
  }
  df_final <- as.data.frame(do.call(rbind,dflist))
  return(df_final)
}

ttest_final <- Result_RMSE_ttest(Result_RMSE_merge_split)


# functions for obtaining imputed value for holdouts only -----------------
# this is for 2l.pan longformat result
Value <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()
unique(imputedDataRepeat$labs)
Result_imputedHoldOut = function(original_file, holdout_file, imputedData_pmm) {
  labs <- unique(imputedData_pmm$labs)
  for (j in 1:50){
    print(j)
    imputedData <- imputedData_pmm %>%
      filter(repeatnum == j) %>%
      filter(XTIME == 2) %>% #only select data for "before" index date. #TIME change to XTIME when using holdoutcase
      dplyr::select(ID, labs, score) %>%  
      spread(labs, score) %>%
      dplyr::select(sort(names(.)))
    for (i in labs){
      print(i)
      transformed_clean_before <- original_file %>%
        dplyr::select(ID, i, XTIME) %>%  #change from XTIME to TIME
        filter(XTIME == 2)  #change from XTIME to TIME
      transformed_clean_holdout_before <- holdout_file %>%
        dplyr::select(ID, i, XTIME) %>%   #change from TIME to XTIME when using holdoutcase
        filter(XTIME == 2)
      merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
      merge <- cbind(merge_test,imputedData[, c(i, "ID")])
      colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", "Score", "Index")
      merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
      merge$labs <- i
      results[[i]] = merge
    }
    dat <- as.data.frame(do.call(rbind,results))
    dat$repeatnum <- j
    datalist[[j]] <- dat
  }
  Value <- as.data.frame(do.call(rbind,datalist))
  return(Value)
}

Result_imputedHoldOut_labs <- Result_imputedHoldOut(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedDataRepeat)

#monotone data imputed result

Value <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()
head(imputedData_pmm_monotone)

Result_imputedHoldOut = function(original_file, holdout_file, imputedData_pmm) {
  imputedData_pmm$ID <- gsub("PT", "", imputedData_pmm$ID)
  labs <- colnames(imputedData_pmm[, 1:(ncol(imputedData_pmm) -3)]) 
  for (j in 1:50){
    print(j)
    imputedData <- imputedData_pmm[imputedData_pmm$repeatnum == j, ]
    for (i in labs){
      print(i)
      transformed_clean_before <- original_file %>%
        dplyr::select(ID, i, XTIME) %>%
        filter(XTIME == 2)
      transformed_clean_holdout_before <- holdout_file %>%
        dplyr::select(ID, i, TIME) %>%
        filter(TIME == 2)
      merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
      merge <- cbind(merge_test,imputedData[, c(i, "ID")])
      colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", "Score", "Index")
      merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
      merge$labs <- i
      results[[i]] = merge
    }
    dat <- as.data.frame(do.call(rbind,results))
    dat$repeatnum <- j
    datalist[[j]] <- dat
  }
  Value <- as.data.frame(do.call(rbind,datalist))
  return(Value)
}

Result_imputedHoldOut_labs <- Result_imputedHoldOut(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedData_pmm_monotone)
dim(Result_imputedHoldOut_labs)
head(Result_imputedHoldOut_labs)
Result_imputedHoldOut_labs$model <- "2LPAN_FCS"
Result_imputedHoldOut_labs$PCA <- "0pcs"
Result_imputedHoldOut_labs$CM <- "0_percent"


# visualization of imputed data for holdouts using scatter plot -----------

Result_imputedHoldOut_labs_summary <- Result_imputedHoldOut_labs %>%
  dplyr::select(-XTIME, -y, -TIME, -Index, -repeatnum, -model, -PCA, -CM) %>%
  dplyr::group_by(ID, labs) %>%
  dplyr::summarise(imputed_median = median(Score), imputed_mean = mean(Score), imputed_sd = sd(Score), original = mean(x)) %>%
  ungroup()
Result_imputedHoldOut_labs_summary$imputed_upper = Result_imputedHoldOut_labs_summary$imputed_mean + (Result_imputedHoldOut_labs_summary$imputed_sd)/sqrt(50)
Result_imputedHoldOut_labs_summary$imputed_lower = Result_imputedHoldOut_labs_summary$imputed_mean - (Result_imputedHoldOut_labs_summary$imputed_sd)/sqrt(50)
Result_imputedHoldOut_labs_summary$imputed_se = Result_imputedHoldOut_labs_summary$imputed_sd/sqrt(50)
Result_imputedHoldOut_labs_summary <- as.data.frame(Result_imputedHoldOut_labs_summary)
#using HgA1c as an example
Result_imputedHoldOut_labs_summary_test <- Result_imputedHoldOut_labs_summary[Result_imputedHoldOut_labs_summary$labs == "X17856.6", ]

library(ggplot2)
library(ggExtra)
stat = cor.test(Result_imputedHoldOut_labs_summary_test$original, Result_imputedHoldOut_labs_summary_test$imputed_mean,  method = "pearson", use = "complete.obs")

tiff("test.tiff", units="in", width=6, height=6, res=600)
p <- ggplot(Result_imputedHoldOut_labs_summary_test, aes(x = original, y = imputed_mean, ymin = imputed_lower, ymax = imputed_upper)) + 
  geom_pointrange(position=position_dodge(width=0.8),  width=1, size=1, color="darkred") + 
  geom_point(aes(x = original, y = imputed_mean, size = 3), col = "red") +
  geom_smooth(method='lm', formula= y~x, cor.coef = TRUE, cor.method = "pearson") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  ggtitle("X17856.6") + annotate("text", x=max(Result_imputedHoldOut_labs_summary_test$original/1.5), y=max(Result_imputedHoldOut_labs_summary_test$imputed_mean), label= "paste(italic(R) ^ 2, \" = .759 (P = 1.65e-10)\")", parse = TRUE)
ggExtra::ggMarginal(p, type = "histogram") 
dev.off()

tiff("test0.tiff", units="in", width=6, height=6, res=600)
p <- ggplot(Result_imputedHoldOut_labs_summary_test, aes(x = original, y = imputed_mean, ymin = imputed_lower, ymax = imputed_upper)) + 
  geom_pointrange(position=position_dodge(width=0.8),  width=1, size=1, color="darkred") + 
  geom_point(aes(x = original, y = imputed_mean, size = 3), col = "red") +
  geom_smooth(method='lm', formula= y~x, cor.coef = TRUE, cor.method = "pearson") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  ggtitle("X17856.6") + annotate("text", x=max(Result_imputedHoldOut_labs_summary_test$original/1.5), y=max(Result_imputedHoldOut_labs_summary_test$imputed_mean), label= "paste(italic(R) ^ 2, \" = .759 (P = 1.65e-10)\")", parse = TRUE)
ggExtra::ggMarginal(p, type = "density") 
dev.off()


tiff("test01.tiff", units="in", width=20, height=20, res=600)
ggplot(Result_imputedHoldOut_labs_summary, aes(x = original, y = imputed_mean)) + 
  geom_point(aes(x = original, y = imputed_mean, size = imputed_se), col = "red") +
  geom_smooth(method='lm', formula= y~x) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(. ~ labs, scales="free") +
  geom_density(colour="red", aes(x = original, y=..density..)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~exp(.), breaks = NULL, name = "density")) + 
  geom_density(colour="blue", aes(x = imputed_mean, y=..density..))
dev.off()

tiff("test3.tiff", units="in", width=18, height=10, res=600)
ggplot(RMSE_plot) + 
  geom_point(aes(x = lab, y = RMSE, color = method, size = -log10(p.value)/10), fill = "white") +
  scale_size(name   = "p.value", 
             breaks = c(0.2, 0.5, 2, 5), 
             labels = c("0.01", "1e-5", "1e-20", "1e-50")) +
  geom_line(aes(x = lab, y = ICC1), size = 1.5, color="red", group = 1, linetype = "dashed")  +
  scale_y_continuous(sec.axis = sec_axis(~., name = "ICC1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))
dev.off()

# visualization of t-test result across different imputation algorithm  --------

#t-test result1
RMSE_plot1 
#t-test result2
RMSE_plot2 

head(RMSE_plot1)
RMSE_plot1A <- RMSE_plot1[, c("V1", "V3", "V6")]
colnames(RMSE_plot1A) <- c("lab", "RMSE", "p.value")
RMSE_plot1A$method <- "PMM_FCS"
RMSE_plot1A$holdout <- "CASE"
RMSE_plot1B <- RMSE_plot1[, c("V1", "V4", "V6")]
colnames(RMSE_plot1B) <- c("lab", "RMSE", "p.value")
RMSE_plot1B$method <- "PMM_MONOTONE"
RMSE_plot1B$holdout <- "CASE"
RMSE_plot1AB <- rbind(RMSE_plot1A, RMSE_plot1B)

head(RMSE_plot2)
RMSE_plot2A <- RMSE_plot2[, c("V1", "V3", "V6")]
colnames(RMSE_plot2A) <- c("lab", "RMSE", "p.value")
RMSE_plot2A$method <- "PMM_FCS"
RMSE_plot2A$holdout <- "VALUE"
RMSE_plot2B <- RMSE_plot2[, c("V1", "V4", "V6")]
colnames(RMSE_plot2B) <- c("lab", "RMSE", "p.value")
RMSE_plot2B$method <- "PMM_MONOTONE"
RMSE_plot2B$holdout <- "VALUE"
RMSE_plot2AB <- rbind(RMSE_plot2A, RMSE_plot2B)
RMSE_plot <- rbind(RMSE_plot1AB, RMSE_plot2AB)

ICC1 <- read.table("ICC1_summary_GNSIS.txt", header = T , stringsAsFactors = F)
ICC1 <- cSplit(ICC1, "labs", ".")
ICC1$lab <- ICC1$labs_1
ICC1$ICC1 <- ICC1$results

RMSE_plot_final <- merge(RMSE_plot, ICC1[, c("ICC1", "lab")], by = "lab")


tiff("test3.tiff", units="in", width=15, height=8, res=600)
ggplot(RMSE_plot_final) + 
  geom_point(aes(x = lab, y = RMSE, color = method, size = -log10(p.value)/10), fill = "white") +
  scale_size(name   = "p.value", 
             breaks = c(0.2, 0.5, 2, 5), 
             labels = c("0.01", "1e-5", "1e-20", "1e-50")) +
  geom_line(aes(x = lab, y = ICC1), size = 1.5, color="red", group = 1, linetype = "dashed")  +
  scale_y_continuous(sec.axis = sec_axis(~., name = "ICC1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_rect(fill = "lightgreen"),
        panel.grid.major = element_line(size = 0.8, linetype = 'solid', colour = "grey"),
        panel.border = element_blank()) + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
  facet_wrap(. ~ holdout, scales="free", ncol = 1) + 
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"))
dev.off()


# density plot for HgA1c before and after imputation for holdouts --------

unique(labData_gs_before_after_transformed_clean$ID)

transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
  dplyr::select(ID, X17856.6, XTIME) %>%
  filter(XTIME == 2)
transformed_clean_holdout_before <- labData_gs_before_after_transformed_clean_holdout %>%
  dplyr::select(ID, X17856.6, TIME) %>%
  filter(TIME == 2)
merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
merge_test$ID <- gsub("PT", "", merge_test$ID) 
merge <- merge(merge_test,imputedData_pmm_monotone[, c("X17856.6", "ID", "repeatnum")], by = "ID")
colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", "Score", "repeatnum")
merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
head(merge)
holdout <- unique(merge$ID)
merge_nonholdout <- imputedData_pmm_monotone[!imputedData_pmm_monotone$ID %in% holdout, c("X17856.6", "ID", "repeatnum")]
transformed_clean_before$ID <- gsub("PT", "", transformed_clean_before$ID)
merge_nonholdout <- merge(transformed_clean_before, merge_nonholdout, by = "ID")
colnames(merge_nonholdout) <- c("ID", "x", "TIME", "y", "repeatnum")
head(merge_nonholdout)

transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
  dplyr::select(ID, X17856.6, XTIME) %>%
  filter(XTIME == 2)
transformed_clean_holdout_case_before <- labData_gs_before_after_transformed_clean_holdout_case %>%
  dplyr::select(ID, X17856.6, XTIME) %>%
  filter(XTIME == 2)
merge_case_test <- merge(transformed_clean_before, transformed_clean_holdout_case_before, by = "ID")
merge_case_test$ID <- gsub("PT", "", merge_case_test$ID) 
merge_case <- merge(merge_case_test,imputedData_pmm_monotone_case[, c("X17856.6", "ID", "repeatnum")], by = "ID")
colnames(merge_case) <- c("ID", "x", "XTIME", "y", "TIME", "Score", "repeatnum")
merge_case <- merge_case[!is.na(merge_case$x) & is.na(merge_case$y), ]
head(merge_case)
holdout_case <- unique(merge_case$ID)
merge_case_nonholdout <- imputedData_pmm_monotone_case[!imputedData_pmm_monotone_case$ID %in% holdout_case, c("X17856.6", "ID", "repeatnum")]
transformed_clean_before$ID <- gsub("PT", "", transformed_clean_before$ID)
merge_case_nonholdout <- merge(transformed_clean_before, merge_case_nonholdout, by = "ID")
colnames(merge_case_nonholdout) <- c("ID", "x", "TIME", "y", "repeatnum")
head(merge_case_nonholdout)
length(is.na(merge_case_nonholdout$x))/nrow(merge_case_nonholdout) #8987
length(is.na(merge_nonholdout$x))/50 #8987
1- (9037-6404)/9037
2633
library(plyr)
mu <- ddply(merge, "repeatnum", summarise, grp.mean=mean(Score))
mu_non <- ddply(merge_nonholdout, "repeatnum", summarise, grp.mean=mean(y))
mu_case <- ddply(merge_case, "repeatnum", summarise, grp.mean=mean(Score)) 
mu_case_non <- ddply(merge_case_nonholdout, "repeatnum", summarise, grp.mean=mean(y)) 

head(mu)
#tiff("test10.tiff", units="in", width=5, height=2, res=300)
p1 <- ggplot(merge, aes(x=Score)) +
  geom_density(aes(group=as.factor(repeatnum)), linetype="dashed", alpha = 0.8, size = 0.1)+
  geom_vline(data=mu, aes(xintercept=grp.mean),
            linetype="dashed", size = 0.1) + 
  geom_density(aes(x=x), color = "red", alpha=.2, fill="#FF6666") +
  geom_vline(xintercept=mean(merge$x), linetype="dashed", 
              color = "red", size=1) + 
  ylab("density") + xlab("holdout value (n = 50; random value)") + xlim(0, 15)
#p1
#dev.off()

#tiff("test11.tiff", units="in", width=5, height=2, res=300)
p2 <- ggplot(merge_nonholdout, aes(x=y)) +
  geom_density(aes(group=as.factor(repeatnum)), linetype="dashed", alpha = 0.8, size = 0.1)+
  geom_vline(data=mu_non, aes(xintercept=grp.mean),
             linetype="dashed", size = 0.1) + 
  geom_density(aes(x=x), color = "red", alpha=.2, fill="#FF6666") +
  geom_vline(xintercept=mean(merge_nonholdout$x), linetype="dashed", 
             color = "red", size=5) + 
  ylab("density") + xlab("nonholdout value (n = 2583; random value)") + xlim(0, 15)
#p2
#dev.off()

#tiff("test12.tiff", units="in", width=5, height=2, res=300)
p3 <- ggplot(merge_case, aes(x=Score)) +
  geom_density(aes(group=as.factor(repeatnum)), linetype="dashed", alpha = 0.8, size = 0.1)+
  geom_vline(data=mu_case, aes(xintercept=grp.mean),
             linetype="dashed", size = 0.1) + 
  geom_density(aes(x=x), color = "red", alpha=.2, fill="#FF6666") +
  geom_vline(xintercept=mean(merge_case$x), linetype="dashed", 
             color = "red", size=1) + 
  ylab("density") + xlab("holdout value (n = 50; random case)") + xlim(0, 15)
#p3
#dev.off()

#tiff("test13.tiff", units="in", width=5, height=2, res=300)
p4 <- ggplot(merge_case_nonholdout, aes(x=y)) +
  geom_density(aes(group=as.factor(repeatnum)), linetype="dashed", alpha = 0.8, size = 0.1)+
  geom_vline(data=mu_case_non, aes(xintercept=grp.mean),
             linetype="dashed", size = 0.1) + 
  geom_density(aes(x=x), color = "red", alpha=.2, fill="#FF6666") +
  geom_vline(xintercept=mean(merge_case_nonholdout$x), linetype="dashed", 
             color = "red", size=2) +
  ylab("density") + xlab("nonholdout value (n = 2583; random case)") + xlim(0, 15)
#p4

tiff("test13.tiff", units="in", width=10, height=10, res=300)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()
#dev.off()

# Multiple plot function (acquired from other)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(plotlist=l, cols = 8)
multiplot(p1, p2, p3, p4, p5, p6, cols = 2)



# function for uncertainty evaluation based on nRMSE from # of repeated imputations --------

repeats <- NULL
RMSE <- data.frame()
merge <- NULL
results <- NULL
datalist <- list()
dat <- data.frame()
merge_test <- data.frame()
imputedData <- data.frame()
RMSE_final <- data.frame()
RMSElist <- NULL

Result_RMSE_NUM = function(original_file, holdout_file, imputedData_pmm){
  labs <- colnames(imputedData_pmm[, 1:(ncol(imputedData_pmm) -3)]) #for GNSIS monotone or nonmonotone
  repeats <- c(5, 10, 20, 30, 40, 50)
  for (n in repeats){
    print(n)
    for (j in 1:n){
      print(j)
      imputedData_pmm$ID <- gsub("PT", "", imputedData_pmm$ID)
      imputedData <- imputedData_pmm[imputedData_pmm$repeatnum == j, ]
      for (i in labs){
        print(i)
        transformed_clean_before <- original_file %>%
          dplyr::select(ID, i, XTIME) %>%   
          filter(XTIME == 2)   
        transformed_clean_holdout_before <- holdout_file %>%
          dplyr::select(ID, i, TIME) %>%  #change TIME to XTIME for holdoutcase
          filter(TIME == 2) #change TIME to XTIME for holdoutcase
        merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
        merge <- cbind(merge_test,imputedData[, c(i, "ID")])
        colnames(merge) <- c("ID", "x", "XTIME", "y", "TIME", i, "Index")
        merge <- merge[!is.na(merge$x) & is.na(merge$y), ]
        results[i] = calculateRMSE(merge[,2], merge[,6]) #calculate RMSE
      }
      dat <- data.frame(results)
      dat$repeatnum <- j
      datalist[[j]] <- dat
      print(datalist)
    }
    RMSE <- as.data.frame(do.call(rbind,datalist))
    RMSE$lab <- rownames(RMSE)
    Result_RMSE_split <- cSplit(RMSE, "lab", ".")
    data <- Result_RMSE_split %>%
      dplyr::select(-lab_2) %>% #add dplyr::
      group_by(lab_1) %>%
      dplyr::summarise(rmse_median = median(results), rmse_mean = mean(results), rmse_sd = sd(results))
    data <- data.frame(data)
    data$numrepeat <- n
    RMSElist[[n]] <- data
    print(RMSElist)
  }
  RMSE_final <- as.data.frame(do.call(rbind, RMSElist))
  return(RMSE_final)
}

Result_RMSE_NUM_pmm <- Result_RMSE_NUM(labData_gs_before_after_transformed_clean, labData_gs_before_after_transformed_clean_holdout, imputedData_pmm_monotone)
Result_RMSE_NUM_pmm$holdout <- "value"
Result_RMSE_NUM_pmm$method <- "2LPAN_FCS"

tiff("test5.tiff", units="in", width=20, height=20, res=600)
ggplot(Result_RMSE_NUM_final_monotone_FCS_pmm, aes(x = numrepeat, y = rmse_mean)) + 
  geom_point(aes(color = as.factor(holdout), size = rmse_sd, shape = as.factor(method)), stat="identity", position = position_dodge(0.8))  +
  scale_color_manual(values=c("#00AFBB", "#E7B800")) +
  scale_size(name   = "rmse_sd", 
             breaks = c(0.05, 0.10, 0.15, 0.20, 0.25), 
             labels = c("0.05", "0.10", "0.15", "0.20", "0.25")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(. ~ lab_1, scales="free")
dev.off()

tiff("test5.tiff", units="in", width=20, height=20, res=600)
ggplot(Result_RMSE_NUM_final_monotone_FCS_2lpan, aes(x = numrepeat, y = rmse_mean)) + 
  geom_point(aes(color = as.factor(holdout), size = rmse_sd, shape = as.factor(method)), stat="identity", position = position_dodge(0.8))  +
  scale_color_manual(values=c("#00AFBB", "#E7B800")) +
  scale_size(name   = "rmse_sd", 
             breaks = c(0.05, 0.10, 0.15, 0.20, 0.25), 
             labels = c("0.05", "0.10", "0.15", "0.20", "0.25")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(. ~ lab_1, scales="free")
dev.off()


# missing pattern (NMAR) related to PCs derived from cormorbidity matrix --------

dir.create("path_to_directory")
result <- NULL
result1 <- NULL
result2 <- NULL
transformed_clean_before <- NULL  #this is for original file without holdouts
transformed_clean_holdout_before <- NULL  #this is for holdout value
transformed_clean_holdoutcase_before <- NULL # this is for holdout case
merge <- NULL
merge_test <- NULL
holdout <- NULL
transformed_clean_before_merge <- NULL
transformed_clean_before_merge_filter <- NULL
result_final <- NULL

all_labs <- colnames(labData_gs_before_after_transformed_clean)[1:ncol(labData_gs_before_after_transformed_clean)] #1 could be replaced with number of the first column of lab variable


for (lab in all_labs) {
  print(lab)
  
  #holdout values
  transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
    dplyr::select(ID, lab, XTIME) %>%
    filter(XTIME == 2)
  transformed_clean_holdout_before <- labData_gs_before_after_transformed_clean_holdout %>%
    dplyr::select(ID, lab, TIME) %>%
    filter(TIME == 2)
  merge_test <- merge(transformed_clean_before, transformed_clean_holdout_before, by = "ID")
  colnames(merge_test) <- c("ID", "x", "XTIME", "y", "TIME")
  merge <- merge_test[!is.na(merge_test$x) & is.na(merge_test$y), ]
  holdout <- merge$ID
  transformed_clean_before <- transformed_clean_before %>%
    dplyr::mutate(MISSING = as.numeric(transformed_clean_before$ID %in% holdout))
  transformed_clean_before_merge <- merge(transformed_clean_before, PC, by = "ID")
  transformed_clean_before_merge_filter <- transformed_clean_before_merge %>%
    filter(!is.na(lab)) %>%
    tidyr::gather(key=DIM, value=score, 5:24) #select top 20 PCs located at 5 to 24 columns
  result <- transformed_clean_before_merge_filter %>%
    group_by(as.factor(DIM)) %>%
    do(tidy(t.test(score ~ MISSING, data = ., 
                   alt = "two.sided", 
                   paired = FALSE, 
                   conf.level = 0.99)))
  result <- result[order(result$`as.factor(DIM)`), ]
  #result$pattern <- result$estimate/abs(result$estimate) #option
  result$Category <- "holdout value"
  
  # holdout cases
  transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
    dplyr::select(ID, lab, XTIME) %>%
    filter(XTIME == 2)
  transformed_clean_holdoutcase_before <- labData_gs_before_after_transformed_clean_holdout_case %>%
    dplyr::select(ID, lab, XTIME) %>%
    filter(XTIME == 2)
  merge_test <- merge(transformed_clean_before, transformed_clean_holdoutcase_before, by = "ID")
  colnames(merge_test) <- c("ID", "x", "XTIME", "y", "TIME")
  merge <- merge_test[!is.na(merge_test$x) & is.na(merge_test$y), ]
  holdout <- merge$ID
  transformed_clean_before <- transformed_clean_before %>%
    mutate(MISSING = as.numeric(transformed_clean_before$ID %in% holdout))
  transformed_clean_before_merge <- merge(transformed_clean_before, PC, by = "ID")
  transformed_clean_before_merge_filter <- transformed_clean_before_merge %>%
    filter(!is.na(lab)) %>%
    tidyr::gather(key=DIM, value=score, 5:24)
  dim(transformed_clean_before_merge_filter)
  result2 <- transformed_clean_before_merge_filter %>%
    group_by(as.factor(DIM)) %>%
    do(tidy(t.test(score ~ MISSING, data = ., 
                   alt = "two.sided", 
                   paired = FALSE, 
                   conf.level = 0.99)))
  result2 <- result2[order(result2$`as.factor(DIM)`), ]
  #result2$pattern <- result2$estimate/abs(result2$estimate) #option
  result2$Category <- "holdout case"
  
  #all data(no simulation of holdouts)
  transformed_clean_before <- labData_gs_before_after_transformed_clean %>%
    dplyr::select(ID, lab, XTIME) %>%
    filter(XTIME == 2)
  transformed_clean_before$MISSING <- as.numeric(is.na(transformed_clean_before[, 2]))
  transformed_clean_before_merge <- merge(transformed_clean_before, PC, by = "ID")
  transformed_clean_before_merge_filter <- transformed_clean_before_merge %>%
    tidyr::gather(key=DIM, value=score, 5:24)
  print(transformed_clean_before_merge_filter)
  result1 <- transformed_clean_before_merge_filter %>%
    group_by(as.factor(DIM)) %>%
    do(tidy(t.test(score ~ MISSING, data = ., 
                   alt = "two.sided", 
                   paired = FALSE, 
                   conf.level = 0.99)))
  #result1$pattern <- -result1$estimate/abs(result1$estimate)
  result1 <- result1[order(result1$`as.factor(DIM)`), ]
  result1$Category <- "all data"
  print(result1)
  result_final <- rbind(result1, result2)
  result_final <- data.frame(rbind(result_final, result))
  p1 <- ggplot(result_final) + 
    geom_point(aes(x = as.factor.DIM., y = -log10(p.value), color = as.factor(Category), size = estimate), fill = "white") +
    ylab(paste("-log10(p.value) for ", lab))
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.background = element_rect(fill = "lightblue"),
          panel.grid.major = element_line(size = 0.8, linetype = 'solid', colour = "grey"),
          panel.border = element_blank()) + 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14)) +
    theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"))
  png(paste('path_to_directory', lab, '.png', sep=''), units="in", width=10, height=8, res=300)
  plot(p1)
  dev.off()
}


# linear mixed effect model to predict lifetime diabetics  ----------------

library(lme4)

head(imputedData) #for short format
head(imputedDataRepeat) #for long format

imputedDataRepeat_50_shortformat <- imputedDataRepeat %>%
  dplyr::select(ID, XTIME, labs, score, repeatnum) %>%
  spread(labs, score) %>%
  dplyr::select(sort(names(.)))
colnames(imputedDataRepeat_50_shortformat)[which(names(imputedDataRepeat_50_shortformat) == "ID")] <- "ID_INDEX"
head(imputedDataRepeat_50_shortformat)
unique(imputedDataRepeat_50_shortformat$ID_INDEX)
length(unique(imputedDataRepeat_50_shortformat$ID_INDEX))
labs <- colnames(imputedDataRepeat_50_shortformat[,4:ncol(imputedDataRepeat_50_shortformat)])

#holdout values
transformed_clean_holdout <- labData_gs_before_after_transformed_clean_holdout %>%
  dplyr::select(ID, X17856.6) #select HgA1c laboratory variable
ID <- transformed_clean_holdout$ID  #no need to remove "PT" prefix

#holdout cases
transformed_clean_holdout <- labData_gs_before_after_transformed_clean_holdout_case %>%
  dplyr::select(ID, X17856.6) 
ID <- transformed_clean_holdout$ID 

#for longformat imputation files imputedDataRepeat
summary <- data.frame()
summaryData <- list()
completeData <- data.frame()
completeDataRepeat <- data.frame()
imputedDataRepeat_add <- data.frame()
Result_glmer = function(imputedDataRepeat, ID) {
  for (j in 1:50) {
    print(j)
    imputedDataRepeat_add <- imputedDataRepeat[imputedDataRepeat$repeatnum == j, ]
    imputedDataRepeat_add <- imputedDataRepeat_add %>%
      arrange(XTIME, ID_INDEX) #change to XTIME for holdout_case
    imputedDataRepeat_add$ID <- ID
    imputedDataRepeat_add <- data.frame(imputedDataRepeat_add[, c("ID", "XTIME", "X17856.6", "repeatnum")])
    features <- c("DIABETES_ALL_TIME", "DIABETES_AT_INDEX", "DIABETES_PREINDEX")
    for (i in 1:length(features)){
      print(i)
      imputedDataRepeat_add_lab <- merge(imputedDataRepeat_add, GNSIS_binary[, c("ID", "PT_SEX", features[i])], by = "ID")
      mod <- glmer(imputedDataRepeat_add_lab[, 6] ~ (1 | imputedDataRepeat_add_lab$ID) + imputedDataRepeat_add_lab$XTIME + as.factor(imputedDataRepeat_add_lab$PT_SEX) + scale(imputedDataRepeat_add_lab$X17856.6), family = binomial, na.action=na.omit, data = imputedDataRepeat_add_lab)
      summary <- data.frame(tidy(mod, conf.int=TRUE,exponentiate=TRUE,effects="fixed"))
      summary$features <- features[i]
      summaryData[[i]] <- summary
    }
    completeData <- as.data.frame(do.call(rbind,summaryData))
    completeData$repeatnum <- j
    completeDataRepeat <- rbind(completeDataRepeat, completeData)
  }
  write.csv(completeDataRepeat, file="path_to_file.csv", col.names = T, row.names = F, quote = F)
  return(completeDataRepeat)
}

Result <- Result_glmer(imputedDataRepeat_50_shortformat, ID)
head(Result)
Result_select <- Result[Result$term == 'scale(imputedDataRepeat_add_lab$X17856.6)', ]
head(Result_select)

##fit pooling for longformat
summary <- data.frame()
summaryData <- list()
completeData <- data.frame()
completeDataRepeat <- data.frame()
imputedDataRepeat_add <- data.frame()
mids_creator = function(imputedDataRepeat, ID) {
  for (j in 1:50) {
    print(j)
    imputedDataRepeat_add <- imputedDataRepeat[imputedDataRepeat$repeatnum == j, ]
    imputedDataRepeat_add <- imputedDataRepeat_add %>%
      arrange(XTIME, ID_INDEX) #change to XTIME for hold_case
    imputedDataRepeat_add$ID <- ID
    imputedDataRepeat_add <- data.frame(imputedDataRepeat_add[, c("ID", "XTIME", "X17856.6", "repeatnum")])
    features <- c("DIABETES_ALL_TIME", "DIABETES_AT_INDEX", "DIABETES_PREINDEX")
    imputedDataRepeat_add_lab <- merge(imputedDataRepeat_add, GNSIS_binary[, c("ID", "PT_SEX", features)], by = "ID")
    summaryData[[j]] <- imputedDataRepeat_add_lab
  }
  return(summaryData)
}

imp <- mids_creator(imputedDataRepeat_50_shortformat, ID)

library(purrr)
fit5_holdoutcase <- imp %>%
  purrr::map(lme4::glmer,
             formula = DIABETES_ALL_TIME ~ (1 | ID) + XTIME + as.factor(PT_SEX) + scale(X17856.6),
             family = binomial) %>%
  pool() %>%
  summary(conf.int = TRUE, exponentiate = TRUE)

fit05_holdoutcase1 <- imp %>%
  purrr::map(lme4::glmer,
             formula = DIABETES_ALL_TIME ~ (1 | ID) + XTIME + as.factor(PT_SEX) + scale(X17856.6),
             family = binomial) %>%
  pool()

#holdout values
pool_result <- read.csv("path_to_poolingresult_file.csv", header =T , stringsAsFactors = F)
head(pool_result)
sort(unique(pool_result$features))
pool_result_select <- pool_result[pool_result$term == "scale(imputedDataRepeat_add_lab$X17856.6)" & pool_result$features == "DIABETES_ALL_TIME", ]
#holdout cases
pool_result1 <- read.csv("path_to_poolingresult_file.csv", header =T , stringsAsFactors = F)
head(pool_result1)
sort(unique(pool_result1$features))
pool_result1_select <- pool_result1[pool_result1$term == "scale(imputedDataRepeat_add_lab$X17856.6)" & pool_result1$features == "DIABETES_ALL_TIME", ]

colnames(pool_result_select)
#[1] "effect"    "term"      "estimate"  "std.error" "statistic" "p.value"   "conf.low"  "conf.high" "features"  "repeatnum"

#within and between imputation variance
within <- mean(pool_result_select$std.error)  #0.328
grandm <- mean(pool_result_select$estimate)
between <- sum((pool_result_select$estimate - grandm)^2)/49   # 0.006594

grandvar <- within + ((1 + (1/50)) * between)

#comparing the effect size (odds ratio) between imputation algorithms
t.test(pool_result_select$estimate, pool_result1_select$estimate, alternative = "two.sided", var.equal = TRUE)

