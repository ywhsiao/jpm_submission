# pan-caner_genesets: HRD score (1) data visualization; (2) cox regression with covariates; (3) survival analysis
rm(list = ls())

set_ID <- c('DDRD_assay_42','mutated_gene_21','HR_PARP_132','DDR_276', 'genomewide')
geneset_id = set_ID[1]

# generate a hrd matrix per geneset
mergeHRDresult <- function(geneset_id){
  dir_loc <- paste0('public_data_2/cnv-based/', geneset_id)
  setwd(dir_loc)
  temp <-  list.files(pattern="*_HRDresults.txt")
  temp_name <- gsub("_HRDresults.txt", "", temp)
  for (i in 1:length(temp)) assign(temp_name[i], read.csv(temp[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  # merge to a big matrix and save it for downstream analysis
  GloEnv <- ls(pattern = "TCGA-")
  df_lst <- list()
  for (i in GloEnv){
    df_lst[[i]] <- get(i)
  }
  
  final_data <- do.call(rbind, df_lst)
  write.csv(final_data, paste0('../',geneset_id,'_HRDresults.csv'), row.names = F) 
}
for (i in 1:5){
  mergeHRDresult(set_ID[i])
}
# generate clinic+hrd matrix per geneset
setwd('public_data_2/')
generateClinicPlusHRD <- function(geneset_id){
  snp_clinic <- read.csv('snp-based/tcga_pan-cancer_snp_annot_clinic.csv', header = T, stringsAsFactors = F)
  snp_clinic <- unique(snp_clinic[,c(1,11:13,16:22)])
  snp_clinic$pathologic_stage_combined <- gsub("Stage ","",snp_clinic$pathologic_stage_combined)
  hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_HRDresults.csv'), header = T, stringsAsFactors = F)
  names(hrd_geneset) <- c('tumor_sample','HRD_LOH', 'HRD_TAI', 'HRD_LST', 'HRD_score')
  clinic_hrd_geneset <- merge(snp_clinic, hrd_geneset, by= 'tumor_sample')
  write.csv(clinic_hrd_geneset,paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), row.names = F)
}
for (i in 1:5){
  generateClinicPlusHRD(set_ID[i])
}

# distribution of HRD score
library(ggplot2)
library(Rmisc)
library(grid)
library(gridExtra)
library(ggpubr)
library(tidyr)
drawStackedBarPlot <- function(geneset_id){
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  to_stacked_bar_plot <- clinic_hrd_geneset[,c(1,2,12:14)]
  to_stacked_bar_plot$acronym <- as.factor(to_stacked_bar_plot$acronym)
  to_stacked_bar_plot <- aggregate(to_stacked_bar_plot[,3:5], by=list(to_stacked_bar_plot$acronym), FUN=mean)
  names(to_stacked_bar_plot)[1] <- 'acronym'
  to_stacked_bar_plot <-gather(to_stacked_bar_plot, HRD_types, HRD_values,HRD_LOH, HRD_TAI, HRD_LST, factor_key=TRUE)
  
  ggp <- ggplot(to_stacked_bar_plot, aes(x = acronym, y =HRD_values, fill = HRD_types)) +  # Create stacked bar chart
    geom_bar(stat = "identity")+
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    ggtitle(geneset_id)
  ggp
  ggsave(paste0('cnv-based/', geneset_id, '_hrd_stackedBarPlot.png'),width = 20, height =10, units = 'cm')
}
for (i in 1:5){
  drawStackedBarPlot(set_ID[i])
}
setwd('../')
# draw a big stacked bar plot faceted by cancer type
lst_toStacked <- list()
for (i in 1:5){
  if(i ==2) next
  geneset_id = set_ID[i]
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  to_stacked_bar_plot <- clinic_hrd_geneset[,c(1,2,12:14)]
  to_stacked_bar_plot$acronym <- as.factor(to_stacked_bar_plot$acronym)
  to_stacked_bar_plot <- aggregate(to_stacked_bar_plot[,3:5], by=list(to_stacked_bar_plot$acronym), FUN=mean)
  names(to_stacked_bar_plot)[1] <- 'acronym'
  to_stacked_bar_plot <-gather(to_stacked_bar_plot, HRD_types, HRD_values,HRD_LOH, HRD_TAI, HRD_LST, factor_key=TRUE)
  to_stacked_bar_plot$geneset <- geneset_id 
  lst_toStacked[[i]] <- to_stacked_bar_plot
}
library(dplyr)
df_toStacked <- Reduce(rbind, lst_toStacked)
names(df_toStacked)
df_toStacked_1 <- df_toStacked %>%
  group_by(HRD_types, geneset) %>% 
  summarise_at(vars("HRD_values"), mean)
df_toStacked_1$acronym <- 'Pan-cancer'
df_toStacked_1 <- df_toStacked_1[,c(4,1,2,3)]
df <- as_tibble(rbind(df_toStacked, df_toStacked_1))
df$geneset <- factor(df$geneset, levels=c('genomewide', 'DDR_276', 'HR_PARP_132', 'DDRD_assay_42'))
base_plot <-  ggplot(data = df,
                     mapping = aes(x = geneset,
                                   y = HRD_values,
                                   fill = HRD_types,
                                   group = acronym)
) 
# add the bar layer
b <- base_plot +
  geom_bar(stat="identity")+
  facet_grid(~acronym,scales = "free") +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1),legend.position="none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
b
ggsave('cnv-based/tcga_pan-cancer_hrd-score_stacked-bar-plot.png',width = 35, height =10, units = 'cm')
# racial comparison of HRD score; failed for geneset 2 due to the relatively low HRD score
library(ggpubr)
library(ggsignif)
`%!in%` = Negate(`%in%`)
drawPairwiseBoxPlot <- function(geneset_id, y, interval){
  #geneset_id=set_ID[5]
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  clinic_hrd_geneset <- clinic_hrd_geneset[clinic_hrd_geneset$acronym %!in% c('GBM','PAAD','PRAD', 'SKCM', 'TGCT'),]
  names(clinic_hrd_geneset)
  cancer_id <- unique(clinic_hrd_geneset$acronym)
  cancer_with3race = c('BLCA', 'BRCA', 'CESC', 'COAD', 'HNSC', 'LIHC', 'STAD', 'THCA', 'UCEC')
  df_lst <- list()
  for (i in 1:length(cancer_with3race)){
    #i=2
    my_comparisons <- list(c('WHITE','ASIAN'), c('ASIAN', 'BLACK OR AFRICAN AMERICAN'),
    c('WHITE', 'BLACK OR AFRICAN AMERICAN'))
    df <- clinic_hrd_geneset
    df <- df[df$acronym == cancer_with3race[i],]
    df <- df[order(df$race, decreasing = T), ]
    df_test <- compare_means(HRD_score ~ race, comparisons = my_comparisons, p.adjust.method = "fdr", method='wilcox.test', data = df)
    df_test <- df_test %>% mutate(y.position = c(y,y+interval,y+interval*2))
    df_test$acronym <- cancer_with3race[i]
    df_lst[[i]] <- df_test
  }
  final_data <- do.call(rbind, df_lst)
  df <- clinic_hrd_geneset 
  race <- unique(df$race)
  
  my_comparisons <- list(c('WHITE','ASIAN'), c('ASIAN', 'BLACK OR AFRICAN AMERICAN'),
                         c('WHITE', 'BLACK OR AFRICAN AMERICAN'))
  ggp <- ggboxplot( df, x ="race", y = "HRD_score", 
                    facet.by = "acronym", ylim = c(0, y+interval*4),
                    color = "race", palette = "jco", xlab = FALSE,legend = "right") +
    stat_pvalue_manual(final_data, label = "p.adj") +
    stat_compare_means(label.y = y+interval*4-3)
  ggp+rremove("x.text")+rremove("x.ticks")+theme(legend.position="bottom")+ggtitle(geneset_id)
  ggsave(paste0("cnv-based/",geneset_id,"_HRDscores_racial_comparison.png"),width = 30, height =30, units = 'cm')
}
drawPairwiseBoxPlot(set_ID[1], 20, 5)
drawPairwiseBoxPlot(set_ID[3], 20, 5)
drawPairwiseBoxPlot(set_ID[4], 80, 8)
drawPairwiseBoxPlot(set_ID[5], 100,8)

# cox regression with the adjustment of covariates for cohort
library("survival")
library("survminer")
library("Hmisc")
setwd('../')
# run multiple cox regression 
# per cancer
for (i in 1:5){
  if (i==2) next
  geneset_id = set_ID[i]
  print(geneset_id)
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  cancer <- unique(clinic_hrd_geneset$acronym)
  cancer <- cancer[-20]
  coef_list <- list()
  pvalue_list <- list()
  for (j in 1:length(cancer)){
    #j =2
    library(dplyr)
    not_all_na <- function(x) any(!is.na(x))
    print(cancer[j])
    df <- clinic_hrd_geneset[clinic_hrd_geneset$acronym == cancer[j],]
    head(df)
    df$overall.survival <- as.numeric(df$overall.survival)
    df$pathologic_T_combined <- as.factor(df$pathologic_T_combined)
    df$pathologic_M_combined <- as.factor(df$pathologic_M_combined)
    df$pathologic_N_combined <- as.factor(df$pathologic_N_combined)
    df$pathologic_stage_combined <- as.factor(df$pathologic_stage_combined)
    df$gender <- as.factor(df$gender)
    df <- df %>% select(where(not_all_na))
    if(length(unique(df$vital_status))==1) next
    if (dim(df)[2] == 15){
      if (length(unique(df$gender)) == 1){
        res.cox <- coxph(Surv(overall.survival,vital_status) ~ HRD_score+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
      }else{
        res.cox <- coxph(Surv(overall.survival,vital_status) ~ HRD_score+gender+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
      }
    }else{
      if (length(unique(df$gender)) == 1){
        res.cox <- coxph(Surv(overall.survival,vital_status) ~ HRD_score+age_at_initial_pathologic_diagnosis, data =  df)
      }else{
        res.cox <- coxph(Surv(overall.survival,vital_status) ~ HRD_score+gender+age_at_initial_pathologic_diagnosis, data =  df)
      }
    }
    x <- summary(res.cox)
    coefs <- cbind(as.data.frame(x$coefficients), as.data.frame(x$conf.int))
    coefs$cancer <- cancer[j]
    coef_list[[j]] <- coefs
    pvalue <- as.data.frame(rbind(x$logtest, x$sctest, x$waldtest))
    row.names(pvalue) <- c('logtest', 'sctest', 'waldtest')
    pvalue$cancer <- cancer[j]
    pvalue_list[[j]] <- pvalue
    
  }
  coef_df <- Reduce(rbind, coef_list)
  coef_df$geneset <- geneset_id
  write.csv(coef_df, paste0('cnv-based/',geneset_id,'_multi-cox_coef.csv'))
  pvalue_df <- Reduce(rbind, pvalue_list)
  write.csv(pvalue,paste0('cnv-based/',geneset_id,'_multi-cox_pvalue.csv'))
}
getwd()
# all cancer
for (i in 1:5){
  if (i==2) next
  #i=1
  not_all_na <- function(x) any(!is.na(x))
  geneset_id = set_ID[i]
  print(geneset_id)
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  df <- clinic_hrd_geneset
  df$overall.survival <- as.numeric(df$overall.survival)
  df$pathologic_T_combined <- as.factor(df$pathologic_T_combined)
  df$pathologic_M_combined <- as.factor(df$pathologic_M_combined)
  df$pathologic_N_combined <- as.factor(df$pathologic_N_combined)
  df$pathologic_stage_combined <- as.factor(df$pathologic_stage_combined)
  df$gender <- as.factor(df$gender)
  df <- df %>% select(where(not_all_na))
  res.cox <- coxph(Surv(overall.survival,vital_status) ~ HRD_score+gender+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
  x <- summary(res.cox)
  coefs <- cbind(as.data.frame(x$coefficients), as.data.frame(x$conf.int))
  coefs$cancer <- 'Pan-cancer'
  coefs$geneset <- geneset_id
  pvalue <- as.data.frame(rbind(x$logtest, x$sctest, x$waldtest))
  row.names(pvalue) <- c('logtest', 'sctest', 'waldtest')
  write.csv(coefs, paste0('cnv-based/',geneset_id,'_pan-cancer_multi-cox_coef.csv'))
  write.csv(pvalue,paste0('cnv-based/',geneset_id,'_pan-cancer_multi-cox_pvalue.csv'))
}

# plot forest plots
setwd('cnv-based/')
library(data.table)
temp <-  list.files(pattern="*_multi-cox_coef.csv")
temp_name <- gsub("_multi-cox_coef.csv", "", temp)
for (i in 1:length(temp)) assign(temp_name[i], read.csv(temp[i], header = TRUE, stringsAsFactors = FALSE))

df_list <- list()
for (i in 1: length(temp_name)){
  #i=8
  df <- get(temp_name[i])
  head(df)
  df <- df[df[,1] %like% "HRD_score", ]
  df <- df[,c(3,6,9,10,11,12)]
  names(df) <- c('OddsRatio', 'P-value', 'CI-low', 'CI-high', 'cancer', 'geneset')
  df_list[[i]] <- df
}
df <- Reduce(rbind, df_list)
df <- df[!is.infinite(df$`CI-high`),]
df <- df[order(df$cancer, decreasing = F), ]
df_pan <- df[df$cancer == 'Pan-cancer',]
rownames(df_pan)<-1:nrow(df_pan)
df <- df[df$cancer != 'Pan-cancer',]
rownames(df)<-1:nrow(df)
df <- as_tibble(rbind(df, df_pan))

df$cancer <- factor(df$cancer, levels=unique(df$cancer))
df$geneset <- factor(df$geneset, levels=c('genomewide', 'DDR_276', 'HR_PARP_132', 'DDRD_assay_42'))
names(df)[1] <- 'HazardRatio'
library(dplyr)
library(ggplot2)
fp <- df %>%
  ggplot(aes(y=HazardRatio, x=geneset, label=geneset)) +
  geom_point(size=1, shape=20) +
  geom_errorbar(aes(ymin=`CI-low`, ymax=`CI-high`),width =0) +
  geom_hline(yintercept=1, linetype='longdash') +
  facet_grid(~cancer, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(angle = 90),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust=1))
fp
ggsave('tcga_pan-cancer_hrd-score_forest-plot.png',width = 35, height =10, units = 'cm')

# figure1

figure <- ggarrange(b, fp,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
ggsave('tcga_pan-cancer_hrd-score_figure1.png',width = 35, height =20, units = 'cm')

# survival analysis; to deal with the count table
rm(list=ls())
setwd('public_data_2/')
library(tibble)
set_ID <- c('DDRD_assay_42','mutated_gene_21','HR_PARP_132','DDR_276', 'genomewide')
geneset_id = set_ID[1]

#percancer
surv_lst <- list()
for (j in 1:5){
  #j=1
  if (j == 2) next
  geneset_id = set_ID[j]
  print(geneset_id)
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  cancer <- unique(clinic_hrd_geneset$acronym)
  df_lst <- list()
  #label 
  for (i in 1:length(cancer)){
    #i =1
    df <- clinic_hrd_geneset[clinic_hrd_geneset$acronym == cancer[i],]
    cutoff <- quantile(df$HRD_score, 0.667)
    df$HRD_score <- ifelse(df$HRD_score>cutoff, 'HR_deficiency', 'Not_HR_deficiency')
    df_lst[[i]] <- df
  }
  clinic_hrd_geneset_new <- Reduce(rbind, df_lst)
  # check quantile
  clinic_hrd_geneset_new$HRD_score <- as.factor(clinic_hrd_geneset_new$HRD_score)
  clinic_hrd_geneset_new$overall.survival <- clinic_hrd_geneset_new$overall.survival/365
  clinic_hrd_geneset_new$geneset <- geneset_id
  surv_lst[[j]] <- clinic_hrd_geneset_new
}
surv_df <- Reduce(rbind, surv_lst)
surv_df <- surv_df[,c(1,2,4,6,11,15,16)]
surv_df <- surv_df[order(surv_df$acronym, decreasing = F), ]
rownames(surv_df)<-1:nrow(surv_df)
#pan-cancer
library(tidyr)
surv_lst <- list()
for (j in 1:5){
  if (j == 2) next
  geneset_id = set_ID[j]
  print(geneset_id)
  clinic_hrd_geneset <- read.csv(paste0('cnv-based/',geneset_id,'_clinic_plus_hrd.csv'), header = T, stringsAsFactors = F)
  cancer <- unique(clinic_hrd_geneset$acronym)
  #label 
  df <- clinic_hrd_geneset
  cutoff <- quantile(df$HRD_score, 0.667)
  df$HRD_score <- ifelse(df$HRD_score>cutoff, 'HR_deficiency', 'Not_HR_deficiency')
  clinic_hrd_geneset_new <- df
  # check quantile
  clinic_hrd_geneset_new$HRD_score <- as.factor(clinic_hrd_geneset_new$HRD_score)
  clinic_hrd_geneset_new$overall.survival <- clinic_hrd_geneset_new$overall.survival/365
  clinic_hrd_geneset_new$geneset <- geneset_id
  clinic_hrd_geneset_new$acronym <- 'Pan-cancer'
  surv_lst[[j]] <- clinic_hrd_geneset_new
}
surv_df_pan <- Reduce(rbind, surv_lst)
surv_df_pan <- surv_df_pan[,c(1,2,4,6,11,15,16)]
surv_df_final <- rbind(surv_df, surv_df_pan)
surv_df_final <- surv_df_final[surv_df_final$acronym !='SKCM',]
toexport <- surv_df_final %>%
  dplyr::count(var1 =geneset, var2 = acronym, var3=HRD_score)
data_wide <- spread(toexport, var3, n)
data_wide$label <- paste0('HRD:', data_wide$HR_deficiency, '; Not:', data_wide$Not_HR_deficiency)
data_wide <- data_wide[,c(1,2,5)]
names(data_wide)[1:2] <- c('geneset', 'acronym')
# select cancer cohorts to present
hrdFocus <- c('BRCA', 'OV', 'PAAD' ,'PRAD','Pan-cancer')
surv_df_final <- surv_df_final[surv_df_final$acronym %in% hrdFocus,]
surv_df_final$overall.survival <- as.numeric(surv_df_final$overall.survival)
surv_df_final$HRD_score <- as.factor(surv_df_final$HRD_score)
df <- as_tibble(surv_df_final)
df$geneset <- factor(df$geneset, levels=c('genomewide', 'DDR_276', 'HR_PARP_132', 'DDRD_assay_42'))
df$acronym <- factor(df$acronym, levels=c('BRCA', 'OV', 'PAAD','PRAD','Pan-cancer'))
fit1 <-  survfit(Surv(overall.survival, vital_status) ~HRD_score, data = df)
names(surv_df)
ggsurv <- ggsurvplot(fit1, df, facet.by = c("geneset","acronym"),
                     palette = "jco", pval = TRUE, pval.size = 30)
ggsurv+ geom_text(
  data    = data_wide[data_wide$acronym %in% hrdFocus, ],
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.7,
  vjust   = -20,
  size = 4)
ggsave("cnv-based/tcga_pan-cancer_hrd-score_survival_curve.png",width = 35, height = 30, units = 'cm')

# Fit a Cox proportional hazards model 
rm(list = ls())
setwd('public_data_2/')
hrd_clinc <- read.csv('cnv-based/genomewide_clinic_plus_hrd.csv', header = T, stringsAsFactors = F)
clinc_pre <- read.csv('tcag_pan-cancer_clinical_preclean.csv', header = T, stringsAsFactors = F)
smp <- read.csv('smp_id.txt', sep= ' ', header = F, stringsAsFactors = F)
smp$V1 <- gsub('\\*','',smp$V1)
clinc_pre <- merge(clinc_pre, smp, by.x = 'bcr_patient_barcode', by.y = 'V1', all.y = T)
clinc_pre <- clinc_pre[,c(1,2,4,13)]
hrd_clinc <- merge(clinc_pre, hrd_clinc, by.x = 'bcr_patient_barcode', by.y = 'tumor_sample', all.x = T)
hrd_clinc <- hrd_clinc[,c(1:4,9,18)]
hrd_clinc$HRD_score[is.na(hrd_clinc$HRD_score)] <- 0
names(hrd_clinc) <- c('tumor_sample', 'acronym', 'vital_status', 'overall.survival', 'race', 'HRD_score')
cancer <- unique(hrd_clinc$acronym)
df_lst <- list()
#each cancer
for (i in 1:length(cancer)){
  #i =1
  df <- hrd_clinc[hrd_clinc$acronym == cancer[i],]
  cutoff <- quantile(df$HRD_score, 0.667)
  df$HRD_score <- ifelse(df$HRD_score>cutoff, 'HR_deficiency', 'Not_HR_deficiency')
  df_lst[[i]] <- df
}
clinic_hrd_geneset_new <- Reduce(rbind, df_lst)
# check quantile
clinic_hrd_geneset_new$HRD_score <- as.factor(clinic_hrd_geneset_new$HRD_score)
clinic_hrd_geneset_new$overall.survival <- clinic_hrd_geneset_new$overall.survival/365
head(clinic_hrd_geneset_new)
surv_df_final <- clinic_hrd_geneset_new
surv_df_final$overall.survival <- as.numeric(surv_df_final$overall.survival)
surv_df_final$HRD_score <- as.factor(surv_df_final$HRD_score)
length(unique(surv_df_final$tumor_sample))
# pancancer 
df <- hrd_clinc
cutoff <- quantile(df$HRD_score, 0.667)
df$HRD_score <- ifelse(df$HRD_score>cutoff, 'HR_deficiency', 'Not_HR_deficiency')
clinic_hrd_geneset_new <- df
# check quantile
clinic_hrd_geneset_new$HRD_score <- as.factor(clinic_hrd_geneset_new$HRD_score)
clinic_hrd_geneset_new$overall.survival <- clinic_hrd_geneset_new$overall.survival/365
clinic_hrd_geneset_new$acronym <- 'Pan-cancer'
length(unique(clinic_hrd_geneset_new$tumor_sample))
df_tmp <- rbind(surv_df_final, clinic_hrd_geneset_new)
length(unique(df_tmp$tumor_sample))
# merge each cancer and pancancer:fit cox model: cancer types
df <- as_tibble(df_tmp)
length(unique(df$tumor_sample))
df$acronym = relevel(factor(df$acronym), ref = "Pan-cancer")
surv_object <- Surv(df$overall.survival, df$vital_status)
fit.coxph <- coxph(surv_object ~ acronym, data = df)
g1 <- ggforest(fit.coxph, data = df, main='')

# not duplicated: fit cox model: race
df <- as_tibble(surv_df_final)
df$race[df$race == 'BLACK OR AFRICAN AMERICAN'] <- 'AFRICAN AMERICAN/BLACK'
df$race = relevel(factor(df$race), ref = "WHITE")
df$HRD_score = relevel(factor(df$HRD_score), ref = "Not_HR_deficiency")
surv_object <- Surv(df$overall.survival, df$vital_status)
fit.coxph <- coxph(surv_object ~ race + HRD_score, data = df)
g2 <- ggforest(fit.coxph, data = df, main='')

ggarrange(g1, g2, ncol = 2, nrow = 1)
ggsave("cnv-based/tcga_pan-cancer_hrd-score_cox_forestpolt.png",width = 50, height = 20, units = 'cm')

