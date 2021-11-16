# global HRD vs global TMB at cancer level:(1) correlation; (2) cindex
rm(list = ls())
setwd('public_data_2/')
# get sample list 
smp_list <- read.table('smp_id.txt', header = F, stringsAsFactors = F)
dim(smp_list)
smp_list$V1 <- gsub('\\*','', smp_list$V1)
names(smp_list) <- c("tumor_sample", "acronym")
# get global TMB
global_tmb <- read.csv('tcga_pan-cancer_snp_annot.csv', header = T, stringsAsFactors = F) #tcga_pan-cancer_snp_annot.csv
genelist <- unique(global_tmb$Hugo_Symbol)
global_tmb$count <- 1
global_tmb <- global_tmb[,c(1,11)]
global_tmb_sum <- aggregate(global_tmb$count, by=list(Category=global_tmb$tumor_sample), FUN=sum)
names(global_tmb_sum) <- c('tumor_sample', 'TMB_value')
global_tmb_sum <- merge(global_tmb_sum, smp_list, by = 'tumor_sample')
dim(global_tmb_sum)
head(global_tmb_sum)

# calculate the total length of whole gene
library("biomaRt")
#setno <- set1
hgnc_list <- genelist
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol", values=hgnc_list, mart=human)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position
gene_coords <- gene_coords[-c(31:37,41,43),]
total_length <- sum(as.numeric(gene_coords$size))

# get global HRD 
global_hrd <- read.csv('cnv-based/genomewide_clinic_plus_hrd.csv', header = T, stringsAsFactors = F)
names(global_hrd)
global_hrd <- global_hrd[,c(1,2,15)]
global_hrd_tmb <- merge(global_hrd, global_tmb_sum, by = 'tumor_sample', all.y = T)
dim(global_hrd_tmb)
head(global_hrd_tmb)
global_hrd_tmb <- global_hrd_tmb[,-2]
names(global_hrd_tmb)[4] <- 'acronym'

global_hrd_tmb$HRD_score[is.na(global_hrd_tmb$HRD_score)] <- 0
global_hrd_tmb_pan <- global_hrd_tmb
global_hrd_tmb_pan$acronym <- 'Pan-cancer'
summary(global_hrd_tmb)
global_hrd_tmb_final <- rbind(global_hrd_tmb, global_hrd_tmb_pan)
global_hrd_tmb_final$TMB_value <- global_hrd_tmb_final$TMB_value/total_length*1000000
global_hrd_tmb_final$acronym <- factor(global_hrd_tmb_final$acronym, levels =unique(global_hrd_tmb_final$acronym))
final_df <- global_hrd_tmb_final
names(final_df)
library(ggpubr)
require(pals)
options(scipen = 100)
p <- ggscatter(final_df, x = "TMB_value", y = "HRD_score", size = 0.3, palette = "stepped",
               add = "reg.line", facet.by= "acronym",conf.int = TRUE, show.legend.text = FALSE) +
  stat_cor(method = "spearman", label.y = 90)+
  rremove("legend")
ggpar(p, ylim = c(0, 100))
ggsave("tcga_pan-cancer_hrd_vs_tmb.png",width = 15, height = 10, dpi = 300)

# cindex using cox regression
clinic <- read.csv('tcag_pan-cancer_clinical_preclean.csv', header = T, stringsAsFactors = F)
names(clinic)[1] <- 'tumor_sample'
clinic <- clinic[clinic$tumor_sample %in% smp_list$tumor_sample,]
hrd_tmb_clinic <- merge(clinic, global_hrd_tmb, by='tumor_sample', all.y = T)
dim(hrd_tmb_clinic)
hrd_tmb_clinic <- hrd_tmb_clinic[,c(-5,-6,-14,-17)]
names(hrd_tmb_clinic)[2] <- 'acronym'
hrd_tmb_clinic_pan <- hrd_tmb_clinic
hrd_tmb_clinic_pan$acronym <- 'Pan-cancer'
hrd_tmb_clinic_final <- rbind(hrd_tmb_clinic, hrd_tmb_clinic_pan)
hrd_tmb_clinic_final$TMB_value <- hrd_tmb_clinic_final$TMB_value/total_length*1000000
head(hrd_tmb_clinic_final)
library(survival)
library(survminer)
library(Hmisc)
cancer <- unique(hrd_tmb_clinic_final$acronym)
cancer <- cancer[-20]
coef_list <- list()
pvalue_list <- list()
for (j in 1:length(cancer)){
  #j=1
  library(dplyr)
  not_all_na <- function(x) any(!is.na(x))
  print(cancer[j])
  df <- hrd_tmb_clinic_final[hrd_tmb_clinic_final$acronym == cancer[j],]
  df$overall.survival <- as.numeric(df$overall.survival)
  df$pathologic_T_combined <- as.factor(df$pathologic_T_combined)
  df$pathologic_M_combined <- as.factor(df$pathologic_M_combined)
  df$pathologic_N_combined <- as.factor(df$pathologic_N_combined)
  df$pathologic_stage_combined <- as.factor(df$pathologic_stage_combined)
  df$gender <- as.factor(df$gender)
  df <- df %>% dplyr::select(where(not_all_na))
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
  cindex <- signif(1-rcorr.cens(predict(res.cox),Surv(df$overall.survival,df$vital_status)) [[1]], digits=3)
  x <- summary(res.cox)
  coefs <- cbind(as.data.frame(x$coefficients), as.data.frame(x$conf.int))
  coefs$cancer <- cancer[j]
  coefs$cindex <- cindex
  coef_list[[j]] <- coefs
  pvalue <- as.data.frame(rbind(x$logtest, x$sctest, x$waldtest))
  row.names(pvalue) <- c('logtest', 'sctest', 'waldtest')
  pvalue$cancer <- cancer[j]
  pvalue_list[[j]] <- pvalue
  
}
coef_df_hrd <- Reduce(rbind, coef_list)
coef_df_hrd$terms <- 'HRD'
coef_df_hrd$coxtypes <- 'uni'
coef_list <- list()
pvalue_list <- list()
for (j in 1:length(cancer)){
  #j =2
  library(dplyr)
  not_all_na <- function(x) any(!is.na(x))
  print(cancer[j])
  df <- hrd_tmb_clinic_final[hrd_tmb_clinic_final$acronym == cancer[j],]
  head(df)
  df$overall.survival <- as.numeric(df$overall.survival)
  df$pathologic_T_combined <- as.factor(df$pathologic_T_combined)
  df$pathologic_M_combined <- as.factor(df$pathologic_M_combined)
  df$pathologic_N_combined <- as.factor(df$pathologic_N_combined)
  df$pathologic_stage_combined <- as.factor(df$pathologic_stage_combined)
  df$gender <- as.factor(df$gender)
  df <- df %>% dplyr::select(where(not_all_na))
  if(length(unique(df$vital_status))==1) next
  if (dim(df)[2] == 15){
    if (length(unique(df$gender)) == 1){
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
    }else{
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+gender+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
    }
  }else{
    if (length(unique(df$gender)) == 1){
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+age_at_initial_pathologic_diagnosis, data =  df)
    }else{
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+gender+age_at_initial_pathologic_diagnosis, data =  df)
    }
  }
  cindex <- signif(1-rcorr.cens(predict(res.cox),Surv(df$overall.survival,df$vital_status)) [[1]], digits=3)
  x <- summary(res.cox)
  coefs <- cbind(as.data.frame(x$coefficients), as.data.frame(x$conf.int))
  coefs$cancer <- cancer[j]
  coefs$cindex <- cindex
  coef_list[[j]] <- coefs
  pvalue <- as.data.frame(rbind(x$logtest, x$sctest, x$waldtest))
  row.names(pvalue) <- c('logtest', 'sctest', 'waldtest')
  pvalue$cancer <- cancer[j]
  pvalue_list[[j]] <- pvalue
  
}
coef_df_tmb <- Reduce(rbind, coef_list)
coef_df_tmb$terms <- 'TMB'
coef_df_tmb$coxtypes <- 'uni'
coef_list <- list()
pvalue_list <- list()
for (j in 1:length(cancer)){
  #j =2
  library(dplyr)
  not_all_na <- function(x) any(!is.na(x))
  print(cancer[j])
  df <- hrd_tmb_clinic_final[hrd_tmb_clinic_final$acronym == cancer[j],]
  head(df)
  df$overall.survival <- as.numeric(df$overall.survival)
  df$pathologic_T_combined <- as.factor(df$pathologic_T_combined)
  df$pathologic_M_combined <- as.factor(df$pathologic_M_combined)
  df$pathologic_N_combined <- as.factor(df$pathologic_N_combined)
  df$pathologic_stage_combined <- as.factor(df$pathologic_stage_combined)
  df$gender <- as.factor(df$gender)
  df <- df %>% dplyr::select(where(not_all_na))
  if(length(unique(df$vital_status))==1) next
  if (dim(df)[2] == 15){
    if (length(unique(df$gender)) == 1){
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+HRD_score+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
    }else{
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+HRD_score+gender+age_at_initial_pathologic_diagnosis+pathologic_T_combined+pathologic_M_combined+pathologic_N_combined+pathologic_stage_combined, data =  df)
    }
  }else{
    if (length(unique(df$gender)) == 1){
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+HRD_score+age_at_initial_pathologic_diagnosis, data =  df)
    }else{
      res.cox <- coxph(Surv(overall.survival,vital_status) ~ TMB_value+HRD_score+gender+age_at_initial_pathologic_diagnosis, data =  df)
    }
  }
  cindex <- signif(1-rcorr.cens(predict(res.cox),Surv(df$overall.survival,df$vital_status)) [[1]], digits=3)
  x <- summary(res.cox)
  coefs <- cbind(as.data.frame(x$coefficients), as.data.frame(x$conf.int))
  coefs$cancer <- cancer[j]
  coefs$cindex <- cindex
  coef_list[[j]] <- coefs
  pvalue <- as.data.frame(rbind(x$logtest, x$sctest, x$waldtest))
  row.names(pvalue) <- c('logtest', 'sctest', 'waldtest')
  pvalue$cancer <- cancer[j]
  pvalue_list[[j]] <- pvalue
  
}
coef_df_two <- Reduce(rbind, coef_list)
coef_df_two$terms <- 'HRD+TMB'
coef_df_two$coxtypes <- 'multi'

library(data.table)
coefs_final <- rbind(coef_df_hrd,coef_df_tmb,coef_df_two)
names(coefs_final)
coefs_final <- unique(coefs_final[,10:13])
write.csv(coefs_final, 'tcga_hrd_vs_tmb_cox_results.csv', row.names = F)

# cindex data visualization
coefs_final <- read.csv('tcga_hrd_vs_tmb_cox_results.csv', header = T, stringsAsFactors = F)
str(coefs_final)
cancer_order <- c(sort(unique(coefs_final$cancer)[-24]),'Pan-cancer')

coefs_final$cancer <- factor(coefs_final$cancer, levels = cancer_order)
library(ggplot2)
g <- ggplot(coefs_final, aes(cindex,terms, colour = terms)) +
  geom_point() +
  facet_grid(~cancer) +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90),axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="bottom") 
g
ggsave("tcga_pan-cancer_hrd_vs_tmb_cox.png",width = 10, height = 3, dpi = 300)

# figure S2
figure <- ggarrange(p, g,
                    labels = c("A", "B"),
                    hjust = 0.0,
                    vjust = 1.1,
                    ncol = 1, 
                    nrow = 2,
                    heights = c(2,1),
                    font.label = list(size = 18, color = "black", face = "bold", family = NULL))
ggsave('tcga_pan-cancer_hrd_vs_tmb_figureS2.png',width = 35, height = 25, units = 'cm', dpi=700)


# racial comparison including pan-cancer; data
# smp check
smp_list <- read.table('smp_id.txt', header = F, stringsAsFactors = F)
dim(smp_list)
smp_list$V1 <- gsub('\\*','', smp_list$V1)
names(smp_list) <- c("tumor_sample", "acronym")
# clinical data
clinic <- read.csv('tcag_pan-cancer_clinical_preclean.csv', header = T, stringsAsFactors = F)
names(clinic)[1] <- 'tumor_sample'
clinic <- clinic[clinic$tumor_sample %in% smp_list$tumor_sample,]
dim(clinic)
names(clinic)
clinic <- clinic[, c(1,2,8)]
# global hrd 
global_hrd <- read.csv('cnv-based/genomewide_clinic_plus_hrd.csv', header = T, stringsAsFactors = F)
global_hrd <- global_hrd[,c(1,2,15)]
global_hrd <- merge(clinic, global_hrd, by = 'tumor_sample', all.x = T)
global_hrd <- global_hrd[,-4]
names(global_hrd)[2] <- "acronym"
global_hrd$HRD_score[is.na(global_hrd$HRD_score)] <- 0
global_hrd_pan <- global_hrd
global_hrd_pan$acronym <- 'Pan-cancer'
global_hrd_final <- rbind(global_hrd, global_hrd_pan)
global_hrd_final$race[global_hrd_final$race == 'WHITE'] <- 'White'
global_hrd_final$race[global_hrd_final$race == 'ASIAN'] <- 'Asian'
global_hrd_final$race[global_hrd_final$race == 'BLACK OR AFRICAN AMERICAN'] <- 'African American/Black'

library(ggpubr)
library(ggsignif)
`%!in%` = Negate(`%in%`) 
global_hrd_final <- global_hrd_final[global_hrd_final$acronym %!in% c('GBM','PAAD','PRAD', 'SKCM', 'TGCT'),]
cancer_id <- unique(global_hrd_final$acronym)
cancer_with3race = c('BLCA', 'BRCA', 'CESC', 'COAD', 'HNSC', 'LIHC', 'STAD', 'THCA', 'UCEC', 'Pan-cancer')
df_lst <- list()
y =100
interval =8
for (i in 1:length(cancer_with3race)){
  #i=2
  my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                         c('White', 'African American/Black'))
  df <- global_hrd_final
  df <- df[df$acronym == cancer_with3race[i],]
  df <- df[order(df$race, decreasing = T), ]
  df_test <- compare_means(HRD_score ~ race, comparisons = my_comparisons, p.adjust.method = "fdr", method='wilcox.test', data = df)
  df_test <- df_test %>% mutate(y.position = c(y,y+interval,y+interval*2))
  df_test$acronym <- cancer_with3race[i]
  df_lst[[i]] <- df_test
}
final_data <- do.call(rbind, df_lst)
df <- global_hrd_final
race <- unique(df$race)
cancer_order <- c(sort(unique(df$acronym))[-14], 'Pan-cancer')
df$acronym <- factor(df$acronym, levels = cancer_order)

# for sig comparisons
final_data_sig <- final_data[final_data$acronym %in% c('BLCA', 'BRCA', 'HNSC', 'KIRP', 'LIHC','Pan-cancer', 'UCEC'),]
df_sig <- df[df$acronym %in% c('BLCA', 'BRCA', 'HNSC', 'KIRP', 'LIHC','Pan-cancer', 'UCEC'),]
my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                       c('White', 'African American/Black'))
ggp_hrd_sig <- ggboxplot( df_sig, x ="race", y = "HRD_score", 
                  facet.by = "acronym", ylim = c(0, y+interval*4),
                  color = "race", palette = "jco", xlab = FALSE,legend = "right") +
  stat_pvalue_manual(final_data_sig, label = "p.adj") +
  stat_compare_means(label.y = y+interval*4-3) +rremove("x.text")+rremove("x.ticks")+theme(legend.position="top")
ggp_hrd_sig
ggsave(paste0("tcga_pan-cancer_HRDscores_racial_comparison_sig.png"),width = 25, height =30, units = 'cm', dpi = 700)

# for non_sig comparisons
final_data_nonsig <- final_data[final_data$acronym %!in% c('BLCA', 'BRCA', 'HNSC', 'KIRP', 'LIHC','Pan-cancer', 'UCEC'),]
df_nonsig <- df[df$acronym %!in% c('BLCA', 'BRCA', 'HNSC', 'KIRP', 'LIHC','Pan-cancer', 'UCEC'),]
my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                       c('White', 'African American/Black'))
ggp_hrd_nonsig <- ggboxplot( df_nonsig, x ="race", y = "HRD_score", 
                          facet.by = "acronym", ylim = c(0, y+interval*4),
                          color = "race", palette = "jco", xlab = FALSE,legend = "right") +
  stat_pvalue_manual(final_data_nonsig, label = "p.adj") +
  stat_compare_means(label.y = y+interval*4-3) +rremove("x.text")+rremove("x.ticks")+theme(legend.position="top")
ggp_hrd_nonsig
ggsave(paste0("tcga_pan-cancer_HRDscores_racial_comparison_nonsig.png"),width = 25, height =30, units = 'cm', dpi = 700)


# global tmb
global_tmb <- read.csv('tcga_pan-cancer_snp_annot.csv', header = T, stringsAsFactors = F)
genelist <- unique(global_tmb$Hugo_Symbol)
global_tmb$count <- 1
global_tmb <- global_tmb[,c(1,11)]
global_tmb_sum <- aggregate(global_tmb$count, by=list(Category=global_tmb$tumor_sample), FUN=sum)
names(global_tmb_sum) <- c('tumor_sample', 'TMB_value')
global_tmb_sum <- merge(clinic, global_tmb_sum, by = 'tumor_sample', )
global_tmb_sum$TMB_value[is.na(global_tmb_sum$TMB_value)] <- 0
global_tmb_sum_pan <- global_tmb_sum
global_tmb_sum$acronym <- 'Pan-cancer'
global_tmb_sum_final <- rbind(global_tmb_sum, global_tmb_sum_pan)

# calculate the total length of whole gene
library("biomaRt")
#setno <- set1
hgnc_list <- genelist
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol", values=hgnc_list, mart=human)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position
gene_coords <- gene_coords[-c(31:37,41,43),]
total_length <- sum(as.numeric(gene_coords$size))
#global_tmb_sum$TMB_value <- global_tmb_sum$TMB_value/total_length*1000000
global_tmb_sum_final$TMB_value <- global_tmb_sum_final$TMB_value/total_length*1000000
global_tmb_sum_final$race[global_tmb_sum_final$race == 'WHITE'] <- 'White'
global_tmb_sum_final$race[global_tmb_sum_final$race == 'ASIAN'] <- 'Asian'
global_tmb_sum_final$race[global_tmb_sum_final$race == 'BLACK OR AFRICAN AMERICAN'] <- 'African American/Black'
summary(global_hrd_tmb_final)

library(ggpubr)
library(ggsignif)
`%!in%` = Negate(`%in%`) 
global_tmb_sum_final <- global_tmb_sum_final[global_tmb_sum_final$acronym %!in% c('GBM','PAAD','PRAD', 'SKCM', 'TGCT'),]
cancer_id <- unique(global_tmb_sum_final$acronym)
cancer_with3race = c('BLCA', 'BRCA', 'CESC', 'COAD', 'HNSC', 'LIHC', 'STAD', 'THCA', 'UCEC', 'Pan-cancer')
df_lst <- list()
y =20
interval =3
for (i in 1:length(cancer_with3race)){
  #i=2
  my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                         c('White', 'African American/Black'))
  df <- global_tmb_sum_final
  df <- df[df$acronym == cancer_with3race[i],]
  df <- df[order(df$race, decreasing = T), ]
  df_test <- compare_means(TMB_value ~ race, comparisons = my_comparisons, p.adjust.method = "fdr", method='wilcox.test', data = df)
  df_test <- df_test %>% mutate(y.position = c(y,y+interval,y+interval*2))
  df_test$acronym <- cancer_with3race[i]
  df_lst[[i]] <- df_test
}
final_data <- do.call(rbind, df_lst)
df <- global_tmb_sum_final
race <- unique(df$race)
cancer_order <- c(sort(unique(df$acronym))[-14], 'Pan-cancer')
df$acronym <- factor(df$acronym, levels = cancer_order)

# for sig comparisons
final_data_sig <- final_data[final_data$acronym %in% c('BLCA', 'BRCA', 'CESC', 'ESCA', 'LIHC', 'LUAD'),]
df_sig <- df[df$acronym %in% c('BLCA', 'BRCA', 'CESC', 'ESCA', 'LIHC', 'LUAD'),]
my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                       c('White', 'African American/Black'))
ggp_tmb_sig <- ggboxplot( df_sig, x ="race", y = "TMB_value", 
                  facet.by = "acronym", ylim = c(0, y+interval*4),
                  color = "race", palette = "jco", xlab = FALSE,legend = "right") +
  stat_pvalue_manual(final_data_sig, label = "p.adj") +
  stat_compare_means(label.y = y+interval*4-3)+rremove("x.text")+rremove("x.ticks")+theme(legend.position="top")
ggsave(paste0("tcga_pan-cancer_TMBvalue_racial_comparison_sig.png"),width = 25, height =30, units = 'cm', dpi = 700)

# for non_sig comparisons
final_data_nonsig <- final_data[final_data$acronym %!in% c('BLCA', 'BRCA', 'CESC', 'ESCA', 'LIHC', 'LUAD'),]
df_nonsig <- df[df$acronym %!in% c('BLCA', 'BRCA', 'CESC', 'ESCA', 'LIHC', 'LUAD'),]
my_comparisons <- list(c('White','Asian'), c('Asian', 'African American/Black'),
                       c('White', 'African American/Black'))
ggp_tmb_nonsig <- ggboxplot( df_nonsig, x ="race", y = "TMB_value", 
                          facet.by = "acronym", ylim = c(0, y+interval*4),
                          color = "race", palette = "jco", xlab = FALSE,legend = "right") +
  stat_pvalue_manual(final_data_nonsig, label = "p.adj") +
  stat_compare_means(label.y =y+interval*4-3)+rremove("x.text")+rremove("x.ticks")+theme(legend.position="top")
ggsave(paste0("tcga_pan-cancer_TMBvalue_racial_comparison_nonsig.png"),width = 25, height =30, units = 'cm', dpi = 700)

# figure 3
figure <- ggarrange(ggp_hrd_sig, ggp_tmb_sig,
                    labels = c("A", "B"),
                    common.legend = T,
                    legend = "bottom",
                    ncol = 1, nrow = 2,
                    heights = c(3.8,2),
                    font.label = list(size = 18, color = "black", face = "bold", family = NULL))
ggsave('tcga_pan-cancer_hrd_vs_tmb_figure3.png',width = 25, height = 35, units = 'cm', dpi = 700)

# figure S3
figure <- ggarrange(ggp_hrd_nonsig, ggp_tmb_nonsig,
                    labels = c("A", "B"),
                    common.legend = T,
                    legend = "bottom",
                    ncol = 1, nrow = 2,
                    heights = c(4.3,4.3),
                    font.label = list(size = 18, color = "black", face = "bold", family = NULL))
ggsave('tcga_pan-cancer_hrd_vs_tmb_figureS3.png',width = 25, height = 45, units = 'cm', dpi = 700)
