# pan-cancer: preclean clinical data ; merge with maf file
rm(list = ls())
setwd('public_data_2/clinical/')

# sample & clinical information
sample_id <-  read.table('Samples.tsv', sep = '\t', header = F, stringsAsFactors = F)
names(sample_id) <- c('submitter_id', 'cancer')
sample_id$submitter_id <- substring(sample_id$submitter_id,1,12)
head(sample_id)
load('clinical_survival_pancancer_atlas.RData')
clinical <- dat
names(clinical)
clinical_list <- as.data.frame(names(clinical))
length(unique(clinical$bcr_patient_barcode))
pheno_race <- as.data.frame(clinical[,c(2,3,4,5,7,8,10,20,741,742,744,745)])
# remove race and OS NA
pheno_race <- pheno_race[pheno_race$race != '[Not Evaluated]',]
length(unique(pheno_race$bcr_patient_barcode))
pheno_race$overall.survival <- NA
for (i in 1:dim(pheno_race)[1]){
  if (is.na(pheno_race$days_to_death[i])){
    pheno_race$overall.survival[i] <- pheno_race$days_to_last_followup[i]
  }else{
    pheno_race$overall.survival[i] <- pheno_race$days_to_death[i]
  }
}

pheno_race <- pheno_race[!is.na(pheno_race$overall.survival),]
pheno_race <- pheno_race[!is.na(pheno_race$vital_status),]
pheno_race$vital_status <- ifelse(pheno_race$vital_status == "Dead", 1, 0)

# sample id with racial info
sample_id_with_race <- merge(pheno_race, sample_id, by.x = 'bcr_patient_barcode', by.y = 'submitter_id')
dim(sample_id_with_race)
table(sample_id_with_race$race)
# keep black or african american, asian and white
sample_id_with_race <- sample_id_with_race[sample_id_with_race$race%in%c('BLACK OR AFRICAN AMERICAN', 'ASIAN', 'WHITE'),]
dim(sample_id_with_race)
head(sample_id_with_race)
cohort_stat <- as.data.frame.matrix(table(sample_id_with_race$acronym, sample_id_with_race$race))
final_cohort <- cohort_stat[rowSums(cohort_stat) >=100, ]
final_cohort_list <- row.names(final_cohort)
sample_id_with_race <- sample_id_with_race[sample_id_with_race$acronym %in% final_cohort_list, ]
for (i in 1:dim(final_cohort)[1]){
  if (final_cohort$ASIAN[i]>=10){
    final_cohort$ASIAN_index[i] <- 1
  }else{
    final_cohort$ASIAN_index[i] <- 0
  }
  if (final_cohort$`BLACK OR AFRICAN AMERICAN`[i]>=10){
    final_cohort$BLACK_index[i] <- 1
  }else{
    final_cohort$BLACK_index[i] <- 0
  }
  if (final_cohort$WHITE[i]>=10){
    final_cohort$WHITE_index[i] <- 1
  }else{
    final_cohort$WHITE_index[i] <- 0
  }
}
asian_remove <- row.names(final_cohort[final_cohort$ASIAN_index == 0,])
black_remove <- row.names(final_cohort[final_cohort$BLACK_index == 0,])
library(dplyr)
sample_id_with_race_over10 <- sample_id_with_race[!(sample_id_with_race$acronym %in% asian_remove & sample_id_with_race$race == 'ASIAN'),]
sample_id_with_race_over10_1 <- sample_id_with_race_over10[!(sample_id_with_race_over10$acronym %in% black_remove &sample_id_with_race_over10$race == 'BLACK OR AFRICAN AMERICAN'),]
sample_id_with_race <- sample_id_with_race_over10_1
table(sample_id_with_race$acronym)
clinical_race <- sample_id_with_race
write.csv(clinical_race, '../tcag_pan-cancer_clinical_preclean.csv', row.names = F)
dim(clinical_race)

# merge with maf file
setwd('public_data_2/snp-based/')
library(maftools)
maf_matrix = read.maf('00_inputs/mc3.v0.2.8.PUBLIC.maf')
tumor_sample <- maf_matrix@data$Tumor_Sample_Barcode
tumor_sample <- substring(tumor_sample,1,12)
snp_annot <- as.data.frame(cbind(tumor_sample, maf_matrix@data$Hugo_Symbol, maf_matrix@data$Chromosome, maf_matrix@data$Start_Position, maf_matrix@data$End_Position, 
                                 maf_matrix@data$Reference_Allele,maf_matrix@data$Tumor_Seq_Allele2 ,as.character(maf_matrix@data$Variant_Classification), maf_matrix@data$SIFT, maf_matrix@data$PolyPhen))
names(snp_annot) <- c('tumor_sample', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'SIFT', 'PolyPhen')
write.csv(snp_annot, '../tcga_pan-cancer_snp_annot.csv', row.names = F)
length(unique(tumor_sample))
names(snp_annot)
snp_annot_clin_race <- merge(snp_annot, clinical_race, by.x='tumor_sample', by.y='bcr_patient_barcode')

snp_annot_clin_race$SIFT <- gsub("\\(\\)", "", (gsub("\\(\\.\\)", "", (gsub("([0-9]+)", "", snp_annot_clin_race$SIFT)))))
snp_annot_clin_race$PolyPhen <- gsub("\\(\\)", "", (gsub("\\(\\.\\)", "", (gsub("([0-9]+)", "",snp_annot_clin_race$PolyPhen)))))

write.csv(snp_annot_clin_race, '../tcga_pan-cancer_snp_annot_clinic.csv', row.names = F)

# sample summary table 
to_count_stats <- unique(snp_annot_clin_race[c("tumor_sample","cancer","race")])
count_stats <- as.data.frame.matrix(table(to_count_stats$cancer,to_count_stats$race))
head(count_stats)
str(count_stats)
count_stats$SUM <- rowSums(count_stats)
count_stats <- rbind(count_stats, Total = colSums(count_stats))
write.csv(count_stats, '../tcga_sample_counts.csv')
