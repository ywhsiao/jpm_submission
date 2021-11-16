# merge all annotation; including gnomAD v2, revel, clinvar and ExAC(no results)
# data visualization 

library(data.table)
library(dplyr)
setwd('public_data_2/ref/gnomad/')
temp <-  list.files(pattern="*.csv")
temp_name <- gsub(".csv", "", temp)
for (i in 1:length(temp)) assign(temp_name[i], fread(temp[i], header = TRUE))
# merge to a big matrix and save it for downstream analysis
GloEnv <- ls(pattern = "gnomad_v2_")
df_lst <- list()
for (i in GloEnv){
  df_lst[[i]] <- get(i)
}

final_data <- do.call(bind_rows, df_lst)
dim(final_data)
fwrite(final_data, 'gnomad_v2.csv', sep = ",")

setwd('public_data_2/')
#import the most sig variants
sig_variant <- read.csv('snp-based/04_groupwise_tests/all_skato_rift_outlier.csv', header = T, stringsAsFactors = F)
names(sig_variant)[2] <- 'SNP_detail'
sig_variant <- sig_variant[,c(2,12:13)]
snp2gene <- read.csv('snp-based/04_groupwise_tests/snp2gene.csv', header = T, stringsAsFactors = F)
names(snp2gene)
snp2gene <- snp2gene[,c(1:3, 6,17,18)]
names(snp2gene) <- c('SNP_detail', 'tumor_sample', 'Hugo_Symbol', 'SNP_info', 'race', 'cancer')
sig_variant_detail <- snp2gene %>% right_join(sig_variant, by=c('SNP_detail', 'race', 'cancer'))
variantList <- unique(sig_variant_detail$SNP_detail)

# query database using dbplyr
# gnomad v2 
setwd('public_data_2/ref/gnomad/')
library(dplyr)
library(dbplyr)
library(DBI)
gnomad <- dbConnect(RSQLite::SQLite(), dbname= 'gnomad_v2.db')
dbListTables(gnomad)
refdb <- tbl(gnomad, "gnomad_v2")
gnomad_results <- as.data.frame(refdb %>% filter(SNP_detail%in% variantList))
gnomad_results_AF <- cbind(gnomad_results %>%
  select(starts_with('AF')),gnomad_results[,32:34])
# revel
setwd('public_data_2/ref/revel/')
library(dplyr)
library(dbplyr)
library(DBI)
revel <- dbConnect(RSQLite::SQLite(), dbname= 'revel_hg19.db')
dbListTables(revel)
refdb <- tbl(revel, "revel_hg19")
revel_results <- as.data.frame(refdb %>% filter(SNP_detail %in% variantList))

#ExAC and clinival provided by tcga
library(data.table)
setwd('public_data_2/snp-based/00_inputs/')
annot <- fread('tcga_variant_ExAC_COMIST_annotation.csv')
annot <- annot[,-c(1:3,21:22)]
names(annot)
# merge all 
merge_tmp1 <- merge(sig_variant_detail, gnomad_results_AF, by = 'SNP_detail', all.x = T)
merge_tmp2 <- merge(merge_tmp1, revel_results, by = 'SNP_detail', all.x = T)
merge_all <- merge(merge_tmp2, annot, by = 'SNP_detail', all.x = T)
write.csv(merge_all, '../04_groupwise_tests/sig_variant_with_annot.csv', row.names = F)
names(merge_all)
selected_df <- unique(merge_all[,c(1,3,4,5,6,18,19,29)])
head(selected_df)
write.csv(selected_df, '../04_groupwise_tests/sig_variant_with_annot_selected.csv', row.names = F)

# summary statistics
library(dplyr)
var_counts <- selected_df[,c(1,4,5)]
var_sum <- var_counts%>% count(race, cancer)
names(var_sum)[3] <- '#sig_variants'

gene_counts <- unique(selected_df[,c(2,4,5)])
gene_sum <- gene_counts%>% count(race, cancer)
names(gene_sum)[3] <- '#sig_genes'
all_sum <- cbind(var_sum, gene_sum)[,-4:-5]

gene_groupby <- as.data.frame(aggregate(Hugo_Symbol ~ race+cancer, data = gene_counts, FUN = paste, collapse = ","))
final_df <- cbind(all_sum,gene_groupby)[,-5:-6]
write.csv(final_df, '../04_groupwise_tests/sig_variant_counts_with_genes.csv', row.names = F)

# data visualization for annotation results
setwd('public_data_2/snp-based/04_groupwise_tests/')
library(ggplot2) 
# revel
revel <- merge_all %>% 
  select(c('race', 'REVEL')) %>%
  filter(!is.na(REVEL))
revel$REVEL <- as.numeric(revel$REVEL)
revel$race[revel$race == 'black'] <- 'African American/Black'
revel$race[revel$race == 'white'] <- 'White'
revel$race[revel$race == 'asian'] <- 'Asian'
revel$race <- factor(revel$race, levels = c("White", "African American/Black", "Asian"))
revel_median <- revel %>%
  group_by(race) %>%
  summarise(Median = median(REVEL))
str(revel)
ggp_revel <- ggplot(revel, aes(x = REVEL)) +
  geom_histogram(binwidth = 0.05, color = "grey30", fill = "white") +
  facet_grid(race ~ .) +
  geom_vline(data = revel_median, aes(xintercept = Median), linetype = "dashed", color = "red") +
  theme_bw()
ggsave('revelscore_byrace.png', width = 7, height = 7, dpi = 300)

# clinivar
clinvar <- merge_all %>% 
  select(c('race', 'CLIN_SIG')) %>%
  filter(CLIN_SIG !='.')
clinvar$race[clinvar$race == 'black'] <- 'African American/Black'
clinvar$race[clinvar$race == 'white'] <- 'White'
clinvar$race[clinvar$race == 'asian'] <- 'Asian'
clinvar$race <- factor(clinvar$race, levels = c("White", "African American/Black", "Asian"))
ggp_clinvar <- ggplot(clinvar, aes(x = CLIN_SIG)) +
  geom_bar(color = "grey30", fill = "white", stat="count") +
  facet_grid(race ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
ggsave('clinvar_byrace.png', width = 7, height = 10, dpi = 300)

# figure S6
library(ggpubr)
figure <- ggarrange(ggp_revel, ggp_clinvar,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
ggsave('tcga_pan-cancer_validation_db_figureS6.png',width = 35, height =20, units = 'cm')

# gnomAD (dismiss)
names(merge_all)
gnomAD <- merge_all %>% 
  select(c('race','tumor_sample', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin'))
#gnomAD <- gnomAD[complete.cases(gnomAD), ]
# white
options(scipen = 999)

df <- merge_all %>%
  select(c('tumor_sample', 'race', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin'))
df$race[df$race == 'black'] <- 'African American/Black'
df$race[df$race == 'white'] <- 'White'
df$race[df$race == 'asian'] <- 'Asian'
df_long <- reshape2::melt(df, id.vars = c('tumor_sample', 'race'), variable.name = "type", value.name = "value")
df_long <- df_long[!is.na(df_long$value),]
df_long$value <- as.numeric(df_long$value)
df_long <- df_long[!is.na(df_long$value),]
ggplot(df_long, aes(x=value)) +
  geom_histogram(bins=1000, color = "grey30", fill = "grey30") +
  facet_grid(type ~ race) +
  xlim(min(df_long$value), max(df_long$value))+
  theme_bw()


ggsave('gnomad_byrace.png', width = 7, height = 10, dpi = 300)

