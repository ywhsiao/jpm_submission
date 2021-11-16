# generate a big reference genelist
rm(list = ls())
setwd('public_data_2/geneset/')
temp <-  list.files(pattern="*_loc_with_genename.txt")
temp_name <- gsub("_loc_with_genename.txt", "", temp)
for (i in 1:length(temp)) assign(temp_name[i], read.csv(temp[i], sep = '\t',header = TRUE, stringsAsFactors = FALSE))

df_list <- list()
for (i in 1: length(temp_name)){
  df <- get(temp_name[i])
  df$geneset <- temp_name[i]
  df_list[[i]] <- df
}
gene_info <- Reduce(rbind, df_list)
gene_info <- gene_info[!duplicated(gene_info$X0),]
head(gene_info)
names(gene_info)[1] <- 'Hugo_Symbol'

# get sample list 
setwd('../')
smp_list <- read.table('smp_id.txt', header = F, stringsAsFactors = F)
dim(smp_list)
smp_list$V1 <- gsub('\\*','', smp_list$V1)
names(smp_list) <- c("tumor_sample", "acronym")

# extract geno & clinic input file from maf file 
library(maftools)
library(dplyr)
maf_matrix <- read.maf('snp-based/00_inputs/mc3.v0.2.8.PUBLIC.maf')
tumor_sample <- maf_matrix@data$Tumor_Sample_Barcode
tumor_sample <- substring(tumor_sample,1,12)
snp_annot <- as.data.frame(cbind(tumor_sample, maf_matrix@data$Hugo_Symbol, maf_matrix@data$Chromosome, maf_matrix@data$Start_Position, maf_matrix@data$End_Position, 
                                 maf_matrix@data$Reference_Allele,maf_matrix@data$Tumor_Seq_Allele2 ,as.character(maf_matrix@data$Variant_Classification), maf_matrix@data$SIFT, maf_matrix@data$PolyPhen, as.character(maf_matrix@data$Variant_Type)))
names(snp_annot) <- c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'SIFT', 'PolyPhen', 'Variant_Type')
selected_maf <- snp_annot[snp_annot$Hugo_Symbol %in% gene_info$Hugo_Symbol,]
selected_maf <- selected_maf[selected_maf$Tumor_Sample_Barcode %in% smp_list$tumor_sample,]
selected_maf$SIFT <- gsub("\\(\\)", "", (gsub("\\(\\.\\)", "", (gsub("([0-9]+)", "", selected_maf$SIFT)))))
selected_maf$PolyPhen <- gsub("\\(\\)", "", (gsub("\\(\\.\\)", "", (gsub("([0-9]+)", "",selected_maf$PolyPhen)))))
selected_maf <-selected_maf[selected_maf$Variant_Classification == "Missense_Mutation" &  selected_maf$SIFT == "deleterious" & selected_maf$PolyPhen %in% c("probably_damaging","possibly_damaging"),]

selected_maf_pan <- selected_maf
selected_maf_pan$acronym <- 'Pan-cancer'
write.table(selected_maf_pan,paste0('snp-based/02_file2skato/maf2skato/Pan-cancer.maf'), quote = F, row.names = F, sep = '\t')
selected_maf_cancer <- merge(smp_list, selected_maf, by.x = 'tumor_sample', by.y = 'Tumor_Sample_Barcode')
cancer_ID <- unique(smp_list$acronym)
for (j in 1:length(cancer_ID)){
  #j=1
  df <- selected_maf_cancer[selected_maf_cancer$acronym== cancer_ID[j],]
  names(df)[1] <- 'Tumor_Sample_Barcode'
  write.table(df,paste0('snp-based/02_file2skato/maf2skato/',cancer_ID[j],'.maf'), quote = F, row.names = F, sep = '\t')
}

# generate clinic input file from clinic file
clinic <- read.csv('tcag_pan-cancer_clinical_preclean.csv', header = T, stringsAsFactors = F)
clinic <- clinic[clinic$bcr_patient_barcode %in% smp_list$tumor_sample,]
clinic <- clinic[,c(1,2,3,7,8:12)]
names(clinic)[1:2] <- c('tumor_sample','cancer')
for (i in 1:length(cancer_ID)){
  tmp <- clinic[clinic$cancer == cancer_ID[i],]
  write.csv(tmp, paste0("snp-based/02_file2skato/clinic2skato/",cancer_ID[i],"_clinic.csv"), row.names = F)
}
clinic$cancer <- 'Pan-cancer'
write.csv(clinic,"snp-based/02_file2skato/clinic2skato/Pan-cancer_clinic.csv", row.names = F)

# run SKAT-O on colab
setwd('public_data_2/')
smp_list <- read.table('smp_id.txt', header = F, stringsAsFactors = F)
dim(smp_list)
smp_list$V1 <- gsub('\\*','', smp_list$V1)
names(smp_list) <- c("tumor_sample", "acronym")

# genotypeMatrix revised function
genotypeMatrix = function(maf, genes = NULL, tsb = NULL, includeSyn = FALSE, vafCol = NULL, vafCutoff = c(0.1, 0.75)){
  
  mdat = subsetMaf(maf = maf, tsb = tsb, genes = genes, fields = vafCol, includeSyn = includeSyn, mafObj = FALSE)
  
  if(nrow(mdat) == 0){
    stop("Zero mutations to make table.")
  }
  
  mdat[,id := paste0(Chromosome, ":", Start_Position, "-", Reference_Allele, '/',Tumor_Seq_Allele2)]
  
  if(!is.null(vafCol)){
    if(vafCol %in% colnames(mdat)){
      colnames(mdat)[which(x = colnames(mdat) == vafCol)] = 't_vaf'
    }else{
      stop(paste0("Column ", vafCol, " not found!"))
    }
    
    if(max(mdat[,t_vaf], na.rm = TRUE) > 1){
      mdat[,t_vaf := as.numeric(as.character(t_vaf))/100]
    }
    
    vafMin = vafCutoff[1]
    vafMax = vafCutoff[2]
    
    mdat.cast = data.table::dcast(data = mdat, formula = id ~ Tumor_Sample_Barcode, value.var = 't_vaf', fill = 0)
    data.table::setDF(x = mdat.cast, rownames = mdat.cast$id)
    mdat.cast = mdat.cast[,-1]
    
    tnumMat = t(mdat.cast) #transposematrix
    mdat.cast = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort
    
    mdat.cast = apply(X = mdat.cast, MARGIN = 2, FUN = function(x){
      ifelse(test = x == 0, yes = "None",
             no = ifelse(test = x > vafMin & x < vafMax, yes = "Het",
                         no = ifelse(test = x > vafMax, yes = "Hom", no = "None")))
    })
  }else{
    mdat.cast = data.table::dcast(data = mdat, formula = id ~ Tumor_Sample_Barcode, value.var = 'id', fill = 0)
    data.table::setDF(x = mdat.cast, rownames = mdat.cast$id)
    mdat.cast = mdat.cast[,-1]
    mdat.cast[mdat.cast != 0] = 1
    
    tnumMat = t(mdat.cast) #transposematrix
    mdat.cast = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort
    
    mdat.cast = apply(X = mdat.cast, MARGIN = 2, FUN = function(x){
      ifelse(test = x == "0", yes = "None",
             no = "Het")
    })
  }
  
  mdat.cast
}

# PC1 for covariate in SKATO
pca_result =read.csv('PCA/ALL_pca.eigenvec', sep = '\t')
names(pca_result)[1] = 'tumor_sample'
pca_result$tumor_sample = substring(pca_result$tumor_sample,1,12)
pca_result = pca_result[,c(1,3)]
length(unique(pca_result$tumor_sample))

library(dplyr) 
cancer_ID <- c(sort(unique(smp_list$acronym)), 'Pan-cancer')
# SKAT-O analysis
for (j in 1:24){
  # convert maf to genotype 
  library(maftools)
  j = 1
  cancer_id = cancer_ID[j]
  print(cancer_id)
  
  maf_file <- read.maf(paste0("snp-based/02_file2skato/maf2skato/",cancer_id,'.maf')) 
  geno <- genotypeMatrix(maf_file, unique(maf_file@data$Hugo_Symbol)) 
  head(geno)# to convert het to 1 whereas none to 0
  geno[geno == 'None'] <- 0
  geno[geno == 'Het'] <- 1
  # to transpose the matrix
  geno <- t(geno)
  # to convert 'Tumor_Sample_Barcode' to 'tumor_sample'
  geno_converted <- cbind(row.names(geno), as.data.frame(geno))
  names(geno_converted)[1] <- 'sampleID'
  # add phenotype information
  library(dplyr)
  pheno_file <- read.csv(file = paste0("snp-based/02_file2skato/clinic2skato/",cancer_id,"_clinic.csv"), header = TRUE, stringsAsFactors = FALSE) # create such file in advance
  table(pheno_file$race)
  pheno_file$pheno_White <- NA
  pheno_file$pheno_White[which(pheno_file$race == "WHITE")] <- '1'
  pheno_file$pheno_White[which(pheno_file$race != "WHITE")] <- '0'
  pheno_file$pheno_Black <- NA
  pheno_file$pheno_Black[which(pheno_file$race == "BLACK OR AFRICAN AMERICAN")] <- '1'
  pheno_file$pheno_Black[which(pheno_file$race != "BLACK OR AFRICAN AMERICAN")] <- '0'
  
  pheno_file$pheno_Asian <- NA
  pheno_file$pheno_Asian[which(pheno_file$race == "ASIAN")] <- '1'
  pheno_file$pheno_Asian[which(pheno_file$race != "ASIAN")] <- '0'
  pheno_file <- pheno_file[,-2]
  head(pheno_file)
  table(pheno_file$pheno_Asian)
  names(pheno_file)[1] <- 'sampleID' 
  pheno_file <- merge(pca_result, pheno_file, by.x='tumor_sample', by.y = 'sampleID', all.y = T)
  names(pheno_file)[1] <- 'sampleID'
  #row.names(pheno_file) <- pheno_file[,1]
  #pheno_file <- pheno_file[,-1]
  pheno_geno <- left_join(pheno_file, geno_converted)
  row.names(pheno_geno) <- pheno_geno[,1]
  pheno_geno[is.na(pheno_geno)] <- 0
  
  # run 'calculateDelta_SKAT0','calculateTukeyFences' and 'calculateMAD' function
  df_lst <- list()
  library(RIFT)
  y <- as.numeric(pheno_geno$pheno_White)
  cov <- as.matrix(pheno_geno$PC1)
  geno <- as.matrix(pheno_geno[, 13:dim(pheno_geno)[2]])
  if (dim(table(y)) !=1 && dim(geno)[2] >=5){
    class(geno) <- "numeric"
    results_white <- calculateDelta_SKAT0(y, geno.mat = geno, cov.mat = cov)
    results_white <- calculateTukeyFences(results_white)
    results_white <- calculateMAD(results_white)
    white <- results_white$ind.stats
    white$race <- 'white'
    df_lst[[1]] <- white
  }else{
    df_lst[[1]] <- NULL
    print('skip white')
  }
  
  y <- as.numeric(pheno_geno$pheno_Black)
  cov <- as.matrix(pheno_geno$PC1)
  geno <- as.matrix(pheno_geno[, 13:dim(pheno_geno)[2]])
  if (dim(table(y)) != 1 && dim(geno)[2] >=5){
    class(geno) <- "numeric"
    results_black <- calculateDelta_SKAT0(y, geno.mat = geno, cov.mat = cov)
    results_black <- calculateTukeyFences(results_black)
    results_black <- calculateMAD(results_black)
    black <- results_black$ind.stats
    black$race <- 'black'
    df_lst[[2]] <- black
  }else{
    df_lst[[2]] <- NULL
    print('skip black')
  }
  
  y <- as.numeric(pheno_geno$pheno_Asian)
  cov <- as.matrix(pheno_geno$PC1)
  geno <- as.matrix(pheno_geno[, 13:dim(pheno_geno)[2]])
  if (dim(table(y)) !=1 && dim(geno)[2] >=5){
    class(geno) <- "numeric"
    results_asian <- calculateDelta_SKAT0(y, geno.mat = geno, cov.mat = cov)
    results_asian <- calculateTukeyFences(results_asian)
    results_asian <- calculateMAD(results_asian) 
    asian <- results_asian$ind.stats
    asian$race <- 'asian'
    df_lst[[3]] <- asian
  }else{
    df_lst[[3]] <- NULL
    print('skip asian')
  }
  # merge 
  df <-do.call('bind_rows', df_lst)
  df$cancer <- cancer_id 
  write.csv(df, paste0("snp-based/04_groupwise_tests/all_skat0_plus_rift_",cancer_id, '.csv'), row.names = F)
}

# find most sig variants; merge with annotation file; count different annotation levels 
rm(list = ls())
setwd('public_data_2/snp-based/04_groupwise_tests/')
temp <-  list.files(pattern="all_skat0_plus_rift_*")
temp_name <- gsub(".csv", "", temp)
temp_name <- gsub("all_skat0_plus_rift_", "", temp_name)
for (i in 1:length(temp)) assign(temp_name[i], read.csv(temp[i],header = TRUE, stringsAsFactors = FALSE))
df_list <- list()
for (i in 1: length(temp_name)){
  df <- get(temp_name[i])
  if(dim(df)[1] == 0){next}
  df_list[[i]] <- df
}

df_list <- df_list[!unlist(lapply(df_list, is.null))]
df <- Reduce(rbind, df_list)
dim(df)
sig_df <- df[df$tukey.extreme == TRUE | df$mad.outlier == TRUE,]
write.csv(sig_df, 'all_skato_rift_outlier.csv', row.names = F)
sig_df$count <- 1
sig_df_count <- as.data.frame.matrix(table(sig_df$cancer,sig_df$race))
write.csv(sig_df_count, 'groupwise_test_summary.csv')

# data visualization
# rearrange input data for heatmap
library(tidyr)
library(ComplexHeatmap)
snp2gene <- read.csv('../tcga_pan-cancer_snp_annot_clinic.csv',header = T, stringsAsFactors = F)
snp2gene$SNP_info <- paste0(snp2gene$Chromosome,':', snp2gene$Start_Position)
snp2gene$SNP_detail <- paste0(snp2gene$Chromosome,':', snp2gene$Start_Position,'-',snp2gene$Reference_Allele, '/', snp2gene$Tumor_Seq_Allele2)
names(snp2gene)
snp2gene <- unique(snp2gene[,c(1:2,17,23:25)])
sig_gene <- merge(snp2gene, sig_df, by.x = 'SNP_detail',by.y = 'SNP_excluded',all.y = T)
write.csv(sig_gene, 'snp2gene.csv', row.names = F)

sig_gene <- sig_gene[,c(3,17:19)]
names(sig_gene)[2:3] <- c('race','cancer')
sig_wide <- as.data.frame(pivot_wider(sig_gene, names_from = Hugo_Symbol, values_from = count, values_fn = sum))
head(sig_wide)
sig_wide[is.na(sig_wide)] <- 0
sig_wide <- sig_wide[, colSums(sig_wide != 0) > 0]
sig_wide$race[sig_wide$race == 'white'] <- 'White'
sig_wide$race[sig_wide$race == 'black'] <- 'AA/B'
sig_wide$race[sig_wide$race == 'asian'] <- 'Asian'

cancer_ID = c(sort(temp_name))
race=unique(sig_wide$race)
head(sig_wide)
data = sig_wide
race = 'asian'
race_spec_forHeatmape <- function(data, race){
  sig_set_wide = data 
  sig_set_wide_spec <- sig_set_wide[sig_set_wide$race == race,]
  library(dplyr)
  cancer = cancer_ID
  sig_set_wide_spec <- left_join(data.frame(cancer), sig_set_wide_spec, by ='cancer')
  sig_set_wide_spec <- sig_set_wide_spec[,-2:-12]
  sig_set_wide_spec[is.na(sig_set_wide_spec)] <- 0
  dim(unique(sig_set_wide_spec))
  return(sig_set_wide_spec)
} # to specific cancer id in advance
black <- race_spec_forHeatmape(sig_wide,'Black')
black <- black[, colSums(black != 0) > 0]
row.names(black) <- black[,1]
black <- black[,-1]
white <- race_spec_forHeatmape(sig_wide,'White')
white <- white[, colSums(white != 0) > 0]
row.names(white) <- white[,1]
white <- white[,-1]
asian <- race_spec_forHeatmape(sig_wide,'Asian')
asian <- asian[, colSums(asian != 0) > 0]
row.names(asian) <- asian[,1]
asian <- asian[,-1]

# draw heatmap
mid_value = 10
max_value = 20
font_size_col = 8
font_size_row = 8


drawHeatmap2 <- function(race_data,race_1,race_2, mid_value, max_value){
  library(circlize)
  #race_data = white
  #race = 'white'
  col_fun = colorRamp2(c(0,mid_value,max_value), c("white","grey", "red"))
  ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points")
  race_pic <- Heatmap(as.matrix(race_data), name = "sig_variant_counts",  col = col_fun, 
                      rect_gp = gpar(col = "white", lwd = 1), 
                      cluster_columns=FALSE, cluster_rows=FALSE,
                      row_names_gp = grid::gpar(fontsize = font_size_row),
                      column_names_gp = grid::gpar(fontsize = font_size_col), 
                      column_title = race_1,
                      column_title_gp = gpar(border = "black"),
                      show_row_names = T,
                      heatmap_legend_param = list(direction = "horizontal",title_position = "lefttop"))
  
  png(paste0(race_2,'_skat_rift_sig.png'),width = 12000, height = 3000, res = 700)
  draw(race_pic, merge_legend = TRUE, heatmap_legend_side = "bottom")
  dev.off()
}
# all
drawHeatmap2(white,'White', 'white',mid_value, max_value)
drawHeatmap2(black,'African American/Black ', 'black',mid_value, max_value)
drawHeatmap2(asian,'Asian', 'asian' ,mid_value, max_value)
# subset
white_subset <- white[,colSums(white)>=10]
write.csv(as.data.frame(colnames(white_subset)), 'white_over10_genelist.csv', row.names = F)
black_subset <- black[,colSums(black)>=10]
write.csv(as.data.frame(colnames(black_subset)), 'black_over10_genelist.csv', row.names = F)
asian_subset <- asian[,colSums(asian)>=10]
write.csv(as.data.frame(colnames(asian_subset)), 'asian_over10_genelist.csv', row.names = F)

drawHeatmap1 <- function(white, black, asian, mid_value, max_value){
  library(circlize)
  col_fun = colorRamp2(c(0,mid_value,max_value), c("white","grey", "red"))
  ht_opt$TITLE_PADDING = unit(c(4.5, 4.5), "points")
  white_pic <- Heatmap(t(as.matrix(white)), name = "mat1",  col = col_fun, 
                       rect_gp = gpar(col = "white", lwd = 1), 
                       cluster_columns=FALSE, cluster_rows=FALSE,
                       column_names_gp = grid::gpar(fontsize = font_size_col),
                       row_names_gp = grid::gpar(fontsize = font_size_row), 
                       row_title = 'White',
                       row_title_gp = gpar(border = "black"),
                       show_column_names = F,
                       show_heatmap_legend = F)
  black_pic <- Heatmap(t(as.matrix(black)), name = "mat2",  col = col_fun, 
                       rect_gp = gpar(col = "white", lwd = 1), 
                       cluster_columns=FALSE, cluster_rows=FALSE,
                       column_names_gp = grid::gpar(fontsize = font_size_col),
                       row_names_gp = grid::gpar(fontsize = font_size_row), 
                       row_title = 'AA/B',
                       row_title_gp = gpar(border = "black"),
                       show_column_names = F,
                       show_heatmap_legend = F)
  asian_pic <- Heatmap(t(as.matrix(asian)), name = "sig_variant_counts",  col = col_fun, 
                       rect_gp = gpar(col = "white", lwd = 1), 
                       cluster_columns=FALSE, cluster_rows=FALSE,
                       column_names_gp = grid::gpar(fontsize = font_size_col),
                       row_names_gp = grid::gpar(fontsize = font_size_row), 
                       row_title = 'Asian',
                       row_title_gp = gpar(border = "black"),
                       show_column_names = T,
                       heatmap_legend_param = list(direction = "horizontal",title_position = "lefttop"))
  
  ht_list = white_pic %v% black_pic %v% asian_pic 
  png(paste0('all_race_skat_rift_sig.png'),width = 5600, height = 5000, res = 700)
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom")
  dev.off()
}

drawHeatmap1(white_subset, black_subset, asian_subset, mid_value, max_value)

