library(dplyr)
library(tidyr)
library(DESeq2)
library(pheatmap)
library(viridis)
library(gridExtra)


#### LOAD DATASETS  ####


#Get pseudobulk counts, delete columns related to HSC and rename all columns to "pseudo_"sample name""

pseudobulk_counts <- read.delim("raw_pseudo_counts")

pseudobulk_counts <- pseudobulk_counts %>%
  dplyr::select(-hsc_d, -hsc_t) %>%
  rename_with(~ c("gene", "pseudo_DD", "pseudo_DL", "pseudo_DLT", "pseudo_DT", "pseudo_TD", 
                  "pseudo_TL", "pseudo_TLT", "pseudo_TT"))



#Get bulk counts and rename columns to "bulk_"sample name""

bulk_counts <- read.delim("raw_bulk_counts")

colnames(bulk_counts) <- gsub("^iMo_", "bulk_", colnames(bulk_counts))


#Get bulk TPMs, rename columns to bulk_"sample name"" and only keep TPMs > 1

bulk_tpm <- read.delim("bulk_tpms")

bulk_tpm <- bulk_tpm[-c(1, 2)]

names(bulk_tpm)[-1] <- paste0("bulk_", names(bulk_tpm)[-1])

bulk_tpm <- bulk_tpm%>%
  filter(apply(., 1, function(x) all(x > 1)))



####  ADD SUM COLUMNS ####


#Add sum columns in the pseudobulk counts dataframe

pseudobulk_counts <- pseudobulk_counts %>%
  mutate(
    sum_pseudo_noLPS = rowSums(dplyr::select(., pseudo_DD, pseudo_DT, pseudo_TD, pseudo_TT)),
    sum_pseudo_LPS = rowSums(dplyr::select(., pseudo_DL, pseudo_DLT, pseudo_TL, pseudo_TLT))
  )


##Add sum columns in the bulk counts dataframe

bulk_counts <- bulk_counts %>%
  mutate(
    sum_bulk_noLPS = rowSums(dplyr::select(., bulk_DD, bulk_DT, bulk_TD, bulk_TT)),
    sum_bulk_LPS = rowSums(dplyr::select(., bulk_DL, bulk_DLT, bulk_TL, bulk_TLT))
  )


#Add sum columns in the bulk tpm dataframe

bulk_tpm <- bulk_tpm %>%
  mutate(
    sum_bulk_noLPS = rowSums(dplyr::select(., bulk_DD, bulk_DT, bulk_TD, bulk_TT)),
    sum_bulk_LPS = rowSums(dplyr::select(., bulk_DL, bulk_DLT, bulk_TL, bulk_TLT))
  )



####  COMPARE RAW COUNTS VS RAW COUNTS  ####


#Merge counts dataframes for bulk and pseudobulk

merged_counts <- inner_join(pseudobulk_counts, bulk_counts, by = "gene")


# Convert all columns except for "gene" to numeric

merged_counts <- merged_counts %>%
  dplyr::mutate(across(-gene, as.numeric))


#Keep only numeric objects and convert genes to rownames

merged_counts_num <- merged_counts %>% column_to_rownames(var = "gene")


#Split again into bulk / pseudobulk

pseudobulk_counts_num <- merged_counts_num %>%
  dplyr::select(1:10)

bulk_counts_num <- merged_counts_num %>%
  dplyr::select(11:20)


# Calculate the Pearson correlation
correlation_raw_counts <- cor(pseudobulk_counts_num, bulk_counts_num, method = "pearson")

correlation_raw_counts <- as.data.frame(correlation_raw_counts)



####  COMPARE NORMALISED COUNTS VS NORMALISED COUNTS


#Create metadata dataframe for the merged counts dataframe

metaData_merged_counts <- data.frame(
  Treatment = c("DD", "DL", "DLT", "DT", "TD", "TL", "TLT", "TT", "noLPS", "LPS",
                "TLT", "TL", "DLT", "DL", "TT", "TD", "DT", "DD", "noLPS", "LPS"),
  Samples = c("pseudo_DD", "pseudo_DL", "pseudo_DLT", "pseudo_DT", "pseudo_TD", 
              "pseudo_TL", "pseudo_TLT", "pseudo_TT", "sum_pseudo_noLPS", "sum_pseudo_LPS",
              "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL", "bulk_TT", "bulk_TD", "bulk_DT",
              "bulk_DD", "sum_bulk_noLPS", "sum_bulk_LPS"))


#Get the normalised counts with DESeq2

dds_merged_counts <- DESeqDataSetFromMatrix(countData = merged_counts_num, colData = metaData_merged_counts, design = ~ Treatment)
dds_merged_counts <- estimateSizeFactors(dds_merged_counts)
merged_norm_counts <- as.data.frame(counts(dds_merged_counts, normalized = TRUE))


#Divide again into pseudobulk / bulk

pseudobulk_norm_counts <- merged_norm_counts %>%
  dplyr::select(1:10)

bulk_norm_counts <- merged_norm_counts %>%
  dplyr::select(11:20)

# Calculate Pearson correlation for normalised counts

correlation_norm_counts <- as.data.frame(cor(pseudobulk_norm_counts, bulk_norm_counts, method = "pearson"))



####  COMPARE CPMs VS CPMs  ####

# Get CPMs
merged_cpm <- sweep(merged_norm_counts, 2, colSums(merged_norm_counts), "/") * 1e6

# Filter out CPMs lower than 2 MADs from the median 

global_median <- median(unlist(merged_cpm), na.rm = TRUE)
global_mad <- mad(unlist(merged_cpm), na.rm = TRUE)
global_threshold <- global_median - 2 * global_mad

# Then filter rows based on threshold
merged_cpm_filtered <- merged_cpm[rowSums(merged_cpm >= global_threshold) == ncol(merged_cpm), ]


#Separate again into pseudobulk / bulk

pseudobulk_cpm <- merged_cpm_filtered[, 1:10]
bulk_cpm <- merged_cpm_filtered[, 11:20]


# Calculate Pearson correlation for CPMs

correlation_cpm <- cor(pseudobulk_cpm, bulk_cpm, method = "pearson")



####  GET NORMALISED COUNTS AND CPMS (ONLY PSEUDOBULK)  ####


#Create metadata

metaData_pseudobulk <- data.frame(
  Treatment = c("DD", "DL", "DLT", "DT", "TD", "TL", "TLT", "TT", "noLPS", "LPS"),
  row.names = c("pseudo_DD", "pseudo_DL", "pseudo_DLT", "pseudo_DT", "pseudo_TD", 
                "pseudo_TL", "pseudo_TLT", "pseudo_TT", "sum_pseudo_noLPS", "sum_pseudo_LPS")) 


#Calculate normalised  with DESeq

dds_pseudo <- DESeqDataSetFromMatrix(countData = pseudobulk_counts_num, colData = metaData_pseudobulk, design = ~ Treatment)
dds_pseudo <- estimateSizeFactors(dds_pseudo)
pseudobulk_counts_only <- counts(dds_pseudo, normalized = TRUE) 
pseudobulk_cpm_only <- sweep(pseudobulk_counts_only, 2, colSums(pseudobulk_counts_only), "/") * 1e6


#Filter high reads for CPM

median <- apply(pseudobulk_cpm_only, 2, median, na.rm = TRUE)  
mads <- apply(pseudobulk_cpm_only, 2, mad, na.rm = TRUE)        


#Calculate the lower threshold per sample

threshold_cpm_only <- median - 2 * mads  


#Filter genes based on threshold for CPMs

pseudobulk_cpm_only_filtered <- pseudobulk_cpm_only[apply(pseudobulk_cpm_only, 1, function(row) {
  all(row >= threshold_cpm_only)  # Keep genes that meet the threshold in all samples
}), ]



####  COMPARE CPMs AND COUNTS  ####


# Convert pseudobulk normalized counts to CPMs

pseudobulk_cpm_vs_counts <- sweep(pseudobulk_norm_counts, 2, colSums(pseudobulk_norm_counts), "/") * 1e6

median <- apply(pseudobulk_cpm_vs_counts, 2, median, na.rm = TRUE)  
mads <- apply(pseudobulk_cpm_vs_counts, 2, mad, na.rm = TRUE)        


#Calculate the lower threshold per sample

threshold_cpm_vs_counts <- median - 2 * mads  


#Filter genes based on threshold for CPMs

pseudobulk_cpm_vs_counts_filtered <- pseudobulk_cpm_vs_counts[apply(pseudobulk_cpm_vs_counts, 1, function(row) {
  all(row >= threshold_cpm_vs_counts)  # Keep genes that meet the threshold in all samples
}), ]


correlation_count_cpm <- cor(pseudobulk_cpm_vs_counts_filtered, bulk_norm_counts, method = "pearson")



####  COMPARE COUNTS (pseudobulk) VS TPMs  ####


#Merge by the "gene" column

merged_counts_tpm <- inner_join(pseudobulk_counts, bulk_tpm, by = "gene")


#Convert the "gene" column back to rownames

merged_counts_tpm <- merged_counts_tpm %>% column_to_rownames(var = "gene")


#Separate pseudobulk / bulk

pseudobulk_counts_tpm <- merged_counts_tpm[, 1:10]
bulk_counts_tpm <- merged_counts_tpm[, 11:20]


# Calculate Pearson correlation for CPMs

correlation_count_tpm <- cor(pseudobulk_counts_tpm, bulk_counts_tpm, method = "pearson")



####  COMPARE CPM VS TPM  ####


#COMPARE COUNTS (pseudobulk) VS TPMs


#Merge by the "gene" column

pseudobulk_cpm_only_temp <- as.data.frame(pseudobulk_cpm_only_filtered) %>% rownames_to_column(var = "gene")

merged_cpm_tpm <- inner_join(pseudobulk_cpm_only_temp, bulk_tpm, by = "gene")


#Convert the "gene" column back to rownames

merged_cpm_tpm <- merged_cpm_tpm %>% column_to_rownames(var = "gene")


#Separate pseudobulk / bulk

pseudobulk_cpm_tpm <- merged_cpm_tpm[, 1:10]
bulk_cpm_tpm <- merged_cpm_tpm[, 11:20]


# Calculate Pearson correlation for CPMs

correlation_cpm_tpm <- cor(pseudobulk_cpm_tpm, bulk_cpm_tpm, method = "pearson")



####  PLOT HEATMAPS  ####


#Find min and max values of all correlations 

global_min <- min(correlation_raw_counts, correlation_norm_counts, correlation_cpm, 
                  correlation_count_tpm, correlation_cpm_tpm)

global_max <- max(correlation_raw_counts, correlation_norm_counts, correlation_cpm, 
                  correlation_count_tpm, correlation_cpm_tpm)


#Define a consistent color scale

color_scale <- viridis(50, option = "magma")  # Use any color palette you prefer

color_scale <- colorRampPalette(c("red", "#FFFDD0", "blue"))(50)


#Plot raw counts 

heatmap_counts <-pheatmap(correlation_raw_counts, 
                          color = color_scale, 
                          breaks = seq(global_min, global_max, length.out = 50),  
                          display_numbers = TRUE, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          main = "Correlation Matrix - raw counts")


#Plot normalised counts 

pheatmap(correlation_norm_counts, 
         color = color_scale, 
         breaks = seq(global_min, global_max, length.out = 50),  
         display_numbers = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Correlation Matrix - normalised counts")


#Plot CPMs 

pheatmap(correlation_cpm, 
         color = color_scale, 
         breaks = seq(global_min, global_max, length.out = 50),  
         display_numbers = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Correlation Matrix - CPM vs CPM")


#Plot TPM vs counts 

heatmap_count_tpm <-pheatmap(correlation_count_tpm, 
                             color = color_scale, 
                             breaks = seq(global_min, global_max, length.out = 50),  
                             display_numbers = TRUE, 
                             clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "euclidean",
                             main = "Correlation Matrix - counts vs TPM")


#Plot TPM vs CPM 

heatmap_tpm_cpm <-pheatmap(correlation_cpm_tpm, 
                           color = color_scale, 
                           breaks = seq(global_min, global_max, length.out = 50),  
                           display_numbers = TRUE, 
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           main = "Correlation Matrix - CPM vs TPM")


#Plot counts vs CPM 

heatmap_count_cpm <- pheatmap(correlation_count_cpm, 
                              color = color_scale, 
                              breaks = seq(global_min, global_max, length.out = 50),  
                              display_numbers = TRUE, 
                              clustering_distance_rows = "euclidean",
                              clustering_distance_cols = "euclidean",
                              main = "Correlation Matrix - CPM (pseudobulk) vs counts")



grid.arrange(heatmap_count_cpm$gtable, heatmap_count_tpm$gtable,
             heatmap_counts$gtable, heatmap_tpm_cpm$gtable,
             ncol = 2, nrow = 2)



####  COMPARE SUM VALUES FOR LPS / noLPS IN BOTH BULK AND PSEUDOBULK ####


correlation_LPS_noLPS_pseudo <- cor(pseudobulk_counts$sum_pseudo_LPS, pseudobulk_counts$sum_pseudo_noLPS, 
                                    method = "pearson")

correlation_LPS_noLPS_bulk <- cor(bulk_tpm$sum_bulk_LPS, bulk_tpm$sum_bulk_noLPS, 
                                  method = "pearson")




#Calculate correlation between LPS / noLPS in the LPS-inducible dataset

merged_counts_LPS_genes <- inner_join(pseudobulk_counts, LPS_inducible_genes, by = c("gene" = "MoMf_LPS_inducible_genes"))

correlation_LPS_genes <- cor(merged_counts_LPS_genes$sum_pseudo_LPS, merged_counts_LPS_genes$sum_pseudo_noLPS, 
                             method = "pearson")


#Same but with the mean instead of the sum, which is the same because of how pearson is calculated

merged_counts_LPS_genes <- merged_counts_LPS_genes %>%
  mutate(
    mean_pseudo_noLPS = rowSums(.[, c("pseudo_DD", "pseudo_DT", "pseudo_TD", "pseudo_TT")], na.rm = TRUE),
    mean_pseudo_LPS = rowSums(.[, c("pseudo_DL", "pseudo_DLT", "pseudo_TL", "pseudo_TLT")], na.rm = TRUE))

correlation_LPS_genes_mean <- cor(merged_counts_LPS_genes$mean_pseudo_LPS, merged_counts_LPS_genes$mean_pseudo_noLPS, 
                                  method = "pearson")

####  OBTAIN THE CORRELATION MATRICES WITHOUT THE SUM COLUMNS  ####


# Remove sum columns

pseudobulk_counts_only <- merged_counts_num %>%
  dplyr::select(starts_with("pseudo_")) 

bulk_counts_only <- merged_counts_num %>%
  dplyr::select(starts_with("bulk_")) 


# Matrix 1: Counts vs Counts
correlation_raw_counts_no_sum <- cor(pseudobulk_counts_only, bulk_counts_only, method = "pearson")


# Matrix 2: Pseudobulk CPM vs Bulk counts

pseudobulk_cpm_vs_counts_only <- pseudobulk_cpm_vs_counts_filtered %>%
  dplyr::select(-sum_pseudo_noLPS, -sum_pseudo_LPS)

bulk_norm_counts_only <- bulk_norm_counts %>%
  dplyr::select(-sum_bulk_noLPS, -sum_bulk_LPS)

correlation_count_cpm_no_sum <- cor(pseudobulk_cpm_vs_counts_only, bulk_norm_counts_only, method = "pearson")


# Matrix 3: Pseudobulk counts vs Bulk TPM

# Use merged_counts_tpm, which includes sums
pseudobulk_counts_tpm_only <- merged_counts_tpm %>%
  dplyr::select(starts_with("pseudo_"))

bulk_counts_tpm_only <- merged_counts_tpm %>%
  dplyr::select(starts_with("bulk_")) 

correlation_count_tpm_no_sum <- cor(pseudobulk_counts_tpm_only, bulk_counts_tpm_only, method = "pearson")


#️Matrix 4: Pseudobulk CPM vs Bulk TPM

pseudobulk_cpm_tpm_only <- merged_cpm_tpm %>%
  dplyr::select(starts_with("pseudo_"))

bulk_cpm_tpm_only <- merged_cpm_tpm %>%
  dplyr::select(starts_with("bulk_")) 

correlation_cpm_tpm_no_sum <- cor(pseudobulk_cpm_tpm_only, bulk_cpm_tpm_only, method = "pearson")


#Plot all four matrices

# Get new global min/max
global_min_no_sum <- min(correlation_raw_counts_no_sum, 
                         correlation_count_cpm_no_sum, 
                         correlation_count_tpm_no_sum, 
                         correlation_cpm_tpm_no_sum)

global_max_no_sum <- max(correlation_raw_counts_no_sum, 
                         correlation_count_cpm_no_sum, 
                         correlation_count_tpm_no_sum, 
                         correlation_cpm_tpm_no_sum)

color_scale <- colorRampPalette(c("red", "#FFFDD0", "blue"))(50)

h1 <- pheatmap(correlation_raw_counts_no_sum,
               color = color_scale,
               breaks = seq(global_min_no_sum, global_max_no_sum, length.out = 50),
               display_numbers = TRUE,
               main = "Counts vs Counts")

h2 <- pheatmap(correlation_count_cpm_no_sum,
               color = color_scale,
               breaks = seq(global_min_no_sum, global_max_no_sum, length.out = 50),
               display_numbers = TRUE,
               main = "Pseudobulk CPM vs Bulk counts", 
               legend = FALSE)

h3 <- pheatmap(correlation_count_tpm_no_sum,
               color = color_scale,
               breaks = seq(global_min_no_sum, global_max_no_sum, length.out = 50),
               display_numbers = TRUE,
               main = "Pseudobulk counts vs Bulk TPM",
               legend = FALSE)

h4 <- pheatmap(correlation_cpm_tpm_no_sum,
               color = color_scale,
               breaks = seq(global_min_no_sum, global_max_no_sum, length.out = 50),
               display_numbers = TRUE,
               main = "Pseudobulk CPM vs Bulk TPM",
               legend = FALSE)

# Arrange in grid
grid.arrange(h1$gtable, h2$gtable, h3$gtable, h4$gtable, ncol = 2)

