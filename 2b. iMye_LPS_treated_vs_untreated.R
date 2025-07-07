library(dplyr)
library(tibble)
library(ggplot2)
library(DESeq2)
library(BioVenn)
library(edgeR)
library(ggvenn)


####  GET DATASETS  ####


##  PSEUDOBULK  ##


#Get pseudo dataset (raw counts of the 8 iMye - and 2 iHSC, which will be eliminated - in pseudobulk format)

pseudobulk_for_diff <- read.delim("counts_pseudo")


#Remove HSC (day12) data

pseudobulk_for_diff <- pseudobulk_for_diff %>%
  dplyr::dplyr::dplyr::select(-hsc_d, -hsc_t)


#Rename the columns from "mono_condition" to "pseudo_condition"

colnames(pseudobulk_for_diff) <- gsub("^mono_", "pseudo_", colnames(pseudobulk_for_diff))


##  BULK  ##


#Get dataset (raw counts of the 8 iMye samples in bulk format)

bulk_for_diff <- read.delim("counts_bulk")

colnames(bulk_for_diff) <- gsub("^iMo_", "bulk_", colnames(bulk_for_diff))


##  GENE LENGTHS  ##


gene_lengths <- read.delim("gene_lengths_from_Seq2Science")


## GET KZnF INVENTORY  ##


krab_krueppel_genes <- read.csv("KZnF_Inventory")


##  LOAD LIST OF KNOWN LPS-INDUCIBLE GENES  ##

LPS_inducible_genes <- read.delim("known_LPS_inducible_genes")



####  GET GENE LENGTH DATASET  ####
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Your genes of interest
genes_of_interest <- c(
  "A1BG", "A1BG-AS1", "A1CF", "A2M", "A2M-AS1", "A2ML1", "A2ML1-AS1", "A2MP1",
  "A3GALT2", "A4GALT", "A4GNT", "AA06", "AAAS", "AACS", "AACSP1", "AADAC",
  "AADACL2", "AADACL2-AS1", "AADACL3", "AADACL4"
)

# 1. Get genomic lengths (gene start to end)
genomic_lengths <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes_of_interest,
  mart = ensembl
) %>%
  mutate(genomic_length = end_position - start_position + 1)

# 2. Get mRNA transcript lengths (sum of exons per transcript)
transcript_lengths <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_length"),
  filters = "hgnc_symbol",
  values = genes_of_interest,
  mart = ensembl
)

# 3. Choose max mRNA length per gene (most complete transcript)
mRNA_lengths <- transcript_lengths %>%
  group_by(hgnc_symbol) %>%
  summarise(mRNA_length = max(transcript_length, na.rm = TRUE), .groups = "drop")

# 4. Combine both lengths
gene_length_summary <- genomic_lengths %>%
  distinct(hgnc_symbol, genomic_length) %>%
  inner_join(mRNA_lengths, by = "hgnc_symbol")

# 5. Optional: add your values for comparison
your_lengths <- data.frame(
  hgnc_symbol = c("A1BG", "A1BG-AS1", "A1CF", "A2M", "A2M-AS1", "A2ML1", "A2ML1-AS1", "A2MP1", 
                  "A3GALT2", "A4GALT", "A4GNT", "AA06", "AAAS", "AACS", "AACSP1", "AADAC", 
                  "AADACL2", "AADACL2-AS1", "AADACL3", "AADACL4"),
  S2S_length = c(3374, 2126, 9466, 4918, 2298, 5445, 817, 1192, 1018, 3037, 1886, 629,
                 1799, 6078, 2814, 1558, 5055, 768, 4051, 2023)
)

# Final comparison table
final_df <- gene_length_summary %>%
  left_join(your_lengths, by = "hgnc_symbol") %>%
  dplyr::select(hgnc_symbol, S2S_length, mRNA_length, genomic_length)

# View result
print(final_df)





gene_list <- unique(gene_lengths$gene_id)

# Query gene start and end positions to calculate gene length
gene_info <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "start_position", "end_position"),
  filters = ifelse(grepl("^ENS", gene_list[1]), "ensembl_gene_id", "hgnc_symbol"),
  values = gene_list,
  mart = ensembl
)

# Calculate gene length = end - start + 1
gene_info <- gene_info %>%
  mutate(gene_length = end_position - start_position + 1)

# Now merge back to your gene_lengths dataframe by gene symbol or Ensembl ID
# Use the appropriate join based on your gene ID type
if(grepl("^ENS", gene_list[1])){
  merged_lengths <- gene_lengths %>%
    left_join(gene_info, by = c("gene_id" = "ensembl_gene_id"))
} else {
  merged_lengths <- gene_lengths %>%
    left_join(gene_info, by = c("gene_id" = "hgnc_symbol"))
}

# View results
head(merged_lengths)
####  PSEUDOBULK-ONLY DAY19 LPS/noLPS DIFF. EXPR.  ####


##  PREPARE FOR DESeq2  ##


#Set gene column as rownames

pseudobulk_for_diff <- pseudobulk_for_diff %>% column_to_rownames(var = "gene")


#Create metadata for the dataset, identifying LPS/noLPS groups

metaData_pseudo_for_diff <- data.frame(
  Treatment = c("noLPS", "LPS", "LPS", "noLPS", "noLPS", "LPS", "LPS", "noLPS" ),
  Samples = c("pseudo_dd", "pseudo_dl", "pseudo_dlt", "pseudo_dt", "pseudo_td","pseudo_tl", "pseudo_tlt", "pseudo_tt"))


#Set sample names as rownames

metaData_pseudo_for_diff <- metaData_pseudo_for_diff %>% column_to_rownames(var = "Samples")


#Make Treatment levels into factors

metaData_pseudo_for_diff$Treatment <- factor(metaData_pseudo_for_diff$Treatment)




##  RUN DESeq2 AND GET THE LPS-BIASED GENES  ##


dds_pseudo <- DESeqDataSetFromMatrix(countData = pseudobulk_for_diff, 
                                     colData = metaData_pseudo_for_diff, 
                                     design = ~ Treatment)

dds_pseudo <- DESeq(dds_pseudo)

results_pseudo <- results(dds_pseudo, contrast = c("Treatment", "LPS", "noLPS"))

results_pseudo <- as.data.frame(results_pseudo)


#Filter genes with log2FC > 1, padj < 0.01, and non-NA padj

LPS_pseudo_only <- results_pseudo %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj))

LPS_pseudo_only <- rownames_to_column(LPS_pseudo_only, var = "gene")



####  BULK-ONLY DAY19 LPS/noLPS DIFF. EXPR.  ####

##  PREPARE FOR DESeq2  ##


#Set gene column as rownames

bulk_for_diff <- bulk_for_diff %>% column_to_rownames(var = "gene")


#Create metadata for the dataset, identifying LPS/noLPS groups

metaData_bulk_for_diff <- data.frame(
  Treatment = c("LPS", "LPS", "LPS", "LPS", "noLPS", "noLPS", "noLPS", "noLPS" ),
  Samples = c("bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL", "bulk_TT", "bulk_TD", "bulk_DT", "bulk_DD"))


#Set sample names as rownames

metaData_bulk_for_diff <- metaData_bulk_for_diff %>% column_to_rownames(var = "Samples")


#Make Treatment levels into factors

metaData_bulk_for_diff$Treatment <- factor(metaData_bulk_for_diff$Treatment)



##  RUN DESeq2 AND GET THE LPS-BIASED GENES  ##


dds_bulk <- DESeqDataSetFromMatrix(countData = bulk_for_diff, 
                                   colData = metaData_bulk_for_diff, 
                                   design = ~ Treatment)

dds_bulk <- DESeq(dds_bulk)

results_bulk <- results(dds_bulk, contrast = c("Treatment", "LPS", "noLPS"))

results_bulk <- as.data.frame(results_bulk)


#Filter genes with log2FC > 1, padj < 0.01, and non-NA padj

LPS_bulk_only <- results_bulk %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj))

LPS_bulk_only <- rownames_to_column(LPS_bulk_only, var = "gene")




####  BULK-PSEUDO DAY19 LPS/noLPS DIFF. EXPR. STARTING FROM COUNTS (pseudo) OR COUNTS/mRNA LENGTH (bULK)  ####


###  MERGE AND PREPARE DATASETS  ###


# Merge datasets and normalize bulk counts by gene length

bulk_for_diff <- rownames_to_column(bulk_for_diff, var = "gene")
pseudobulk_for_diff <- rownames_to_column(pseudobulk_for_diff, var = "gene")

merged_counts_pseudo_and_counts_norm_bulk <- inner_join(pseudobulk_for_diff, bulk_for_diff, by = "gene") %>% 
  inner_join(merged_lengths, by = join_by("gene" == "gene_id")) %>% 
  mutate(across(10:17, ~ .x / .data[["gene_lengths"]])) %>%
  dplyr::select(-gene_lengths) %>% 
  distinct(gene, .keep_all = TRUE) %>%   # keeps only first occurrence of each gene
  column_to_rownames(var = "gene")


#Get dataset with pseudobulk LPS samples and bulk noLPS samples (normalised by gene length)

pseudo_LPS_bulk_noLPS_norm <- merged_counts_pseudo_and_counts_norm_bulk %>% 
  dplyr::select(-c("pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
                   "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL"))

#Get dataset with pseudobulk LPS samples and bulk noLPS samples (normalised by gene length)

pseudo_noLPS_bulk_LPS_norm <- merged_counts_pseudo_and_counts_norm_bulk %>% 
  dplyr::dplyr::select(c("pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
                         "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL"))




### PSEUDO LPS VS BULK noLPS  ###


#Create metadata and data matrix

meta1 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("pseudo_dl", "pseudo_dlt", "pseudo_tl", "pseudo_tlt",
              "bulk_TT", "bulk_TD", "bulk_DT", "bulk_DD")
) %>% column_to_rownames("Samples")

meta1$Treatment <- factor(meta1$Treatment, levels = c("noLPS", "LPS"))

matrix1 <- as.matrix(pseudo_LPS_bulk_noLPS_norm)


#Run normalisation

dge1 <- DGEList(counts = matrix1)

dge1 <- calcNormFactors(dge1)

design1 <- model.matrix(~ Treatment, data = meta1)

v1 <- voom(dge1, design1)

fit1 <- lmFit(v1, design1)

fit1 <- eBayes(fit1)


#Get results and rename it like DESeq2

results_pseudo_LPS_bulk_noLPS_norm <- topTable(fit1, coef = "TreatmentLPS", number = Inf, sort.by = "p")

colnames(results_pseudo_LPS_bulk_noLPS_norm)[colnames(results_pseudo_LPS_bulk_noLPS_norm) == "logFC"] <- "log2FoldChange"
colnames(results_pseudo_LPS_bulk_noLPS_norm)[colnames(results_pseudo_LPS_bulk_noLPS_norm) == "P.Value"] <- "pvalue"
colnames(results_pseudo_LPS_bulk_noLPS_norm)[colnames(results_pseudo_LPS_bulk_noLPS_norm) == "adj.P.Val"] <- "padj"


#Get LPS-inducible genes

LPS_pLPS_bnoLPS <- results_pseudo_LPS_bulk_noLPS_norm %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj))

LPS_pLPS_bnoLPS <- rownames_to_column(LPS_pLPS_bnoLPS, var = "gene")



###  BULK LPS VS PSEUDO noLPS  ###


#Create metadata and data matrix

meta2 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
              "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL")
) %>% column_to_rownames("Samples")

meta2$Treatment <- factor(meta2$Treatment, levels = c("noLPS", "LPS"))

matrix2 <- as.matrix(pseudo_noLPS_bulk_LPS_norm)


#Run normalisation

dge2 <- DGEList(counts = matrix2)

dge2 <- calcNormFactors(dge2)

design2 <- model.matrix(~ Treatment, data = meta2)

v2 <- voom(dge2, design2)

fit2 <- lmFit(v2, design2)

fit2 <- eBayes(fit2)


#Get results and rename it like DESeq2

results_pseudo_noLPS_bulk_LPS_norm <- topTable(fit2, coef = "TreatmentLPS", number = Inf, sort.by = "p")

colnames(results_pseudo_noLPS_bulk_LPS_norm)[colnames(results_pseudo_noLPS_bulk_LPS_norm) == "logFC"] <- "log2FoldChange"
colnames(results_pseudo_noLPS_bulk_LPS_norm)[colnames(results_pseudo_noLPS_bulk_LPS_norm) == "P.Value"] <- "pvalue"
colnames(results_pseudo_noLPS_bulk_LPS_norm)[colnames(results_pseudo_noLPS_bulk_LPS_norm) == "adj.P.Val"] <- "padj"


#Get LPS-inducible genes

LPS_bLPS_pnoLPS <- results_pseudo_noLPS_bulk_LPS_norm %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj))

LPS_bLPS_pnoLPS <- rownames_to_column(LPS_bLPS_pnoLPS, var = "gene")


####  GET OVERLAP OF LPS-INDUCIBLE GENES IDENTIFIED WITH BULK, PSEUDO OR CROSS  ####


#Make venn diagram with the four groups of different methods

LPS_inducible_found_sets <- list(
  Bulk = LPS_bulk_only$gene,
  Pseudo = LPS_pseudo_only$gene,
  pLPS_bnoLPS = LPS_pLPS_bnoLPS$gene,
  pnoLPS_bLPS = LPS_bLPS_pnoLPS$gene)

ggvenn(LPS_inducible_found_sets, 
       show_percentage = FALSE,
       set_name_size = 3.5, 
       text_size = 5)


#calculate ratios of intersecting genes

#Make function for ratios

calculate_intersection_ratios <- function(genes_1, genes_2, label_1, label_2) {
  genes_1 %>%
    intersect(genes_2) %>%
    length() %>%
    { 
      intersection_count <- .
      print(paste("Genes in the intersection of", label_1, "and", label_2, ":", intersection_count))  
      genes_1_ratio <- round(intersection_count / length(genes_1) * 100, 2)
      genes_2_ratio <- round(intersection_count / length(genes_2) * 100, 2)
      print(paste("Which is a ratio of:", genes_1_ratio, "% for", label_1, "and", genes_2_ratio, "% for", label_2))
    }
}


#Ratio of intersecting genes bulk-bulk / pseudo-pseudo
calculate_intersection_ratios(LPS_bulk_only$gene, LPS_pseudo_only$gene, label_1 = "bulk/bulk", label_2 = "pseudo/pseudo")






####  BULK-PSEUDO DAY19 LPS/noLPS DIFF. EXPR. STARTING FROM COUNTS (pseudo) OR COUNTS/GENE LENGTH (bULK)  ####


###  MERGE AND PREPARE DATASETS  ###


# Add gene column
if (!"gene" %in% colnames(pseudobulk_for_diff)) {
  pseudobulk_for_diff <- rownames_to_column(pseudobulk_for_diff, var = "gene")}
if (!"gene" %in% colnames(bulk_for_diff)) {
  bulk_for_diff <- rownames_to_column(bulk_for_diff, var = "gene")}

# Merge gene lengths with bulk data
bulk_for_diff <- inner_join(bulk_for_diff, merged_lengths, by = c("gene" = "gene_id"))

# Normalize bulk by gene length (adjust column indexes to your bulk sample columns!)
bulk_for_diff_norm <- bulk_for_diff %>%
  mutate(across(6:9, ~ .x / gene_lengths)) %>%
  dplyr::dplyr::select(-gene_lengths)

# Merge pseudobulk and normalized bulk
merged_counts <- inner_join(pseudobulk_for_diff, bulk_for_diff_norm, by = "gene") %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene")


### Create matrices for differential expression ###

# Pseudobulk LPS vs Bulk noLPS
pseudo_LPS_bulk_noLPS_norm <- merged_counts %>% 
  dplyr::dplyr::select(pseudo_dl, pseudo_dlt, pseudo_tl, pseudo_tlt, bulk_TT, bulk_TD, bulk_DT, bulk_DD)

# Pseudobulk noLPS vs Bulk LPS
pseudo_noLPS_bulk_LPS_norm <- merged_counts %>% 
  dplyr::dplyr::select(pseudo_tt, pseudo_td, pseudo_dt, pseudo_dd, bulk_TLT, bulk_TL, bulk_DLT, bulk_DL)

# Clean matrices (remove negative values)
matrix1 <- as.matrix(pseudo_LPS_bulk_noLPS_norm)
matrix1 <- matrix1[!apply(matrix1 < 0, 1, any), ]

matrix2 <- as.matrix(pseudo_noLPS_bulk_LPS_norm)
matrix2 <- matrix2[!apply(matrix2 < 0, 1, any), ]


### Create metadata ###

meta1 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("pseudo_dl", "pseudo_dlt", "pseudo_tl", "pseudo_tlt",
              "bulk_TT", "bulk_TD", "bulk_DT", "bulk_DD")
) %>% column_to_rownames("Samples")
meta1$Treatment <- factor(meta1$Treatment, levels = c("noLPS", "LPS"))

meta2 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
              "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL")
) %>% column_to_rownames("Samples")
meta2$Treatment <- factor(meta2$Treatment, levels = c("noLPS", "LPS"))


### limma-voom differential expression ###

# DE: Pseudobulk LPS vs Bulk noLPS
matrix1_clean <- matrix1[complete.cases(matrix1) & !apply(matrix1 < 0, 1, any), ]
dge1 <- DGEList(counts = matrix1_clean)
dge1 <- calcNormFactors(dge1)
design1 <- model.matrix(~ Treatment, data = meta1)
v1 <- voom(dge1, design1)
fit1 <- lmFit(v1, design1)
fit1 <- eBayes(fit1)

results_pseudo_LPS_bulk_noLPS_norm <- topTable(fit1, coef = "TreatmentLPS", number = Inf, sort.by = "p") %>%
  rownames_to_column("gene") %>%
  rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

LPS_pLPS_bnoLPS <- results_pseudo_LPS_bulk_noLPS_norm %>%
  filter(log2FoldChange > 1, padj < 0.01, !is.na(padj))


# DE: Bulk LPS vs Pseudobulk noLPS
matrix2_clean <- matrix2[complete.cases(matrix2) & !apply(matrix2 < 0, 1, any), ]
dge2 <- DGEList(counts = matrix2_clean)
dge2 <- calcNormFactors(dge2)
design2 <- model.matrix(~ Treatment, data = meta2)
v2 <- voom(dge2, design2)
fit2 <- lmFit(v2, design2)
fit2 <- eBayes(fit2)

results_pseudo_noLPS_bulk_LPS_norm <- topTable(fit2, coef = "TreatmentLPS", number = Inf, sort.by = "p") %>%
  rownames_to_column("gene") %>%
  rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

LPS_bLPS_pnoLPS <- results_pseudo_noLPS_bulk_LPS_norm %>%
  filter(log2FoldChange > 1, padj < 0.01, !is.na(padj))


### Venn diagram ###

LPS_inducible_found_sets <- list(
  Bulk = LPS_bulk_only$gene,
  Pseudo = LPS_pseudo_only$gene,
  pLPS_bnoLPS = LPS_pLPS_bnoLPS$gene,
  pnoLPS_bLPS = LPS_bLPS_pnoLPS$gene
)

ggvenn(LPS_inducible_found_sets, 
       show_percentage = FALSE,
       set_name_size = 3.5, 
       text_size = 5)



####  HISTOGRAM AND STATS BEFORE / AFTER SCALING ####

#dplyr::select columns
dplyr::selected_cols <- names(merged_counts)[grepl("^(bulk_|pseudo_)", names(merged_counts))]

#Convert to long format
long_data <- merged_counts %>%
  dplyr::dplyr::select(all_of(dplyr::selected_cols)) %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "expression")

#Add method of each sample
long_data <- long_data %>%
  mutate(sample_type = ifelse(grepl("^bulk_", sample), "bulk", "pseudo"))

#Add unscaled and scaled values
long_data <- long_data %>%
  mutate(
    log1p_expression = log1p(expression),
    scaled_expression = expression * 100,               
    log1p_scaled_expression = log1p(scaled_expression)
  )

#Plot histograms of initial values
ggplot(long_data, aes(x = log1p_expression)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  facet_wrap(~sample, scales = "fixed") +
  labs(
    title = "log1p(Expression) Histogram by Sample",
    x = "log1p(counts)",
    y = "Frequency"
  ) +
  theme_minimal()

#Plot histograms of scaled values
ggplot(long_data, aes(x = log1p_scaled_expression)) +
  geom_histogram(bins = 100, fill = "darkorange", color = "black") +
  facet_wrap(~sample, scales = "fixed") +
  labs(
    title = "log1p(Counts × 10) Histogram by Sample",
    x = "log1p(Expression × 10)",
    y = "Frequency"
  ) +
  theme_minimal()

#Plot an overlayed histogram
ggplot(long_data, aes(x = log1p_scaled_expression, fill = sample_type)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity", color = "black") +
  facet_wrap(~sample_type, scales = "fixed") +
  labs(
    title = "Overlayed log1p(counts × 10) by Sample Type",
    x = "log1p(counts × 10)",
    y = "Frequency"
  ) +
  scale_fill_manual(values = c("bulk" = "#0072B2", "pseudo" = "#D55E00")) +
  theme_minimal()




#Summary statistics 
summary_stats <- long_data %>%
  summarise(
    mean = mean(log1p_expression, na.rm = TRUE),
    median = median(log1p_expression, na.rm = TRUE),
    sd = sd(log1p_expression, na.rm = TRUE),
    min = min(log1p_expression, na.rm = TRUE),
    max = max(log1p_expression, na.rm = TRUE),
    n = n()
  )

print(summary_stats)




# Reshape to long format for expression type
comparison_data <- long_data %>%
  dplyr::dplyr::select(sample, sample_type, log1p_expression, log1p_scaled_expression) %>%
  tidyr::pivot_longer(cols = c(log1p_expression, log1p_scaled_expression),
                      names_to = "expression_type",
                      values_to = "value") %>%
  dplyr::mutate(expression_type = dplyr::recode(expression_type,
                                                log1p_expression = "Before scaling",
                                                log1p_scaled_expression = "After scaling"))

ggplot(comparison_data, aes(x = value, fill = expression_type)) +
  geom_histogram(bins = 100, alpha = 0.5, position = "identity", color = "black") +
  facet_wrap(~sample_type) +
  scale_fill_manual(values = c("Before scaling" = "steelblue", "After scaling" = "darkorange")) +
  labs(
    title = "Comparison of Distributions Before and After Scaling",
    x = "log1p(Expression)",
    y = "Frequency",
    fill = "Dostribution"
  ) +
  theme_minimal()


####  GET DESEQ2 SEPARATE CONTRASTS  ####

#merge counts and gene lengths

# Prepare separate data frames
bulk_df <- bulk_for_diff %>%
  as.data.frame() %>%
  rownames_to_column("gene")

pseudo_df <- pseudobulk_for_diff %>%
  as.data.frame() %>%
  rownames_to_column("gene")

lengths_df <- merged_lengths %>%
  group_by(gene_id) %>%
  summarise(length = mean(length, na.rm = TRUE)) %>%
  rename(gene = gene_id)

#Merge by gene
merged_counts_gene_length <- bulk_df %>%
  inner_join(pseudo_df, by = "gene") %>%
  inner_join(lengths_df, by = "gene") %>%
  arrange(gene)


#Identify sample groups

bulk_samples <- merged_counts_gene_length %>%
  dplyr::select(dplyr::starts_with("bulk")) %>%
  colnames()

pseudo_samples <- merged_counts_gene_length %>%
  dplyr::select(dplyr::starts_with("pseudo")) %>%
  colnames()


#Normalise only bulk by gene length

merged_counts_gene_length <- merged_counts_gene_length %>%
  dplyr::mutate(across(all_of(bulk_samples), ~ .x / length)) %>%
  dplyr::select(-length)


#Scale and round to integers

merged_counts_gene_length <- merged_counts_gene_length %>%
  dplyr::mutate(across(all_of(c(bulk_samples, pseudo_samples)), ~ round(.x * 100))) %>%
  dplyr::mutate(across(all_of(c(bulk_samples, pseudo_samples)), ~ tidyr::replace_na(.x, 0)))


#Split the dataset for the two contrasts

pseudo_LPS_bulk_noLPS_for_DESeq2 <- merged_counts_gene_length %>%
  dplyr::select(-c("pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
                   "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL")) %>%
  column_to_rownames(var = "gene")

pseudo_noLPS_bulk_LPS_for_DESeq2 <- merged_counts_gene_length %>%
  dplyr::select(c("gene", "pseudo_tt", "pseudo_td", "pseudo_dt", "pseudo_dd",
                  "bulk_TLT", "bulk_TL", "bulk_DLT", "bulk_DL")) %>%
  column_to_rownames(var = "gene")


#Run DESeq2 for pseudo_LPS/bulk_noLPS

#Create metadata for the dataset, identifying LPS/noLPS groups

metaData_pLPS_bnoLPS_for_DESeq2 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("bulk_TT", "bulk_TD", "bulk_DT", "bulk_DD", 
              "pseudo_dl","pseudo_dlt", "pseudo_tl",  "pseudo_tlt")) %>%
  column_to_rownames(var = "Samples")


##  RUN DESeq2 AND GET THE LPS-BIASED GENES  ##


dds_pLPSbnoLPS <- DESeqDataSetFromMatrix(countData = pseudo_LPS_bulk_noLPS_for_DESeq2, 
                                         colData = metaData_pLPS_bnoLPS_for_DESeq2, 
                                         design = ~ Treatment)

dds_pLPSbnoLPS <- DESeq(dds_pLPSbnoLPS)

results_pLPSbnoLPS <- results(dds_pLPSbnoLPS, contrast = c("Treatment", "LPS", "noLPS")) %>%
  as.data.frame()


#Filter genes with log2FC > 1, padj < 0.01, and non-NA padj

LPS_pLPS_bnoLPS_DESeq2 <- results_pLPSbnoLPS %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj)) %>%
  rownames_to_column(var = "gene")


# Run DESeq2 for pseudo_noLPS/bulk_LPS


#Create metadata for the dataset, identifying LPS/noLPS groups

metaData_pnoLPS_bLPS_for_DESeq2 <- data.frame(
  Treatment = c("noLPS", "noLPS", "noLPS", "noLPS", "LPS", "LPS", "LPS", "LPS"),
  Samples = c("pseudo_tt","pseudo_td","pseudo_dt","pseudo_dd",
              "bulk_TLT","bulk_TL","bulk_DLT","bulk_DL" )) %>%
  column_to_rownames(var = "Samples")



##  RUN DESeq2 AND GET THE LPS-BIASED GENES  ##


dds_pnoLPSbLPS <- DESeqDataSetFromMatrix(countData = pseudo_noLPS_bulk_LPS_for_DESeq2, 
                                         colData = metaData_pnoLPS_bLPS_for_DESeq2, 
                                         design = ~ Treatment)

dds_pnoLPSbLPS <- DESeq(dds_pnoLPSbLPS)

results_pnoLPSbLPS <- results(dds_pnoLPSbLPS, contrast = c("Treatment", "LPS", "noLPS")) %>%
  as.data.frame()


#Filter genes with log2FC > 1, padj < 0.01, and non-NA padj

LPS_pnoLPS_bLPS_DESeq2 <- results_pnoLPSbLPS %>%
  filter(log2FoldChange > 1 & padj < 0.01 & !is.na(padj)) %>%
  rownames_to_column(var = "gene")


####  MAKE VENN DIAGRAMS FOR BULK/PSEUDO AND CROSS-CONTRASTS  ####

# Prepare named lists for the three venn diagrams

venn1 <- list(
  "bulk / bulk"   = LPS_bulk_only$gene,
  "pseudo / pseudo" = LPS_pseudo_only$gene
)

venn2 <- list(
  `bulk LPS / pseudo noLPS` = as.character(LPS_pLPS_bnoLPS_DESeq2$gene),
  `pseudo LPS / bulk noLPS` = as.character(LPS_pnoLPS_bLPS_DESeq2$gene)
)

venn3 <- list(
  "bulk / bulk"      = LPS_bulk_only$gene,
  "pseudo / pseudo"  = LPS_pseudo_only$gene,
  "bulk LPS / pseudo noLPS" = sig_genes_A,
  "pseudo LPS / bulk noLPS" = sig_genes_B
)

# 2. Function to calculate % confirmed genes from the known LPS-inducible

pct_confirmed <- function(found, known){
  n_found <- length(found)
  n_intersect <- length(intersect(found, known))
  round(n_intersect / n_found * 100, 2)
}

known <- unique(na.omit(LPS_inducible_genes$MoMf_LPS_inducible_genes))

# 3. Plot the Venn diagrams

# Venn 1: bulk vs pseudo
ggvenn(venn1, 
       fill_color = c("blue", "red"),
       show_percentage = FALSE,
       set_name_size = 6,
       text_size = 6) +
  ggtitle("Bulk vs Pseudobulk DESeq2")

# Venn 2: cross DESeq2
ggvenn(
  venn2,
  fill_color = c("green", "magenta"),
  show_percentage = FALSE,
  set_name_size = 6,
  text_size = 8,
)



####  MAKE SUMMARY TABLE WITH % CONFIRMED GENES PER GROUP  ####

summary_tbl <- tibble(
  `LPS‐biased genes` = c("bulk / bulk",
                         "pseudo / pseudo",
                         "bulk LPS / pseudo noLPS",
                         "pseudo LPS / bulk noLPS"),
  `% confirmed` = c(
    pct_confirmed(LPS_bulk_only$gene,   known),
    pct_confirmed(LPS_pseudo_only$gene, known),
    pct_confirmed(sig_genes_A,          known),
    pct_confirmed(sig_genes_B,          known)
  )
)

print(summary_tbl)

