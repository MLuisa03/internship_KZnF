library(dplyr)
library(tidyverse)
library(tibble)
library(DESeq2)
library(viridis)
library(ComplexHeatmap)
library(grid)
library(circlize)



####  HIERARCHICAL HEATMAP OF VARIANCE-STABILISED COUNTS  ####

#Get all counts

merged_counts <- read.delim("counts_of_all_bulk_samples")

colnames(merged_counts) <- c("gene", "iHSC_1", "adult_mono_1", "adult_mono_2", "adult_mono_3", "adult_mono_4", 
                             "adult_mono_5", "adult_mono_6", "adult_mono_7", "adult_mono_8", "adult_mono_9", 
                             "adult_mono_10", "adult_mono_11", "iMye_1", "iMye_2", "iMye_3",
                             "iMye_4", "iMye_5", "iMye_6", "iMye_7", "iPSC_1", "iPSC_2", "iPSC_3", "iPSC_4",
                             "adult_mono_12", "iMye_8", "iHSC_2", "adult_mono_13", "iPSC_5", "iHSC_3", 
                             "iHSC_4", "iPSC_6", "iHSC_5", "adult_mono_14", "iHSC_6", "iHSC_7", "iHSC_8", 
                             "iHSC_9", "iHSC_10")

merged_counts <- column_to_rownames(merged_counts, var = "gene")

# Create DESeq2 dataset


metadata_merged_counts <- data.frame(
  Sample = c("iHSC_1", "adult_mono_1", "adult_mono_2", "adult_mono_3", "adult_mono_4", 
             "adult_mono_5", "adult_mono_6", "adult_mono_7", "adult_mono_8", "adult_mono_9", 
             "adult_mono_10", "adult_mono_11", "iMye_1", "iMye_2", "iMye_3",
             "iMye_4", "iMye_5", "iMye_6", "iMye_7", "iPSC_1", "iPSC_2", "iPSC_3", "iPSC_4",
             "adult_mono_12", "iMye_8", "iHSC_2", "adult_mono_13", "iPSC_5", "iHSC_3", 
             "iHSC_4", "iPSC_6", "iHSC_5", "adult_mono_14", "iHSC_6", "iHSC_7", "iHSC_8", 
             "iHSC_9", "iHSC_10"),
  Stage = c("iHSC", "adult_mono", "adult_mono", "adult_mono", "adult_mono", 
            "adult_mono", "adult_mono", "adult_mono", "adult_mono", "adult_mono", 
            "adult_mono", "adult_mono", "iMye", "iMye", "iMye",
            "iMye", "iMye", "iMye", "iMye", "iPSC", "iPSC", "iPSC", "iPSC",
            "adult_mono", "iMye", "iHSC", "adult_mono", "iPSC", "iHSC", 
            "iHSC", "iPSC", "iHSC", "adult_mono", "iHSC", "iHSC", "iHSC", 
            "iHSC", "iHSC"),
  Treatment = c("DMSO", "TA", "DMSO", "DMSO", "DMSO",
                "DMSO", "TA", "DEX", "PRED", "DMSO",
                "TA", "DEX", "DMSO", "TA", "DEX",
                "PRED", "DMSO", "TA", "DEX", "DMSO", "TA", "TA", "TA",
                "TA", "DLT", "DMSO", "DMSO", "DMSO", "DMSO",
                "DMSO", "TA", "DMSO", "DMSO", "TA", "TA", "DMSO",
                "DMSO", "DMSO")) %>% 
  column_to_rownames(var = "Sample")


dds <- DESeqDataSetFromMatrix(countData = merged_counts,
                              colData = metadata_merged_counts,
                              design = ~ Stage)


# Run DESeq to estimate dispersions 
dds <- DESeq(dds)

# Variance stabilizing transformation

#blind = TRUE, unguided by design

vsd_blind <- vst(dds, blind = TRUE)

vsd_mat_blind <- assay(vsd_blind)

KZnF_vsd_blind <- as.data.frame(vsd_mat_blind) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(rowSums(.[, -1]) > 0) %>%  
  inner_join(krab_krueppel_genes, by = "gene") %>% 
  column_to_rownames(var = "gene") %>% 
  { .[apply(., 1, var) > 0, ] } %>%  
  as.matrix()



# Save pheatmap to file
png("heatmap_hierarchical_tall.png", width = 8, height = 10, units = "in", res = 600, bg = "transparent")

pheatmap(KZnF_vsd_blind,
         scale = "row",
         main = "Heatmap KZnF",
         show_rownames = FALSE,
         color = viridis(100))

dev.off()


####  KMEANS HEATMAP WITH ROW SCALING DONE ONLY WITH iPSCs AND ADULT MONO  ####


##Define the subset of columns to use for clustering
cols_for_scale <- c("iPSC_1", "iPSC_2", "iPSC_3", "iPSC_4", "iPSC_5", "iPSC_6",
                    "adult_mono_1", "adult_mono_2", "adult_mono_3", "adult_mono_4", "adult_mono_5", "adult_mono_6", "adult_mono_7",
                    "adult_mono_8", "adult_mono_9", "adult_mono_10", "adult_mono_11", "adult_mono_12", "adult_mono_13", "adult_mono_14")

##Compute row means and SDs based on the subset
row_means <- rowMeans(KZnF_vsd_blind[, cols_for_scale])
row_sds   <- apply(KZnF_vsd_blind[, cols_for_scale], 1, sd)
invariant <- row_sds == 0
row_sds[invariant] <- 1

## Z-score scaling across all columns using means/SDs from the subset
mat_scaled <- sweep(KZnF_vsd_blind, 1, row_means, "-")
mat_scaled <- sweep(mat_scaled, 1, row_sds, "/")
mat_scaled[invariant, ] <- 0

##Determine optimal k using Elbow Method
wss <- sapply(1:10, function(k) {
  set.seed(123)
  kmeans(mat_scaled[, cols_for_scale], centers = k, nstart = 10)$tot.withinss
})

# Optional: plot elbow curve
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares")

##Pick optimal k (based on elbow plot)
optimal_k <- 2  

## Perform k-means clustering
set.seed(123)
km <- kmeans(mat_scaled[, cols_for_scale], centers = optimal_k)
row_clusters <- km$cluster

##Column reordering according to stage
desired_column_order <- c(
  "iPSC_1", "iPSC_2", "iPSC_3", "iPSC_4", "iPSC_5", "iPSC_6",
  "iHSC_1", "iHSC_2", "iHSC_3", "iHSC_4", "iHSC_5", "iHSC_6", "iHSC_7", "iHSC_8", "iHSC_9", "iHSC_10",
  "iMye_1", "iMye_2", "iMye_3", "iMye_4", "iMye_5", "iMye_6", "iMye_7", "iMye_8",
  "adult_mono_1", "adult_mono_2", "adult_mono_3", "adult_mono_4", "adult_mono_5", "adult_mono_6", "adult_mono_7",
  "adult_mono_8", "adult_mono_9", "adult_mono_10", "adult_mono_11", "adult_mono_12", "adult_mono_13", "adult_mono_14"
)

##Color palette and annotation 

breaks_values <- seq(min(mat_scaled), max(mat_scaled), length.out = 256)

# Create the color function with specified breaks and the viridis palette
col_fun <- colorRamp2(
  breaks = breaks_values,  
  colors = viridis(256) 
)
cluster_cols <- RColorBrewer::brewer.pal(3, "Set2")[1:optimal_k]  # Only take the first 2 colors

##Ensure that 'Cluster' is a named vector
row_ha <- rowAnnotation(
  Cluster = factor(row_clusters),
  col     = list(Cluster = setNames(cluster_cols, 1:optimal_k)),
  show_annotation_name = TRUE
)

## heatmap
ht <- Heatmap(
  mat_scaled[, desired_column_order],
  name           = "z-score",
  col            = col_fun,
  show_row_names   = FALSE,
  show_column_names = TRUE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_split        = row_clusters,
  left_annotation  = row_ha,
  heatmap_legend_param = list(title = "Scaled\nexpression\n(z-score)")
)

grid::grid.newpage()
draw(ht)

