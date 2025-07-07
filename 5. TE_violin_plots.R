library(ggplot2)
library(GGally)
library(dplyr)



####  LOAD TE TARGETS AND EXAMINE VARIABLES  ####

#Load the dataset 

TE_targets <- read.delim("Supplemental_Table_S3_de_Tribolet_Hardy_2023")


#Examine co-variance of the variables: n_peaks, expected_peaks and weight

pairs(TE_targets[, c("n.peaks.on.subfamily", "expected.number.of.peaks", "weight")],
      main = "Pair Plot of: n.peaks, expected peaks, and weight",
      pch = 19,          
      col = "steelblue") 

ggpairs(TE_targets[, c("n.peaks.on.subfamily", "expected.number.of.peaks", "weight")])



####  GET SELECTED GENES AND CALCULATE NORMALISED PEAKS  ####

## Get important genes as vectors (from the previous script )  ### NEEDS TO BE RUN WITH PREVIOUS SCRIPT

adult_b_vs_iPSC <- as.data.frame(adult_biased_vs_iPSC)
adult_b_vs_iPSC <- adult_b_vs_iPSC$gene
iPSC_b_vs_adult <- as.data.frame(iPSC_biased_vs_adult)
iPSC_b_vs_adult <- iPSC_b_vs_adult$gene


# Calculate the normalised peaks

TE_targets$TE_peaks <- TE_targets$n.peaks.on.subfamily / TE_targets$expected.number.of.peaks
TE_targets$TE_peaks_weight <- (TE_targets$n.peaks.on.subfamily / TE_targets$expected.number.of.peaks) * TE_targets$weight 


# Start with empty list to collect data for each group

grouped_rows <- list()


# For each group, subset relevant rows and assign group label

grouped_rows[["iHSC_biased_vs_iPSC"]] <- TE_targets[TE_targets$KZFP %in% iHSC_biased_vs_iPSC, ]
grouped_rows[["iMye_biased_vs_iPSC"]] <- TE_targets[TE_targets$KZFP %in% iMye_biased_vs_iPSC, ]
grouped_rows[["adult_biased_vs_iPSC"]] <- TE_targets[TE_targets$KZFP %in% adult_biased_vs_iPSC, ]
grouped_rows[["iPSC_biased_vs_adult"]] <- TE_targets[TE_targets$KZFP %in% iPSC_biased_vs_adult, ]
grouped_rows[["iPSC_biased_vs_iHSC"]] <- TE_targets[TE_targets$KZFP %in% iPSC_biased_vs_iHSC, ]
grouped_rows[["iPSC_biased_vs_iMye"]] <- TE_targets[TE_targets$KZFP %in% iPSC_biased_vs_iMye, ]

grouped_rows[["iHSC_biased_vs_iMye"]] <- TE_targets[TE_targets$KZFP %in% iHSC_biased_vs_iMye, ]
grouped_rows[["iMye_biased_vs_iHSC"]] <- TE_targets[TE_targets$KZFP %in% iMye_biased_vs_iHSC, ]
grouped_rows[["iMye_biased_vs_adult"]] <- TE_targets[TE_targets$KZFP %in% iMye_biased_vs_adult, ]
grouped_rows[["adult_biased_vs_iMye"]] <- TE_targets[TE_targets$KZFP %in% adult_biased_vs_iMye, ]
grouped_rows[["iHSC_biased_vs_adult"]] <- TE_targets[TE_targets$KZFP %in% iHSC_biased_vs_adult, ]
grouped_rows[["adult_biased_vs_iHSC"]] <- TE_targets[TE_targets$KZFP %in% adult_biased_vs_iHSC, ]



# For each subset, add a 'group' column

for (grp in names(grouped_rows)) {
  if (nrow(grouped_rows[[grp]]) > 0) {
    grouped_rows[[grp]]$group <- grp}}


# Combine all subsets into one long-format dataframe

TE_targets_long <- do.call(rbind, grouped_rows)

TE_targets_long_new <- TE_targets_long %>%
  filter(is.finite(TE_peaks), TE_peaks > 0)

TE_targets_long_new$TE_peaks <- TE_targets_long_new$TE_peaks + 1e-6


####  VISUALISE PEAKS AS VIOLIN PLOTS  #### 

# Set order of groups

TE_targets_long$group <- factor(
  TE_targets_long$group,
  levels = c(
    "iPSC_biased_vs_iHSC",
    "iPSC_biased_vs_iMye",
    "iPSC_biased_vs_adult",
    "iHSC_biased_vs_iPSC",
    "iMye_biased_vs_iPSC",
    "adult_biased_vs_iPSC"))


# Define colours

group_colors <- c(
  "iHSC_biased_vs_iPSC" = "#007090",
  "iMye_biased_vs_iPSC" = "#007F5F",
  "adult_biased_vs_iPSC" = "#B35794",
  "iPSC_biased_vs_adult" = "#CC8400",
  "iPSC_biased_vs_iHSC" = "#D98800",
  "iPSC_biased_vs_iMye" = "#F0B94E")


# Create violin plot

p <- ggplot(TE_targets_long, aes(x = group, y = TE_peaks, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal(base_size = 20) +
  labs(
    y = "Log10(TE peaks / expected)",
    x = "Group",
    title = "Log-Scaled Comparison of TE peaks Across Groups"
  ) +
  scale_fill_manual(values = group_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "bold"),
    panel.grid.major = element_line(color = "black"),
    panel.grid.minor = element_line(color = "black"),
    legend.position = "none")


# Save plot

ggsave("TE_violin_plot_blackgrid.png", plot = p, dpi = 600, width = 8, height = 10, bg = "transparent")
