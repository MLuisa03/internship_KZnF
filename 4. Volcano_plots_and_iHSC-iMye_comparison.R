library(dplyr)
library(tibble)
library(ggplot2)



#### LOAD DIFFERENTIAL EXPRESSION RESULTS FROM SEQ2SCIENCE  ####

iMye_iPSC <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_iMyebulk_iPSC.diffexp.tsv",
                        sep = "\t", header = TRUE)

iMo_iHSCpseudo <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage_iMo_iHSCpseudo.diffexp.tsv",
                             sep = "\t", header = TRUE)

iMo2_iMo6 <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage3_iMo2_iMo6.diffexp.tsv",
                        sep = "\t", header = TRUE)

adult_iHSCbulk <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_adult_iHSCbulk.diffexp.tsv",
                             sep = "\t", header = TRUE)

adult_iMyebulk <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_adult_iMyebulk.diffexp.tsv",
                             sep = "\t", header = TRUE)

adult_iPSC <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_adult_iPSC.diffexp.tsv",
                         sep = "\t", header = TRUE)

iHSCbulk_iPSC <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_iHSCbulk_iPSC.diffexp.tsv",
                            sep = "\t", header = TRUE)

iMo_iMye <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage_iMo_iMye.diffexp.tsv",
                       sep = "\t", header = TRUE)

iMyebulk_iHSCbulk <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-stage1_iMyebulk_iHSCbulk.diffexp.tsv",
                                sep = "\t", header = TRUE)

iMye_iHSC_pseudo <- read.delim("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Seq2science output/All_contrasts_latest/hg38-Day19id_iMye_iHSCpseudo.diffexp.tsv",
                               sep = "\t", header = TRUE)

krab_krueppel_genes <- read.csv("C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Inventory datasets/overlap_databases_Trono2023_KRABopedia.csv")



####  CLEAN DATASETS AND EXTRAPOLATE KZNF  ####

#iMye_iPSC
iMye_iPSC <- iMye_iPSC %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMye_iPSC <- dplyr::inner_join(iMye_iPSC, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iMo_iHSCpseudo
iMo_iHSCpseudo <- iMo_iHSCpseudo %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMo_iHSCpseudo <- dplyr::inner_join(iMo_iHSCpseudo, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iMo2_iMo6
iMo2_iMo6 <- iMo2_iMo6 %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMo2_iMo6 <- dplyr::inner_join(iMo2_iMo6, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#adult_iHSCbulk
adult_iHSCbulk <- adult_iHSCbulk %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_adult_iHSCbulk <- dplyr::inner_join(adult_iHSCbulk, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#adult_iMyebulk
adult_iMyebulk <- adult_iMyebulk %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_adult_iMyebulk <- dplyr::inner_join(adult_iMyebulk, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#adult_iPSC
adult_iPSC <- adult_iPSC %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_adult_iPSC <- dplyr::inner_join(adult_iPSC, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iHSCbulk_iPSC
iHSCbulk_iPSC <- iHSCbulk_iPSC %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iHSCbulk_iPSC <- dplyr::inner_join(iHSCbulk_iPSC, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iMo_iMye
iMo_iMye <- iMo_iMye %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMo_iMye <- dplyr::inner_join(iMo_iMye, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iMyebulk_iHSCbulk
iMyebulk_iHSCbulk <- iMyebulk_iHSCbulk %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMyebulk_iHSCbulk <- dplyr::inner_join(iMyebulk_iHSCbulk, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")

#iMye_iHSCpseudo
iMye_iHSC_pseudo <- iMye_iHSC_pseudo %>%
  na.omit() %>%
  filter(rowSums(. != 0) > 0)
KZnF_iMye_iHSC_pseudo <- dplyr::inner_join(iMye_iHSC_pseudo, krab_krueppel_genes, by = join_by("X" == "gene")) %>%
  column_to_rownames(var = "X")



####  ESTABLISH VOLCANO PLOT FUNCTION  ####

volcano_plot <- function(data, 
                         baseline_label, 
                         up_biased_label, 
                         baseline_color, 
                         up_biased_color, 
                         baseline_shape, 
                         up_biased_shape, 
                         title) {
  
  #Replace pvalues of 0 to avoid infinite values when taking -log10(pvalue)
  data$pvalue[data$pvalue == 0] <- 1e-10
  
  #Define categories
  data$category <- "Not Significant"
  data$category[data$padj < 0.01] <- "Significant"
  data$category[data$padj < 0.01 & data$log2FoldChange > 1] <- up_biased_label
  data$category[data$padj < 0.01 & data$log2FoldChange < -1] <- baseline_label
  
  # Create the combined factor with specified levels
  data$combined <- factor(data$category, 
                          levels = c("Significant", up_biased_label, 
                                     baseline_label, "Not Significant"))
  
  #Define colors and shapes of the datapoints
  colors <- setNames(
    c("#0072B2", up_biased_color, baseline_color, "grey"),
    c("Significant", up_biased_label, baseline_label, "Not Significant")
  )
  
  shapes <- setNames(
    c(16, up_biased_shape, baseline_shape, 1),
    c("Significant", up_biased_label, baseline_label, "Not Significant")
  )
  
  #Get category counts and create labels 
  category_counts <- data %>%
    group_by(category) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(category, " (", n, ")"))
  labels_with_counts <- setNames(category_counts$label, category_counts$category)
  
  #Make the volcano plot
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), 
                   color = combined, shape = combined)) +
    geom_point(alpha = 0.8, size = 2.5) +  
    scale_color_manual(values = colors, labels = labels_with_counts) + 
    scale_shape_manual(values = shapes, labels = labels_with_counts) +
    labs(title = title,
         subtitle = "Genes were considered biased for a cell type if l2FC>1",
         x = "log2 Fold Change",
         y = "-log10(padj)",
         color = "Genes",   
         shape = "Genes") +  
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 60)) + 
    theme_minimal() +  
    theme(aspect.ratio = 1)
}

####  ESTABLISH ALTERNATIVE AESTHETIC VOLCANO FUNCTION  ####

volcano_plot_alt <- function(data, 
                             baseline_label, 
                             up_biased_label, 
                             baseline_color, 
                             up_biased_color, 
                             baseline_shape, 
                             up_biased_shape, 
                             title) {
  
  #Replace pvalues of 0 to avoid infinite values when taking -log10(pvalue)
  data$pvalue[data$pvalue == 0] <- 1e-10
  
  #Define categories
  data$category <- "Not Significant"
  data$category[data$padj < 0.01] <- "Invariant"
  data$category[data$padj < 0.01 & data$log2FoldChange > 1] <- up_biased_label
  data$category[data$padj < 0.01 & data$log2FoldChange < -1] <- baseline_label
  
  # Create the combined factor with specified levels
  data$combined <- factor(data$category, 
                          levels = c("Invariant", up_biased_label, 
                                     baseline_label, "Not Significant"))
  
  #Define colors and shapes of the datapoints
  colors <- setNames(
    c("#4D4D4D", up_biased_color, baseline_color, "grey"),
    c("Invariant", up_biased_label, baseline_label, "Not Significant")
  )
  
  shapes <- setNames(
    c(16, up_biased_shape, baseline_shape, 1),
    c("Invariant", up_biased_label, baseline_label, "Not Significant")
  )
  
  #Get category counts and create labels 
  category_counts <- data %>%
    group_by(category) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(category, " (", n, ")"))
  labels_with_counts <- setNames(category_counts$label, category_counts$category)
  
  #Make the volcano plot
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), 
                   color = combined, shape = combined)) +
    geom_point(alpha = 0.8, size = 4) +  
    scale_color_manual(values = colors, labels = labels_with_counts) + 
    scale_shape_manual(values = shapes, labels = labels_with_counts) +
    labs(title = title,
         x = "log2 Fold Change",
         y = "-log10(padj)",
         color = "Genes",   
         shape = "Genes") +  
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 60)) + 
    theme_minimal() +  
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 24),
      plot.title = element_text(size = 30),
      axis.title = element_text(size = 24),       
      axis.text = element_text(size = 20),
      panel.grid.major = element_line(color = "black", linewidth = 0.2),
      panel.grid.minor = element_line(color = "black", linewidth = 0.1)
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    )
}

####  VOLCANO PLOTS  ####

#iHSC_bulk / iPSC
volcano_plot(KZnF_iHSCbulk_iPSC, 
             baseline_label = "iPSC-biased", 
             up_biased_label = "iHSC-biased", 
             baseline_color = "#E69F00", 
             up_biased_color = "#4C8BF5", 
             baseline_shape = 15, 
             up_biased_shape = 25, 
             title = "iHSC vs iPSC (bulk) - KZnF genes")


#iMye / iPSC
volcano_plot(KZnF_iMye_iPSC, 
             baseline_label = "iPSC-biased", 
             up_biased_label = "iMye-biased", 
             baseline_color = "#E69F00", 
             up_biased_color = "#009E73", 
             baseline_shape = 15, 
             up_biased_shape = 17, 
             title = "iMye vs iPSC (bulk) - KZnF genes")

#adult / iPSC
volcano_plot(KZnF_adult_iPSC, 
             baseline_label = "iPSC-biased", 
             up_biased_label = "adult_mono-biased", 
             baseline_color = "#E69F00", 
             up_biased_color = "#CC79A7", 
             baseline_shape = 15, 
             up_biased_shape = 8, 
             title = "adult monocytes vs iPSC (bulk) - KZnF genes")

#iMye / iHSC (bulk)
volcano_plot(KZnF_iMyebulk_iHSCbulk, 
             baseline_label = "iHSC-biased", 
             up_biased_label = "iMye-biased", 
             baseline_color = "#4C8BF5", 
             up_biased_color = "#009E73", 
             baseline_shape = 25, 
             up_biased_shape = 17, 
             title = "iHSC vs iMye (bulk) - KZnF genes")

#iMye / iHSC (pseudobulk)
volcano_plot(KZnF_iMye_iHSC_pseudo, 
             baseline_label = "iHSC-biased", 
             up_biased_label = "iMye-biased", 
             baseline_color = "#4C8BF5", 
             up_biased_color = "#009E73", 
             baseline_shape = 25, 
             up_biased_shape = 17, 
             title = "iHSC vs iMye (pseudo) - KZnF genes")

#adult / iHSC
volcano_plot(KZnF_adult_iHSCbulk, 
             baseline_label = "iHSC-biased", 
             up_biased_label = "adult_mono-biased", 
             baseline_color = "#4C8BF5", 
             up_biased_color = "#CC79A7", 
             baseline_shape = 25, 
             up_biased_shape = 8, 
             title = "adult monocytes vs iHSC (bulk) - KZnF genes")

#adult / iMye
volcano_plot(KZnF_adult_iMyebulk, 
             baseline_label = "iMye-biased", 
             up_biased_label = "adult_mono-biased", 
             baseline_color = "#009E73", 
             up_biased_color = "#CC79A7", 
             baseline_shape = 17, 
             up_biased_shape = 8, 
             title = "adult monocytes vs iMye (bulk) - KZnF genes")

#iMo / iHSC
volcano_plot(KZnF_iMo_iHSCpseudo, 
             baseline_label = "iHSC-biased", 
             up_biased_label = "iMo-biased", 
             baseline_color = "#4C8BF5", 
             up_biased_color = "turquoise3", 
             baseline_shape = 25, 
             up_biased_shape = 18, 
             title = "iMo vs iHSC (pseudo) - KZnF genes")

#iMo / iMye - KZnF genes
volcano_plot(KZnF_iMo_iMye, 
             baseline_label = "iMye-biased", 
             up_biased_label = "iMo-biased", 
             baseline_color = "#009E73", 
             up_biased_color = "turquoise3", 
             baseline_shape = 17, 
             up_biased_shape = 18, 
             title = "iMo vs iMye (pseudo) - KZnF genes")

#iMo / iMye - all genes
volcano_plot(iMo_iMye, 
             baseline_label = "iMye-biased", 
             up_biased_label = "iMo-biased", 
             baseline_color = "#009E73", 
             up_biased_color = "turquoise3", 
             baseline_shape = 17, 
             up_biased_shape = 18, 
             title = "iMo vs iMye (pseudo) - all genes")

#iMo2 / iMo6
volcano_plot(KZnF_iMo2_iMo6, 
             baseline_label = "iMo_6-biased", 
             up_biased_label = "iMo_2-biased", 
             baseline_color = "turquoise4", 
             up_biased_color = "turquoise", 
             baseline_shape = 18, 
             up_biased_shape = 18, 
             title = "iMo group 2 vs iMo group 6 (pseudo) - KZnF genes")



####  ALTERNATIVE VOLCANO PLOTS  ####

#iHSC_bulk / iPSC
volcano_plot_alt(KZnF_iHSCbulk_iPSC, 
                 baseline_label = "iPSC-biased", 
                 up_biased_label = "iHSC-biased", 
                 baseline_color = "#D98800", 
                 up_biased_color = "#005086", 
                 baseline_shape = 15, 
                 up_biased_shape = 17, 
                 title = "iHSC vs iPSC (bulk) - KZnF genes")

ggsave("volcano_iHSC_iPSC.png", dpi = 600, width = 8, height = 8, bg = "transparent")

#adult / iPSC
volcano_plot_alt(KZnF_adult_iPSC, 
                 baseline_label = "iPSC-biased", 
                 up_biased_label = "adult_mono-biased", 
                 baseline_color = "#CC8400", 
                 up_biased_color = "#B35794", 
                 baseline_shape = 15, 
                 up_biased_shape = 18, 
                 title = "aMo vs iPSC (bulk) - KZnF genes")

ggsave("volcano_adult_iPSC.png", dpi = 600, width = 8.5, height = 8.5, bg = "transparent")


#iMye / iHSC (bulk)
volcano_plot_alt(KZnF_iMyebulk_iHSCbulk, 
                 baseline_label = "iHSC-biased", 
                 up_biased_label = "iMye-biased", 
                 baseline_color = "#66A9D3", 
                 up_biased_color = "#33B88A", 
                 baseline_shape = 17, 
                 up_biased_shape = 25, 
                 title = "iMye vs iHSC (bulk) - KZnF genes")

ggsave("volcano_iMye_iHSC.png", dpi = 600, width = 8, height = 8, bg = "transparent")


#iMo / iMye - KZnF genes
volcano_plot_alt(KZnF_iMo_iMye, 
                 baseline_label = "iMye-biased", 
                 up_biased_label = "iMo-biased", 
                 baseline_color = "#A6761D", 
                 up_biased_color = "turquoise3", 
                 baseline_shape = 25, 
                 up_biased_shape = 19, 
                 title = "iMo vs iMye (pseudo) - KZnF genes")

ggsave("volcano_iMo_iMye.png", dpi = 600, width = 8.2, height = 8.2, bg = "transparent")



#adult / iMye
volcano_plot_alt(KZnF_adult_iMyebulk, 
                 baseline_label = "iMye-biased", 
                 up_biased_label = "adult_mono-biased", 
                 baseline_color = "#339980", 
                 up_biased_color = "#A94C82", 
                 baseline_shape = 20, 
                 up_biased_shape = 18, 
                 title = "aMo vs iMye (bulk) - KZnF genes")

ggsave("volcano_aMo_iMye.png", dpi = 600, width = 8.5, height = 8.5, bg = "transparent")




####  ALTERNATIVE AESTHETIC VOLCANO PLOTS  ####

#iHSC_bulk / iPSC
volcano_plot_alt(KZnF_iHSCbulk_iPSC, 
                 baseline_label = "iPSC-biased", 
                 up_biased_label = "iHSC-biased", 
                 baseline_color = "#D98800", 
                 up_biased_color = "#005086", 
                 baseline_shape = 15, 
                 up_biased_shape = 17, 
                 title = "iHSC vs iPSC (bulk) - KZnF genes")

ggsave("volcano_iHSC_iPSC.png", dpi = 600, width = 8, height = 8, bg = "transparent")

#adult / iPSC
volcano_plot_alt(KZnF_adult_iPSC, 
                 baseline_label = "iPSC-biased", 
                 up_biased_label = "adult_mono-biased", 
                 baseline_color = "#CC8400", 
                 up_biased_color = "#B35794", 
                 baseline_shape = 15, 
                 up_biased_shape = 18, 
                 title = "aMo vs iPSC (bulk) - KZnF genes")

ggsave("volcano_adult_iPSC.png", dpi = 600, width = 8.5, height = 8.5, bg = "transparent")


#iMye / iHSC (bulk)
volcano_plot_alt(KZnF_iMyebulk_iHSCbulk, 
                 baseline_label = "iHSC-biased", 
                 up_biased_label = "iMye-biased", 
                 baseline_color = "#66A9D3", 
                 up_biased_color = "#33B88A", 
                 baseline_shape = 17, 
                 up_biased_shape = 25, 
                 title = "iMye vs iHSC (bulk) - KZnF genes")

ggsave("volcano_iMye_iHSC.png", dpi = 600, width = 8, height = 8, bg = "transparent")


#iMo / iMye - KZnF genes
volcano_plot_alt(KZnF_iMo_iMye, 
                 baseline_label = "iMye-biased", 
                 up_biased_label = "iMo-biased", 
                 baseline_color = "#A6761D", 
                 up_biased_color = "turquoise3", 
                 baseline_shape = 25, 
                 up_biased_shape = 19, 
                 title = "iMo vs iMye (pseudo) - KZnF genes")

ggsave("volcano_iMo_iMye.png", dpi = 600, width = 8.2, height = 8.2, bg = "transparent")



#adult / iMye
volcano_plot_alt(KZnF_adult_iMyebulk, 
                 baseline_label = "iMye-biased", 
                 up_biased_label = "adult_mono-biased", 
                 baseline_color = "#339980", 
                 up_biased_color = "#A94C82", 
                 baseline_shape = 20, 
                 up_biased_shape = 18, 
                 title = "aMo vs iMye (bulk) - KZnF genes")

ggsave("volcano_aMo_iMye.png", dpi = 600, width = 8.5, height = 8.5, bg = "transparent")



####  RETRIEVE ALL GENES OF INTEREST AND MERGE THEM TO THE INVENTORY  ####

iPSC_biased_vs_iHSC <- KZnF_iHSCbulk_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1) %>%
  pull(gene)

iHSC_biased_vs_iPSC <- KZnF_iHSCbulk_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1) %>%
  pull(gene)


iPSC_biased_vs_iMye <- KZnF_iMye_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1) %>%
  pull(gene)

iMye_biased_vs_iPSC <- KZnF_iMye_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1) %>%
  pull(gene)


iPSC_biased_vs_adult <- KZnF_adult_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1) %>%
  pull(gene)

adult_biased_vs_iPSC <- KZnF_adult_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1)%>%
  pull(gene)


iHSC_biased_vs_adult <- KZnF_adult_iHSCbulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1)%>%
  pull(gene)

adult_biased_vs_iHSC <- KZnF_adult_iHSCbulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1)%>%
  pull(gene)


iMye_biased_vs_adult <- KZnF_adult_iMyebulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1)%>%
  pull(gene)

adult_biased_vs_iMye <- KZnF_adult_iMyebulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1)%>%
  pull(gene)


iHSC_biased_vs_iMye <- KZnF_iMyebulk_iHSCbulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange < -1)%>%
  pull(gene)

iMye_biased_vs_iHSC <- KZnF_iMyebulk_iHSCbulk %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > 1)%>%
  pull(gene)


# Merge to the inventory

gene_sets <- list(
  iPSC_biased_vs_iHSC   = iPSC_biased_vs_iHSC,
  iHSC_biased_vs_iPSC   = iHSC_biased_vs_iPSC,
  iPSC_biased_vs_iMye   = iPSC_biased_vs_iMye,
  iMye_biased_vs_iPSC   = iMye_biased_vs_iPSC,
  iPSC_biased_vs_adult  = iPSC_biased_vs_adult,
  adult_biased_vs_iPSC  = adult_biased_vs_iPSC,
  iHSC_biased_vs_adult  = iHSC_biased_vs_adult,
  adult_biased_vs_iHSC  = adult_biased_vs_iHSC,
  iMye_biased_vs_adult  = iMye_biased_vs_adult,
  adult_biased_vs_iMye  = adult_biased_vs_iMye,
  iHSC_biased_vs_iMye   = iHSC_biased_vs_iMye,
  iMye_biased_vs_iHSC   = iMye_biased_vs_iHSC
)

annotated_inventory <- krab_krueppel_genes
annotated_inventory$tags <- sapply(annotated_inventory$gene, function(gene) {
  matches <- names(gene_sets)[sapply(gene_sets, function(vec) gene %in% vec)]
  if (length(matches) == 0) "" else paste(matches, collapse = ", ")
})

write.csv(annotated_inventory,
          "C:/Users/MLuis/OneDrive/Desktop/UNI/Year 3/INTERNSHIP/Inventory datasets/annotated_KZnF_inventory.csv")


#################################################################################
#################################################################################
####  COMPARE iMye/iPSC and iHSC/iPSC  ####

## venn of  iMye-biaed and iHSC-biased vs iPSC

png("venn_biased_vs_iPSC.png", width = 6000, height = 5000, res = 600, bg = "transparent")

grid.newpage()

vp <- viewport(width = 0.5, height = 0.6)

pushViewport(vp)
grid.draw(draw.pairwise.venn(
  area1 = length(iHSC_biased_vs_iPSC),
  area2 = length(iMye_biased_vs_iPSC),
  cross.area = length(intersect(iHSC_biased_vs_iPSC, iMye_biased_vs_iPSC)),
  category = c("iMye-biased", "iHSC-biased"),
  fill = c("#005086", "#007F5F"),
  lty = "blank",
  cex = 2,
  cat.cex = 3,
  cat.pos = c(-50, 50),
  cat.dist = c(-0.5, -0.5),
  ind = FALSE
))
popViewport()

grid.text(paste(setdiff(iHSC_biased_vs_iPSC, iMye_biased_vs_iPSC), collapse = "\n"), x = 0.67, y = 0.43, 
          gp = gpar(cex = 0.7))
grid.text(paste(setdiff(iMye_biased_vs_iPSC, iHSC_biased_vs_iPSC), collapse = "\n"), x = 0.35, y = 0.40, 
          gp = gpar(cex = 0.7))
grid.text(paste(intersect(iMye_biased_vs_iPSC, iHSC_biased_vs_iPSC), collapse = "\n"), x = 0.56, y = 0.47, 
          gp = gpar(cex = 0.7))
grid.text("iHSC-biased and iMye-biased\n KZnFs compared to iPSC", y = unit(0.9, "npc"),
          gp = gpar(fontsize = 30, fontface = "bold"))

dev.off()



## Venn of iPSC-biased vs iHSC or vs iMye

# Open high-res PNG device
png("venn_iPSC_biased.png", width = 6000, height = 5000, res = 600, bg = "transparent")

# Draw Venn diagram
grid.newpage()

vp <- viewport(width = 0.5, height = 0.6)
pushViewport(vp)

grid.draw(draw.pairwise.venn(
  area1 = length(iPSC_biased_vs_iHSC),
  area2 = length(iPSC_biased_vs_iMye),
  cross.area = length(intersect(iPSC_biased_vs_iHSC, iPSC_biased_vs_iMye)),
  category = c("iMye/iPSC", "iHSC/iPSC"),
  fill = c("#D98800", "#F0B94E"),
  lty = "blank",
  cex = 2,
  cat.cex = 3,
  cat.pos = c(-40, 40),
  cat.dist = c(-0.48, -0.48),
  ind = FALSE
))

popViewport()

# Add title
grid.text("iPSC-biased KZnFs in\n iHSC/iPSC and iMye/iPSC", y = unit(0.9, "npc"),
          gp = gpar(fontsize = 30, fontface = "bold"))

# Close device
dev.off()



## Venn of significantly unchanged genes

unchanged_iPSC_iHSC <- KZnF_iHSCbulk_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > -1, log2FoldChange < 1) %>%
  pull(gene)


unchanged_iPSC__iMye <- KZnF_iMye_iPSC %>%
  rownames_to_column(var = "gene") %>%
  filter(padj < 0.01, log2FoldChange > -1, log2FoldChange < 1) %>%
  pull(gene)


png("venn_invariant.png", width = 6000, height = 5000, res = 600, bg = "transparent")

grid.newpage()

vp <- viewport(width = 0.5, height = 0.6)

pushViewport(vp)
grid.draw(draw.pairwise.venn(
  area1 = length(unchanged_iPSC_iHSC),
  area2 = length(unchanged_iPSC__iMye),
  cross.area = length(intersect(unchanged_iPSC_iHSC, unchanged_iPSC__iMye)),
  category = c("iMye/iPSC", "iHSC/iPSC"),
  fill = c("lightblue", "lightblue4"),
  lty = "blank",
  cex = 2,
  cat.cex = 3,
  cat.pos = c(-20, 20),
  cat.dist = c(0.06, 0.06),
  ind = FALSE
))
popViewport()

grid.text("Invariant KZnFs in\n iHSC/iPSC and iMye/iPSC", y = unit(0.9, "npc"),
          gp = gpar(fontsize = 30, fontface = "bold"))

dev.off()

draw.venn(unchanged_iPSC__iMye, unchanged_iPSC_iHSC, NULL)
