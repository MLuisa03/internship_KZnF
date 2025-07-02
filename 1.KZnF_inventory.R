library(biomaRt)
library(dplyr)
library(VennDiagram)
library(grid)
library(BioVenn)



#Get Ensembl Biomart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Get genes annotated with GO:0003677 (DNA-binding proteins)
dna_binding_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
  filters = "go",
  values = "GO:0003677",
  mart = ensembl
)

#Only keep GO:0003677 
dna_binding_genes <- dna_binding_genes %>% filter(go_id == "GO:0003677")

#Remove genes without names 
dna_binding_genes <- dna_binding_genes %>% filter(external_gene_name != "")

#Remove duplicates
dna_binding_genes <- dna_binding_genes %>% distinct(ensembl_gene_id, .keep_all = TRUE)


#Get genes annotated with IPR036236 (C2H2 Zinc Finger domain)
c2h2_znfs <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "interpro"),
  filters = "interpro",
  values = "IPR036236",  # C2H2 Zn fingers
  mart = ensembl
)

#Filter the original list of DNA-binding genes to keep only those with C2H2 ZnFs
c2h2_genes <- dna_binding_genes[dna_binding_genes$ensembl_gene_id %in% c2h2_znfs$ensembl_gene_id, ]


#Get genes annotated with all family IDs related to Krueppel-family ZnFs
krueppel_znfs <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "interpro"),
  filters = "interpro",
  values = c("IPR050169", "IPR050527", "IPR050826", "IPR051967", "IPR052274"),  # Krueppel C2H2 Zn fingers
  mart = ensembl
)

#Filter the original list of C2H2 ZnF genes to keep only those in Krueppel-family ZnFs
krueppel_genes <- c2h2_genes[c2h2_genes$ensembl_gene_id %in% krueppel_znfs$ensembl_gene_id, ]


# Get genes annotated with IPR001909 (KRAB domain)
krab_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "interpro"),
  filters = "interpro",
  values = "IPR001909",  # KRAB domain
  mart = ensembl
)

#Filter the original list of DNA-binding genes to keep only those that also have KRAB domains
krab_krueppel_genes <- krueppel_genes[krueppel_genes$ensembl_gene_id %in% krab_genes$ensembl_gene_id, ]


krab_krueppel_genes <- krab_krueppel_genes %>% distinct(external_gene_name, .keep_all = TRUE)



#Load the Uniprot query: human + reviewd + evidence at protein level + DNA-binding + krueppel C2H2-type zinc-finger protein family + KRAB

UniProt_KKZnF <- read.delim("downloaded_from_UniProt")


#Load RefSeq query: KRAB + ZnF

RefSeq <- read.delim("downloaded_from_RefSeq")



##  Overlap KKZnF from ensembl, UniProt and RefSeq  ##

grid.newpage()
pushViewport(viewport(width = 0.7, height = 0.7))

grid.draw(draw.triple.venn(
  area1 = length(unique(UniProt_KKZnF$Gene.Names..primary.)),
  area2 = length(unique(krab_krueppel_genes$external_gene_name)),
  area3 = length(unique(RefSeq$Symbol)),
  n12 = length(intersect(UniProt_KKZnF$Gene.Names..primary., krab_krueppel_genes$external_gene_name)),
  n13 = length(intersect(UniProt_KKZnF$Gene.Names..primary., RefSeq$Symbol)),
  n23 = length(intersect(krab_krueppel_genes$external_gene_name, RefSeq$Symbol)),
  n123 = length(Reduce(intersect, list(UniProt_KKZnF$Gene.Names..primary., krab_krueppel_genes$external_gene_name, RefSeq$Symbol))),
  category = c("UniProt", "Ensembl", "RefSeq"),
  fill = c("#1f77b4", "#ff7f0f", "#5ca02c"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.3,
  cat.col = c("#1f77b4", "#ff7f0f", "#5ca02c"),
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  ind = FALSE
))

popViewport()

# Title
grid.text("KRAB Krueppel ZnFs in RefSeq, UniProt, Ensembl", y = unit(0.93, "npc"),
          gp = gpar(fontsize = 16, fontface = "bold"))


##  Same but in the BioVenn version (still useful to get overlap lists quickly)

venn_RS_UP_en <- draw.venn(
  list_x = UniProt_KKZnF$Gene.Names..primary.,
  list_y = krab_krueppel_genes$x,
  list_z = RefSeq$Symbol,
  xtitle = "Uniprot",
  ytitle = "ensembl",
  ztitle = "RefSeq",
  title = "KKZnFs RefSeq + UniProt + ensembl", 
  subtitle = NULL, 
  xt_s = 1,
  yt_s = 1,
  zt_s = 1,
  nr_s = 3,
  x_c = "#1f77b4",  
  y_c = "#ff7f0f",  
  z_c = "#5ca02c")


#Get intersections

intersection_RS_UP_only <- venn_RS_UP_en$xz_only 

intersection_RS_en_only <- venn_RS_UP_en$yz_only  

intersection_UP_en_only <- venn_RS_UP_en$xy_only 

intersection_RS_UP_en <- venn_RS_UP_en$xyz 



#Files from Trono paper and KRABopedia

Trono_2023 <- read.delim("Supplementary_table_3_de_Tribolet-Hardy_2023")


KRABopedia <- read.delim("list_of_KRAB_genes_on_the_website_of_KRABopedia")



grid.newpage()
pushViewport(viewport(width = 0.7, height = 0.7))

grid.draw(draw.triple.venn(
  area1 = length(unique(intersection_RS_UP_en)),
  area2 = length(unique(Trono_2023$assigned_gene)),
  area3 = length(unique(KRABopedia$gene)),
  n12 = length(intersect(intersection_RS_UP_en, Trono_2023$assigned_gene)),
  n13 = length(intersect(intersection_RS_UP_en, KRABopedia$gene)),
  n23 = length(intersect(Trono_2023$assigned_gene, KRABopedia$gene)),
  n123 = length(Reduce(intersect, list(intersection_RS_UP_en, Trono_2023$assigned_gene, KRABopedia$gene))),
  category = c("Databases", "de T.-H. 2023", "KRABopedia"),
  fill = c("#929E51", "purple", "yellow"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.3,
  cat.col = c("darkgreen", "purple", "brown"),
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.04),
  ind = FALSE
))

popViewport()

# Add title
grid.text("KRAB Krueppel ZnFs from Databases, \n de Tribolet-Hardy (2023) and KRABopedia", 
          y = unit(0.93, "npc"), 
          gp = gpar(fontsize = 15, fontface = "bold"))



#Obtain the final overlap list 

venn_database_Tr_KPedia <- draw.venn(
  list_x = intersection_RS_UP_en,
  list_y = Trono_2023$assigned_gene,
  list_z = KRABopedia$gene,
  xtitle = NULL,
  ytitle = NULL,
  ztitle = NULL,
  title = "KKZnFs from databases (UniProt, ensembl and RefSeq), de Tribolet-Hardy (2023), and KRABopedia", 
  subtitle = NULL, 
  xt_s = 1.5,
  yt_s = 1.5,
  zt_s = 1.5,
  nr_s = 2,
  x_c = "#2ca03b",
  y_c = "purple",
  z_c = "yellow")

text(0.3, 1.0, "Databases", col = "darkgreen", font = 2, cex = 1.2)
text(0.2, 0.5, "de T. - H.", col = "purple", font = 2, cex = 1.2)
text(0.7, 0.5, "KRABopedia", col = "brown", font = 2, cex = 1.2)


intersection_database_Trono_KRABopedia <- venn_database_Tr_KPedia$xyz 
