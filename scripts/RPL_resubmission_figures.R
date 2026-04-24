# Final Figure Plots/Analysis for RPL each vs. Controls Analysis for resubmission

library(targets)
library(tarchetypes)
library(tidyverse)
library(knitr)
library(cowplot)
library(factoextra)
library(FactoMineR)
library(PCAtools)
library(ggsci)
library(ggpubr)
library(DESeq2)
library(RColorBrewer)
library(UpSetR)
library(ComplexUpset)
library(msigdbr)
library(limma)

while (dev.cur() > 1) dev.off()

# PCA Plot ----------------------------------------------------------------

#PCAtools takes in the vst object from DESeq2. Load that from the targets pipeline.
targets::tar_load(vst_norm_object)

#we want the 3rd object in the list, which is ~ sex + condition, with each RPL separately tested against combined CTRLs
vst_obj <- assay(vst_norm_object[[3]])

#let's swap out with symbols
ids <- rownames(vst_obj)
grch38 <- annotables::grch38

grch38 <- grch38 |> 
  dplyr::mutate(symbol = dplyr::if_else(symbol == "", ensgene, symbol)) |> 
  dplyr::distinct(ensgene, .keep_all = TRUE)

swap <-  data.frame(ensgene = ids) |> 
  dplyr::left_join(grch38) |> 
  dplyr::mutate(symbol = make.unique(symbol))

rownames(vst_obj) <- swap$symbol

p <- pca(vst_obj, metadata = colData(vst_norm_object[[3]]), removeVar = 0.9)

p <- biplot(p, colby = "cell_line", 
            colkey = c('CT27' = 'grey75', 'CT29' = 'grey30', 'R002' = 'red', 'R003' = 'green' , 'R004' = 'blue', 'R005' = 'magenta'), 
            labSize = 5, pointSize = 5)

cowplot::ggsave2(plot = p, "./results/plots/PCA_plot_resubmission.pdf", height = 6, width = 6, units = "in")

# Distance Matrix ---------------------------------------------------------

#Get sample distances
sampleDists <- dist(t(vst_obj))
sampleDistMatrix <- as.matrix(sampleDists)

#Set color palette
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#Set annotations
coldata <- as.data.frame(colData(vst_norm_object[[3]])) |> 
  dplyr::select(rpl, sex)

#set annotation colors for each column separately
ann_colors <- list(
  sex = c(Male = "navy", Female = "salmon"),
  rpl = c(Control = "#E69F00", RPL = "#009E73")
)

while (dev.cur() > 1) dev.off()
pdf("./results/plots/Distance_matrix_RPL_Control_Sex.pdf", height=6, width = 7)
ph <- pheatmap::pheatmap(sampleDistMatrix, 
                   clustering_distance_rows = sampleDists, 
                   clustering_distance_cols = sampleDists, 
                   col = colors, 
                   annotation_col = coldata,
                   annotation_row = coldata,
                   annotation_colors = ann_colors
)
print(ph)
dev.off()

# UpSet Plots DEGs --------------------------------------------------------

#read in the DE results object
tar_load(results_list)
de_results_list <- results_list[[3]]
names(de_results_list) <- c("RPL002_vs_Control", "RPL003_vs_Control", "RPL004_vs_Control", "R005_vs_Control")

#get up significant hits
up_hits <- lapply(de_results_list, function(x) x |> dplyr::filter(padj <= 0.05, log2FoldChange > 0) |> tibble::rownames_to_column(var = "ensgene") |> dplyr::pull(ensgene))

down_hits <- lapply(de_results_list, function(x) x |> dplyr::filter(padj <= 0.05, log2FoldChange < 0) |> tibble::rownames_to_column(var = "ensgene") |> dplyr::pull(ensgene))

up_df <- UpSetR::fromList(up_hits) 
down_df <- UpSetR::fromList(down_hits)

up_plot <- ComplexUpset::upset(
  up_df,
  intersect = names(up_hits),    # <-- REQUIRED by ComplexUpset
  set_sizes = FALSE,  # nice default bars on the left
  n_intersections = 8, 
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
) + 
  theme_cowplot() +
  ggtitle("Overlap of Up-regulated DEGs", subtitle = "Adj. p-val <= 0.05")

down_plot <- ComplexUpset::upset(
  down_df,
  intersect = names(down_hits),    # <-- REQUIRED by ComplexUpset
  set_sizes = FALSE,  # nice default bars on the left
  n_intersections = 8, 
  base_annotations = list(
    '# Features' = intersection_size()
  )
) + 
  theme_cowplot() +
  ggtitle("Overlap of Down-regulated DEGs", subtitle = "Adj. p-val <= 0.05")

cowplot::plot_grid(up_plot, down_plot)
ggsave2("./results/plots/UpSet_plots_resubmission.pdf", height = 6, width = 12, units = "in")

# GSEA with limma::camera() -----------------------------------------------

#get the quickomics DESeq2 test results
de_data <- data.table::fread("./results/Quickomics/sex_condition/quickomics_tests.csv") |> 
  dplyr::mutate(test = stringr::str_replace(test, "_vs_", " vs "),
                test = as.factor(test))

#split by test
de_list <- de_data |>
  dplyr::rename(ensembl_gene_id = UniqueID,
                log2FoldChange = logFC) |> 
  dplyr::select(ensembl_gene_id, log2FoldChange, test) |> 
  jmvtools::named_group_split(test)

#convert split into named vectors, values = logFC, names = ENSGENE
de_list_vectors <- lapply(de_list, FUN = function(x) { x |> dplyr::select(-test) |> tibble::deframe() } )

#get the geneset information from msigdbr

#pull the "C2" group which includes WikiPathways, KEGG, and others.
c2_genesets <- msigdbr(db_species = "HS", collection = "C2")

#subset the C2 gene sets to only keep gene sets that include senescence terms.
senescence_gs <- c2_genesets |> dplyr::filter(stringr::str_detect(gs_name, "SENESCENCE"))

#Pull the Hallmark 50 gene sets
hallmark_gs <- msigdbr(db_species = "HS", collection = "H")

#Now, combine Hallmark with Senescence Pathways
comb_gs <- dplyr::bind_rows(hallmark_gs, senescence_gs)

#HM + Senescence
gs_list <- split(x = comb_gs$ensembl_gene, f = comb_gs$gs_name)

run_camera <- function(gs_list, stat_vector) {
  
  #make index
  cam_index <- limma::ids2indices(gs_list, names(stat_vector))
  
  #run cameraPR
  cam_res <- limma::cameraPR(statistic = stat_vector, index = cam_index, use.ranks = TRUE)
}

#Run camera on Hallmark + 23 Senescence Pathways
camera_res_full_sen <- lapply(de_list_vectors, run_camera, gs_list = gs_list)
camera_res_full_sen <- lapply(camera_res_full_sen, function(x) x |> tibble::rownames_to_column(var = "GeneSet"))

#Combine into one output data frame
camera_res_hm_comb_sen <- camera_res_full_sen |> 
  dplyr::bind_rows(.id = "Contrast") |> 
  dplyr::filter(FDR <= 0.01)

write.csv(camera_res_hm_comb_sen, "./results/Combined_Hallmark_plus_senescence_GSEA_with_camera_significant_FDR_0pt01_sex_condition.csv", row.names = F)

#Make plots
camera_plot <- function(x, n_path = 10, colors_use = NULL, x_axis_min = 0, x_axis_max = 16) {
  
  #Check that Contrast is a factor
  if(is.factor(x$Contrast) == FALSE) {
    print("Contrast column in input data frame not passed as a factor. Coercing to factor with default level names. Change levels in input and re-run if different levels are desired.")
    x$Contrast <- factor(x$Contrast, levels = unique(x$Contrast))
  }
  
  #Specifying the color mapping to levels in the contrast column (which is a factor)
  contrast_lvls <- levels(x$Contrast)
  print(paste0("Contrast Levels: ", contrast_lvls))
  n_contrasts<- length(contrast_lvls)
  
  #automatically generate colors for contrasts using ggsci npg (less than 10 contrasts) or ggsci igv (10+ contrasts) if no user specific colors are provided.
  if(is.null(colors_use)) {
    if(n_contrasts <= 10) {
      #npg palette handles up to 10 colors
      contrast_cols <- ggsci::pal_npg("nrc")(n_contrasts)
    } else {
      #igv palette handles up to 51 colors
      contrast_cols <- ggsci::pal_igv("default")(n_contrasts)
    }
  }
  
  #Else, if user provides colors, then use those. Make sure that the number of colors provided match the number of levels in the contrasts.
  #If the numbers don't match, either subset or pad with grey and alert user.
  if(!is.null(colors_use)) {
    
    n_colors_input <- length(colors_use)
    
    #check that length of vector of colors provided matches the number of contrast levels, if not stop.
    if(n_contrasts < n_colors_input) {
      message("Too many colors provided with colors_use for number of contrast levels, using subset of input colors.")
      contrast_cols <- colors_use[n_contrasts]
    } else if (n_contrasts > n_colors_input) {
      message("Not enough colors provided with colors_use for number of contrast levels, padding with grey80. Fix input colors as needed.")
      contrast_cols <- c(colors_use, rep("grey60", n_contrasts - n_colors_input))
    } else {
      contrast_cols <- colors_use
    }
  }
  
  #link the colorst to contrast factor levels for consistency across panels.
  names(contrast_cols) <- contrast_lvls
  
  #get the range of values for NGenes for Up and Down
  ext_fun    <- scales::breaks_pretty(n = 4)
  size_breaks <- ext_fun(c(min(x$NGenes), max(x$NGenes)))
  
  #modify the pathway labels
  x <- x |> 
    dplyr::mutate(GeneSet = stringr::str_replace_all(GeneSet, "_", " ")) |> 
    #GeneSet = stringr::str_to_sentence(GeneSet))
    dplyr::group_by(Direction, GeneSet) |> 
    dplyr::mutate(MaxSig = max(-log10(FDR)))
  
  #Get top 10 pathways Up and Down
  top_gs <- x |> dplyr::ungroup() |> 
    dplyr::group_by(Direction, GeneSet) |> 
    dplyr::summarise(MaxSigFilt = max(MaxSig)) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(Direction) |> 
    dplyr::arrange(desc(MaxSigFilt)) |> 
    dplyr::slice_head(n=n_path) |> 
    jmvtools::named_group_split(Direction)
  #want dot plot with x = FDR, size = NGenes
  up <- x |> 
    dplyr::filter(Direction == "Up", GeneSet %in% top_gs$Up$GeneSet) |> 
    ggplot(aes(x=-log10(FDR), y = reorder(GeneSet, MaxSig), size = NGenes, fill = Contrast)) +
    geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_point(alpha = 0.7, color = "black", pch = 21) +
    scale_size_area(
      name   = "# Genes in Set",
      breaks = size_breaks,
      limits = c(min(size_breaks), max(size_breaks)),
      max_size = 10,
      oob = scales::oob_squish
    ) +
    ggtitle("Enriched Pathways, Up Regulated DEGs") +
    xlim(x_axis_min,x_axis_max) +
    ylab("") +
    scale_fill_manual(values = contrast_cols) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = 12)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 6), order = 1))
  
  down <- x |> 
    dplyr::filter(Direction == "Down", GeneSet %in% top_gs$Down$GeneSet) |> 
    ggplot(aes(x=-log10(FDR), y = reorder(GeneSet, MaxSig), size = NGenes, fill = Contrast)) +
    geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_point(alpha = 0.7, color = "black", pch = 21) +
    scale_size_area(
      name   = "# Genes in Set",
      breaks = size_breaks,
      limits = c(min(size_breaks), max(size_breaks)),
      max_size = 10,
      oob = scales::oob_squish
    ) +
    ggtitle("Enriched Pathways, Down Regulated DEGs") +
    xlim(x_axis_min,x_axis_max) +
    ylab("") +
    scale_fill_manual(values = contrast_cols) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black", size = 12)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 6), order = 1))
  
  #combine plots into one output plot
  out <- cowplot::plot_grid(up, down, ncol = 1, align = "v", axis = "lr")
  
  return(out)
  
}

camera_plots_sen <- camera_plot(camera_res_hm_comb_sen, n_path = 15)
camera_plots_sen
cowplot::ggsave2("./results/plots/Camera_pathway_dotplots_sex_condition_with_senescence.pdf", height = 10, width = 12, units = "in")

#try matching colors
camera_plots_sen2 <- camera_plot(camera_res_hm_comb_sen, n_path = 15, colors_use = c("red", "green", "blue", "magenta"))
camera_plots_sen2
cowplot::ggsave2("./results/plots/Camera_pathway_dotplots_sex_condition_with_senescence_v2.pdf", height = 10, width = 12, units = "in")


# Annotated DE Results ----------------------------------------------------

tar_load(annotated_significant_results)
de_significant <- annotated_significant_results[[3]]
names(de_significant) <- c("RPL002.vs.Control", "RPL003.vs.Control", "RPL004.vs.Control", "RPL005.vs.Control")
openxlsx::write.xlsx(de_significant, file = "./results/Significant_DEGs_sex_vs_condition.xlsx")

#Add columns in the gene set list to add if it was significant for each contrast
comb_geneset_with_de <- comb_gs |> 
  dplyr::mutate(R002 = ensembl_gene %in% de_significant$RPL002.vs.Control$ensgene,
                R003 = ensembl_gene %in% de_significant$RPL003.vs.Control$ensgene,
                R004 = ensembl_gene %in% de_significant$RPL004.vs.Control$ensgene,
                R005 = ensembl_gene %in% de_significant$RPL005.vs.Control$ensgene)

write.csv(comb_geneset_with_de, "./results/Geneset_genes_reference_with_DE_annotations.csv", row.names=F)

# Individual Gene Bar Plots -----------------------------------------------

#Per discussion at bottom here from Mike Love - use vst counts to visualize
#https://support.bioconductor.org/p/112214/

# Get data for ggplot2
genes_to_plot <- c("CDKN1A", "CDKN2A", "CDKN2B", "TGFB1", "TGFB2", "SERPINE1", "CXCL2", "IGFBP7", "IL1A", "INHBA", "IL6", "JUN")
ensemble_ids_to_plot <- annotables::grch38 |> 
  dplyr::filter(symbol %in% genes_to_plot) |> 
  dplyr::pull(ensgene)

# pull the vst normalized counts matrix to use for plotting.
vst_norm_mat <- assay(vst_norm_object[[3]])

#subset to get only the genes we want
plot_df <- vst_norm_mat[ensemble_ids_to_plot, ] |> as.data.frame() |> 
  tibble::rownames_to_column(var = "ensgene") |> 
  dplyr::left_join(select(annotables::grch38, ensgene, symbol)) |> 
  dplyr::select(-ensgene) |> 
  dplyr::select(symbol, everything()) |> 
  tidyr::pivot_longer(cols = -symbol, names_to = "Sample", values_to = "expression") |> 
  tidyr::separate_wider_delim(cols = "Sample", names = c("Line", "Replicate"), delim = "-") |> 
  dplyr::mutate(Group = dplyr::if_else(stringr::str_detect(Line, "CT"), "Control", Line))

#Set colors to match other figures
p <- ggplot(plot_df, aes(x=Group, y=expression, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c('Control' = 'grey75', 'R002' = 'red', 'R003' = 'green' , 'R004' = 'blue', 'R005' = 'magenta')) +
  facet_wrap(~symbol, scale = "free_y") +
  theme_bw(base_size = 14) +
  xlab("") +
  ylab("Normalized Counts") +
  theme(
    legend.position = "none",
    text = element_text(color = "black"),          
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 16)
  )

cowplot::ggsave2(plot = p, "./results/plots/Gene_boxplots.pdf", height = 9, width = 11, units = "in")
