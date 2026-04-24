# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline


options(bitmapType = "cairo")
knitr::opts_chunk$set(dev = "png", dev.args = list(type = "cairo"))
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.


# Set target options:
tar_option_set(
  packages = c("tibble", "DESeq2", "tidyverse", "cowplot", "ggsci", "ggpubr", "jmvtools", "RColorBrewer", "UpSetR", "data.table", "quarto", "WebGestaltR", "DEGreport", "SPIA") # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("./scripts/targets_pipeline_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_file_read(
    name = sample_sheet_tracked,
    command = "./data/trophoblast_sample_sheet.csv",
    read = data.table::fread(!!.x),
    description = "Tracks samplesheet to monitor for changes, will invalidate if updated, also reads the file in."
  ),
  tar_file_read(
    name = gene_counts_tracked,
    command = "../Varberg.trophoblast.cell.culture.bulkRNAseq/results/nfcore_rnaseq/nfcore_output/star_rsem/rsem.merged.gene_counts.tsv",
    read = data.table::fread(!!.x),
    description = "Tracks gene counts matrix with results to monitor for changes, also reads in the matrix as a target."
  ),
  tar_target(
    counts_matrix,
    A01_make_counts(gene_counts_tracked),
    description = "Counts matrix formatted to use to create DESeq2 object. Removes all rows that have zero counts in all samples."
  ),
  tar_target(
    name = column_data,
    command = A02_make_colData(sample_sheet_tracked, counts_matrix),
    description = "Modified sample sheet/metadata to use for colData to create DESeq2 object."
  ),
  tar_target(
    name = design_lrt,
    command = A03_LRT_model_check(counts_matrix, column_data, full_model = "~ sex + rpl", reduced_model = "~ rpl"),
    description = "LRT test to see if the model is better with including sex as a covariate."
  ),
  tar_target(
    name = design_setup,
    command = list(
      list(
        design = "condition",
        label = "cond",
        comparisons = list(
          c("condition", "R002", "Control"),
          c("condition", "R003", "Control"),
          c("condition", "R004", "Control"),
          c("condition", "R005", "Control")
        )
      ),
      list(
        design = "rpl",
        label = "rpl",
        comparisons = list(
          c("rpl", "RPL", "Control")
        )
      ),
      list(
        design = "sex + condition",
        label = "sex_condition",
        comparisons = list(
          c("condition", "R002", "Control"),
          c("condition", "R003", "Control"),
          c("condition", "R004", "Control"),
          c("condition", "R005", "Control")
        )
      )
    ),
    iteration = "list",
    description = "Use this to create a list of specific pair-wise comparisons/contrasts you want to get results for from your model."
  ),
  # 2) Dynamically branch over designs to build each dds
  tar_target(
    name = dynamic_dds,
    command = {
      
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData   = column_data,
        design    = as.formula(paste0("~ ", design_setup$design))
      )
      
      S4Vectors::metadata(dds)$branch_label <- design_setup$label
      S4Vectors::metadata(dds)$comparisons <- design_setup$comparisons
      dds
    },
    pattern   = map(design_setup),
    iteration = "list"
  ),
  # 3) Downstream targets can now also branch automatically
  tar_target(
    name = dynamic_dds_fit,
    command = {
      dds_fit <- DESeq2::DESeq(dynamic_dds, test = "Wald")
      S4Vectors::metadata(dds_fit)$branch_label <- S4Vectors::metadata(dynamic_dds)$label
      S4Vectors::metadata(dds_fit)$comparisons <- S4Vectors::metadata(dynamic_dds)$comparisons
      dds_fit
    },
    pattern   = map(dynamic_dds),
    iteration = "list"
  ),
  tar_target(
    name = deseq_results,
    command = {
      results <- DESeq2::results(dynamic_dds_fit)
      S4Vectors::metadata(results)$branch_label <- S4Vectors::metadata(dynamic_dds_fit)$label
      S4Vectors::metadata(results)$comparisons <- S4Vectors::metadata(dynamic_dds_fit)$comparisons
    },
    pattern = map(dynamic_dds_fit),
    iteration = "list",
    description = "Pulls the results from the modeling performed in the DESeq command."
  ),
  tar_target(
    name = results_list,
    command = {
      comps <- S4Vectors::metadata(dynamic_dds_fit)$comparisons
      lab <- S4Vectors::metadata(dynamic_dds_fit)$branch_label
      
      output <- purrr::map(comps, function(comp) {
        res <- DESeq2::results(dynamic_dds_fit, contrast = comp, test = "Wald")
        lfc <- DESeq2::lfcShrink(dynamic_dds_fit, contrast = comp, res = res, type = "ashr")
        
        #return data frame
        res_df <- as.data.frame(lfc)
        res_df$comparison <- paste(comp, collapse = "_")
        res_df
      })
      
      output
    },
    pattern = map(dynamic_dds_fit),
    iteration = "list",
    description = "Extracts results specified in the contrasts object. Returns a list of data frames. Sets names of list items."
  ),
  tar_target(
    name = annotated_results_list,
    command = lapply(results_list, A04_annotate_results, species = "human", id_type = "ensgene"),
    pattern = map(results_list),
    iteration = "list",
    description = "Uses annotables package reference tables to add gene annotations to results data frames."
  ),
  tar_target(
    name = annotated_significant_results,
    command = lapply(annotated_results_list, function(x) x |> dplyr::filter(padj <= 0.05, abs(log2FoldChange) >= 0.58)),
    pattern = map(annotated_results_list),
    iteration = "list",
    description = "Apply FDR and/or log2FoldChange filtering/thresholds to define significant hits. Modify as desired."
  ),
  tar_target(
    name = vst_norm_object,
    command = vst(dynamic_dds_fit),
    pattern = map(dynamic_dds_fit),
    iteration = "list",
    description = "VST normalization should be used for PCA, heatmaps, boxplots etc., but not for any DE testing. Can use rlog() instead, vst faster for datasets with many (>30) samples."
  ),
  tar_target(
    name = vst_norm_counts_matrix,
    command = assay(vst_norm_object),
    pattern = map(vst_norm_object),
    iteration = "list",
    description = "VST normalized counts in matrix form, to use for plots, save as CSV, etc."
  ),
  tar_target(
    name = quickomics_files,
    command = A05_quickomics_export(normalized_dds_object = vst_norm_object,
                                    de_results_list = results_list,
                                    model_name = design_setup$label,
                                    outDir = "./results/Quickomics"),
    pattern = map(vst_norm_object, results_list, design_setup),
    iteration = "list",
    description = "Creates, saves CSV files in outDir, and returns list object of data.frames needed for Quickomics file upload."
  )
)
