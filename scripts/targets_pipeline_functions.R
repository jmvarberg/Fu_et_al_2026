#scripts for targets pipeline


# Analysis Functions ------------------------------------------------------


A01_make_counts <- function(gene_counts) {
  
  counts_matrix <- gene_counts %>% 
    dplyr::select(-`transcript_id(s)`) %>% 
    tibble::column_to_rownames(var = "gene_id") %>% 
    dplyr::mutate(across(where(is.numeric), as.integer))
  
  #remove features/rows where there are zero counts in all columns/samples
  counts_matrix <- counts_matrix[rowSums(counts_matrix != 0) > 0, , drop = F]
  
}

A02_make_colData <- function(sample_sheet, counts_matrix) {
  
  # tar_load(sample_sheet)
  # tar_load(counts_matrix)
  #make the colData as needed for grouping and sample information
  metadata <- sample_sheet %>% 
    dplyr::select(sample, cell_line) %>% 
    dplyr::distinct() |> 
    dplyr::filter(sample %in% colnames(counts_matrix)) |> #this removes any samples that were in the sample sheet but didn't make it through the analysis to the counts matrix.
    dplyr::mutate(condition = dplyr::if_else(stringr::str_detect(cell_line, "CT2"), "Control", cell_line),
                  rpl = dplyr::if_else(condition == "Control", "Control", "RPL"),
                  sex = dplyr::if_else(cell_line %in% c("R002", "CT29"), "Male", "Female"),
                  group = dplyr::if_else(cell_line == "R003", "R003", rpl),
                  ) |> 
    dplyr::arrange(match(sample, colnames(counts_matrix))) |> #this makes sure columns are ordered in the same order as columns in counts matrix
    DataFrame() #this last step needed in order to get the DESeq2 object creation to work
}

A03_LRT_model_check <- function(counts_matrix, column_data, full_model = "~ sex + rpl", reduced_model = "~ rpl") {
  
  # tar_load(counts_matrix)
  # tar_load(column_data)

  #Build the DESeq object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData   = column_data,
    design    = as.formula(full_model)
  )
  
  
  # LRT: compare FULL (~ condition + batch) vs REDUCED (~ condition)
  dds_lrt <- DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced_model))
  
  res_lrt <- results(dds_lrt)
  head(res_lrt)
  
  
  f <- function() {
    cat("Genes improved by adding sex as co-variate (padj < 0.05):",
        sum(res_lrt$padj < 0.05, na.rm=TRUE), "\n")
    
    cat("Total number of genes tested: ",
        sum(!is.na(res_lrt$padj)),"\n")
    
    cat("Fraction total genes improved by adding sex as co-variate: ",
        round(sum(res_lrt$padj < 0.05, na.rm=TRUE)/sum(!is.na(res_lrt$padj)), digits = 3),"\n")
    
    cat("Median deviance reduction (LRT stat):",
        median(res_lrt$stat, na.rm=TRUE), "\n")
  }
  
  # Capture as a character vector (one element per printed line)
  txt_vec <- capture.output(f())
  txt_vec
 
  # Model WITHOUT batch (for comparison)
  dds_without <- DESeqDataSetFromMatrix(counts_matrix, column_data, design = as.formula(reduced_model))
  dds_without <- DESeq(dds_without)
  
  # Model WITH batch (Wald)
  dds_with <- dds
  dds_with <- DESeq(dds_with)
  
  res_no  <- results(dds_without, contrast = c("rpl","RPL","Control"))
  res_yes <- results(dds_with, contrast = c("rpl","RPL","Control"))
  
  table(no = res_no$padj < 0.05, yes = res_yes$padj < 0.05)

  
  # 1) Put both results into data.frames and align by common rownames
  res_no_df  <- as.data.frame(res_no)
  res_yes_df <- as.data.frame(res_yes)
  
  common <- intersect(rownames(res_no_df), rownames(res_yes_df))
  res_no_c  <- res_no_df [common, , drop = FALSE]
  res_yes_c <- res_yes_df[common, , drop = FALSE]
  
  # 2) Define lost / gained on the aligned objects (same order, same set)
  alpha <- 0.05  # (or 0.1 if that’s your project threshold; be consistent)
  lost_ids   <- common[(res_no_c$padj < alpha) & !(res_yes_c$padj < alpha)]
  gained_ids <- common[!(res_no_c$padj < alpha) & (res_yes_c$padj < alpha)]
  
  # 3) Summarize numerics *after* removing rows with NA values
  lost_summary <- summary(na.omit(res_yes_c[lost_ids, c("baseMean","log2FoldChange")]))
  gained_summary <- summary(na.omit(res_yes_c[gained_ids, c("baseMean","log2FoldChange")]))
  lost_summary; gained_summary
  
  
  lost_logical   <- (res_no_c$padj < alpha) & !(res_yes_c$padj < alpha)  # sig -> not sig
  gained_logical <- !(res_no_c$padj < alpha) & (res_yes_c$padj < alpha)  # not sig -> sig
  stable_sig     <- (res_no_c$padj < alpha) & (res_yes_c$padj < alpha)
  stable_nonsig  <- !(res_no_c$padj < alpha) & !(res_yes_c$padj < alpha)
  
  status <- rep("stable_nonsig", length(common))
  status[stable_sig]     <- "stable_sig"
  status[lost_logical]   <- "lost"
  status[gained_logical] <- "gained"
  
  # ---- 3) (Optional but recommended) Shrink LFCs for nicer visualization ----
  # The DESeq2 vignette recommends LFC shrinkage (e.g., apeglm) for plotting. [1](https://bioconductor.org/packages//release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
  # Skip this block if you prefer raw LFCs.
  suppressPackageStartupMessages(library(apeglm))
  # Adjust the 'coef' to your contrast name in the WITH-sex model:
  # Check with: resultsNames(dds_yes)
  # Example: "condition_treated_vs_control"
  coef_name <- resultsNames(dds_with)[grep("rpl", resultsNames(dds_with))[1]]
  
  res_yes_shr <- lfcShrink(dds_with, coef = coef_name, type = "apeglm")
  res_yes_shr_df <- as.data.frame(res_yes_shr)[common, , drop = FALSE]
  plotMA(res_yes_shr)
  # ---- 4) Build plotting data (using WITH-sex model’s baseMean + LFC) ----
  # Choose either shrunken or raw for plotting by setting use_shrinkage
  use_shrinkage <- TRUE
  
  plot_df <- if (use_shrinkage) {
    data.frame(
      gene      = common,
      baseMean  = res_yes_shr_df$baseMean,
      log2FC    = res_yes_shr_df$log2FoldChange,
      padj_yes  = res_yes_c$padj,
      status    = factor(status, levels = c("stable_nonsig","stable_sig","lost","gained"))
    )
  } else {
    data.frame(
      gene      = common,
      baseMean  = res_yes_c$baseMean,
      log2FC    = res_yes_c$log2FoldChange,
      padj_yes  = res_yes_c$padj,
      status    = factor(status, levels = c("stable_nonsig","stable_sig","lost","gained"))
    )
  }
  
  #add annotations from annotables
  plot_df <- plot_df |> dplyr::left_join(annotables::grch38, by = c("gene" = "ensgene"))
  
  # Base color mapping: grey (non-sig both), blue (sig both), red (lost), green (gained)
  cols <- c(stable_nonsig = "grey80", stable_sig = "lightblue", lost = "firebrick", gained = "navy")
  
  # ---- 5) MA Plot with ggplot2 ----
  suppressPackageStartupMessages(library(ggplot2))
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean), y = log2FC, color = status)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.4, color = "black") +
    geom_point(data = plot_df |> dplyr::filter(status %in% c("stable_nonsig", "stable_sig")), alpha = 0.5, size = 1) +
    geom_point(data = plot_df |> dplyr::filter(status %in% c("lost", "gained")), alpha = 0.9, size = 1.5) +
    scale_color_manual(values = cols, name = "Status vs. sex covariate") +
    labs(
      title = sprintf("MA Plot (WITH sex; alpha=%.2g)%s",
                      alpha, if (use_shrinkage) " — LFC shrinkage (apeglm)" else ""),
      x = "mean of normalized counts (baseMean, log10)",
      y = "log2 fold change",
      caption = "Color indicates whether adding 'sex' retained, lost, or gained significance"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "right")
  p_ma
  
  # ---- 6) Save high-res if needed ----
  cowplot::ggsave2("./results/plots/MA_RPL_with_sex_highlight_lost_gained.pdf", width = 8, height = 8, units = "in")
  
  # ---- 7) (Optional) Label a few extreme points for context ----
  # Pick top N by |log2FC| among gained/lost to annotate
  suppressPackageStartupMessages(library(dplyr))
  topN <- 15
  to_annotate <- bind_rows(
    plot_df %>% filter(status == "gained") %>% slice_max(order_by = abs(log2FC), n = topN),
    plot_df %>% filter(status == "lost")   %>% slice_max(order_by = abs(log2FC), n = topN)
  )
  
  p_ma_lab <- p_ma +
    ggrepel::geom_label_repel(
      data = to_annotate,
      aes(label = symbol),
      size = 3, max.overlaps = 50, 
      box.padding = 0.3, 
      point.padding = 0.2, 
      show.legend = FALSE
    )
  
  p_ma_lab
  cowplot::ggsave2("./results/plots/MA_RPL_with_sex_highlight_lost_gained_labeled.pdf", width = 8, height = 8, units = "in")
  
  return(txt_vec)
}

#extract all results and perfrom lfcShrink with 'ashr' method. User provides list of coefficients or contrasts to pull. Returns as data.frame object.
# A03_extract_results <- function(comps, dds, type = c("coef", "contrasts")) {
#   
#   if(type == "coef") {
#     res <- results(dds, name = comps, test = "Wald")
#     res_lfc <- lfcShrink(dds, coef = comps, res = res, type = "ashr")
#   }
#   
#   if(type == "contrasts") {
#     res <- results(dds, contrast = comps, test = "Wald")
#     res_lfc <- lfcShrink(dds, contrast = comps, res = res, type = "ashr")
#     
#   }
#   
#   return(as.data.frame(res_lfc))
# }

#annotate results data frames using annotables
A04_annotate_results <- function(results, species = c("human", "mouse"), id_type = c("ensgene", "entrez", "symbol")) {
  
  #get the reference annotation table
  if(species == "human") {
    ref <- annotables::grch38
  } else if(species == "mouse") {
    ref <- annotables::grcm38
  }
  
  out <- results |> 
    tibble::rownames_to_column(var = id_type) |> 
    dplyr::left_join(ref) |> 
    dplyr::distinct(symbol, .keep_all = TRUE)
}

#normalized_dds_object - either rlog() or vst(); de_results_list - full DE results for each contrast/comparison
A05_quickomics_export <- function(normalized_dds_object, de_results_list, model_name = NULL, outDir = "./results/Quickomics") {
  
  # #Troubleshooting sex+group error in quickomics export
  # normalized_dds_object = vst_norm_object$vst_norm_object_32697ca9f0f5662a
  # de_results_list = results_list$results_list_32697ca9f0f5662a
  # model_name = design_setup[[5]]$label 
  # outDir = "./results/Quickomics"

  output_location <- file.path(outDir, model_name)
  
  #Create the Quikomics output directory as specified in outDir path
  if(!dir.exists(output_location)) {
    dir.create(output_location, recursive = T)
    print(paste0("Creating output directory for Quickomics files at: ", output_location))
  }
  
  #metadata: sampleid, group, additional columns. sampleid must match expresion data file column names
  quickomics_md <- as.data.frame(colData(normalized_dds_object)) |> 
    dplyr::rename(sampleid = sample,
                  genotype = cell_line) |> 
    dplyr::select(sampleid, group, everything(), -sizeFactor)
  quickomics_md
  
  write.csv(quickomics_md, file = file.path(output_location,"quickomics_md.csv"), row.names = F)
  
  #expression data - UniqueID, sampleids, use vst/rlog normalized counts object
  quickomics_expression <- as.data.frame(assay(normalized_dds_object)) |> 
    tibble::rownames_to_column(var = "UniqueID") |> 
    dplyr::filter(rowSums(across(where(is.numeric))) != 0)
  
  write.csv(quickomics_expression, file = file.path(output_location, "quickomics_expression.csv"), row.names = F)
  
  #Comparison data: UniqueID, test, Adj.P.Value, P.Value, logFC. test labels must match group MD labels with 'vs' separator no spaces
  quickomics_tests <- de_results_list |> 
    #dplyr::bind_rows(.id = "test") |> 
    dplyr::bind_rows() |> 
    tibble::rownames_to_column(var = "UniqueID") |> 
    dplyr::mutate(UniqueID = stringr::str_remove(UniqueID, pattern = "\\.\\.\\..*$"),
                  parts = stringr::str_split_fixed(comparison, "_", n = 3),
                  test = stringr::str_c(parts[, 2], "vs", parts[, 3], sep = "_")) |> 
    dplyr::filter(UniqueID %in% quickomics_expression$UniqueID) |> 
    dplyr::rename(P.Value = pvalue,
                  Adj.P.Value = padj, 
                  logFC = log2FoldChange) |> 
    dplyr::select(-baseMean, -lfcSE, -parts) |> 
    na.omit()
  
  write.csv(quickomics_tests, file = file.path(output_location, "quickomics_tests.csv"), row.names = F)
  
  return(list(meta = quickomics_md, 
              norm_expr = quickomics_expression,
              de_res = quickomics_tests))
  
}
