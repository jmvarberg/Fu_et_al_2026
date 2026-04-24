This repository contains the code used to process the bulk RNAseq data for Fu et al., 2026.

# Primary Analysis - nf-core/rnaseq

The nextflow nf-core/rnaseq pipeline was used to analyze the raw sequencing files, and to generate genes x samples counts matrix using STAR and RSEM. The nextflow submission script is located at __scripts/run_nfcore_rnaseq.sbatch__.

# Downstream Analysis - DESeq2

Secondary analysis using the RSEM counts matrix was performed using DESeq2 (v. 1.46.0) in R 4.4.0. This analysis was run using the {targets} package. The structure of the pipeline is outlined in __scripts/_targets.R__, and the functions that are used in this pipeline are in __/scripts/targets_pipeline_functions.R__.

# Figures

The code used to generate the final figures for the manuscript is located in __/scripts/RPL_resubmission_figures.R__. PDFs of the generated plots are located in __results/plots/__.

# Quickomics

The results from the DESeq2 analysis testing using the '~ sex + condition' model (testing each RPL line against the combined CT27/CT29 controls) were extracted and are provided as CVS files as outlined below:

- __Expression Matrix:__ VST-normalized counts were extracted and are provided in the file __results/Quickomics/sex_condition/quickomics_expression.csv__
- __Differential Expression:__ DESeq2 results were extracted and lfcShrink() was applied using the 'ashr' method. The results for each comparison are combined and provided in a single CSV file at __results/Quickomics/sex_conditon/quickomics_tests.csv__
- __Metadata:__ The colData from the DESeq2 object was extracted and is provided in the file __results/Quickomics/sex_condition/quickomics_md.csv__

These files have been formated to allow them to be uploaded for interactive data exploration using the Quickomics Shiny dashboard, which is publicly available here: https://quickomics.bxgenomics.com
