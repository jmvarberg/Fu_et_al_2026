#R Script to create the md5sum hash info for GEO upload

#read in the sample sheet
ss <- data.table::fread("./data/trophoblast_sample_sheet.csv")

#pivot longer to get all fastq files in a column and create basenames
ss_long <- ss |> 
  tidyr::pivot_longer(cols = c("fastq_1", "fastq_2"), names_to = "read", values_to = "path") |> 
  dplyr::mutate(file_name = basename(path),
                md5 = unname(tools::md5sum(path)))

#pivot rows back to wide to match their expected format
ss_wide <- ss_long |> 
  tidyr::pivot_wider(values_from = c("path", "file_name", "md5"), names_from = "read")

write.csv(ss_wide, "./results/tables/GEO_submission_information.csv", row.names = F)


#Check md5sum after copying fastq files locally for GEO upload
geo_copy <- list.files("./data/GEO_raw", full.names = TRUE)
md5_copy <- sapply(geo_copy, tools::md5sum)
