# final master script automated
# this is the master script for survival analysis
# this contains all the scripts for full analysis but many are not part of the main pipeline
# run code line by line

# required libraries
load_packages <- c("BiocManager","stringi", "XML", "stringr", "methods", "survival", "survminer",
                   "cowplot", "dplyr", "ggpubr", "GenomicFeatures", "Homo.sapiens",
                   "TxDb.Hsapiens.UCSC.hg19.knownGene", "rstatix", "ggplot2", "biomaRt",
                   "karyoploteR", "BSgenome.Hsapiens.UCSC.hg19", "regioneR", "GenomicRanges",
                   "ChIPseeker", "clusterProfiler", "plyr", "stringr", "rentrez",
                   "tidyr", "IRanges", "factoextra", "NbClust", "lattice", "wesanderson", "ggsci",
                   "gridExtra", "VariantAnnotation", "COSMIC.67")
lapply(load_packages, require, character.only = TRUE)

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_ordered")
#if this script doesn't run its fine move on- just means there is no missing clinical data
#source("1_reading_aUPD_data_final_script.R") # checks to see if there are an missing clinical data and if so deletes the aUPD file as cannot be analysed
#reading_data_output <- dataset_details("COAD") # runs for specified dataset
#complete_analysis <- function(dataset){
setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("1.2_reading_aUPD_data_final_script_missing_data_removed.R") # extracts dataset details from clinical and aUPD data
reading_data_output <- dataset_details("READ") # runs for specified dataset
reading_data_output
setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("2_generate_survival_df_final.R")
aUPD_surv_analysis <- survival_dataframe(reading_data_output$filenames_clinical_ID, reading_data_output$filenames_aUPD_results, "READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("3_generate_upd_location_df_final.R")
#tcn_location_dataframe_READ <- total_copy_number_location(filenames_clinical_ID, filenames_aUPD_results, dataset)
tcn_location_dataframe <- total_copy_number_location(filenames_clinical_ID, filenames_aUPD_results, "READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("3.2_number_upd_tables.R")
number_upd_chrom_number_patients_upd("READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("4.1_survival_by_chrom_final.R")
surv_plot("READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("5_cox_hazards_chrom_level_final.R")
coxph_chrom("READ")


setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("6_aUPD_percentage_length_final.R")
perc_aupd_coverage("READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("9_karyotype_plotting_final.R")
karyoplot_distribution_upd("READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("12_sliding_windows_survival_final.R")
regionaldf_list_patient_id <- sliding_window_analysis("READ")
regionaldf_list_patient_id
setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("16_segmental_vs_whole_chrom_surv_final.R")
segmental_whole_chrom_surv("READ", "Segmental", tcn_location_dataframe_segmental)
segmental_whole_chrom_surv("READ", "Whole_chrom", tcn_location_dataframe_whole_chrom)
#}

#complete_analysis("READ")

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("15_pulling_COSMIC_mutations_inc_whole_chrom.R")
# a = chrom_number; b = start_window_number; c = end_window_number; d = region_start; e = region_end)
dataset_genes_cosmic <- pull_cosmic_genes(as.numeric(gsub("chr", "", (mydata$chrom)))[1],
                                          rownames(mydata)[1], rownames(mydata)[length(rownames(mydata))],
                                          mydata$start[1], mydata$end[length(mydata$end)])
dataset_genes_cosmic
setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/scripts/final_tcn_all_v2")
source("15.2_saving_positive_negative_COSMIC.R")
saving_COSMIC(dataset_genes_cosmic, "READ", paste(mydata$start[1], mydata$end[length(mydata$end)], sep ="_"))
positive_patient_mutations
oncogenes$Genome.Location <- paste("chr", oncogenes$Genome.Location, sep = "")
oncogene_df <- NULL
for(i in 1:length(oncogenes$Genome.Location)){
if(nchar(oncogenes$Genome.Location[i]) > 7){
  oncogene_vector <- rbind(oncogene_vector, data.frame(gene = oncogenes$Gene.Symbol[i], location= oncogenes$Genome.Location[i]))
}
}

