# now we must pull out the genes from COSMIC
# these contain sequencing files for aUPD- positive and -negative patients at these regions

data(cosmic_67, package = "COSMIC.67")
setwd("/Users/uw20204/Downloads/split_COSMIC_data")
filenames_COSMIC_split <- list.files("/Users/uw20204/Downloads/split_COSMIC_data", 
                                     recursive=TRUE) #filters to xml files

ID_aUPD_COSMIC_format <- paste("TCGA-", reading_data_output$ID_aUPD, "-01", sep = "")
colnames_COSMIC <- colnames(read.csv(file = filenames_COSMIC_split[1], header = TRUE))

COSMIC_dataframe <- NULL
for(i in 1:length(filenames_COSMIC_split)){
  df_split <- read.csv(file = filenames_COSMIC_split[i], header = FALSE)
  if(nrow(df_split %>%
          filter(V5 %in% ID_aUPD_COSMIC_format)) > 0){
    df_split <- df_split %>%
      filter(V5 %in% ID_aUPD_COSMIC_format)
  }
  COSMIC_dataframe <- rbind(COSMIC_dataframe, df_split)
}

colnames(COSMIC_dataframe) <- colnames_COSMIC
COSMIC_dataframe$MUTATION_GENOME_POSITION <- paste("chr", COSMIC_dataframe$MUTATION_GENOME_POSITION, sep = "")
COSMIC_dataframe
# create a function to pull out genes for UPD-positive and negative patients for specified regions
##############################
#this will create a dataframe of GRanges that can be used to merge to original df and then use regions of overlap.
#analysis for chr1
pull_cosmic_genes <- function(chrom_number, start_window_number, end_window_number, region_start, region_end){
COSMIC_dataframe_chr1 <- COSMIC_dataframe[grepl(paste("chr", chrom_number, ":", sep = ""), COSMIC_dataframe$MUTATION_GENOME_POSITION, fixed = TRUE), ]

# changes the format of genomic ranges to be in proper formating of start end and chrom
chr1_location_df <- NULL
for(region_row in 1:length(COSMIC_dataframe_chr1$MUTATION_GENOME_POSITION)){
  splitting_location <- unlist(str_split(COSMIC_dataframe_chr1$MUTATION_GENOME_POSITION[region_row], ":"))
  splitting_location <- unlist(str_split(splitting_location, "-"))
  chr1_location_df <- rbind(chr1_location_df, splitting_location)
}

colnames(chr1_location_df) <- c("chromosome", "start", "end")

# changes the original chr5 cosmic df to include start and end etc sepaate columns
chr1_location_df <- cbind(chr1_location_df, COSMIC_dataframe_chr1[, -which(colnames(COSMIC_dataframe_chr1) 
                                                                           == "MUTATION_GENOME_POSITION")])
# tidy up gene name column in COSMIC df
genes_without_ensembl_id <- NULL
for(gene_name_col in 1:length(chr1_location_df$GENE_NAME)){
  split <- strsplit(chr1_location_df$GENE_NAME, "_")[[gene_name_col]][1]
  genes_without_ensembl_id <- c(genes_without_ensembl_id, split)
}
chr1_location_df$GENE_NAME <- genes_without_ensembl_id

# finding positive patients ids for region of interest
positive_patient_total <- NULL
for(window_number in start_window_number:end_window_number){
  patient_id_pos <- regionaldf_list_patient_id[[2]][[chrom_number]][which(regionaldf_list_patient_id[[2]][[chrom_number]][, (7 + window_number)] == "Positive"), ]$patient_id
  tcn_location_dataframe_whole_chrom_17 <- tcn_location_dataframe_whole_chrom %>%
    filter(Chrom == chrom_number)
  tcn_location_dataframe_whole_chrom_17 <- unique(tcn_location_dataframe_whole_chrom_17$patient_id)
  positive_patient_total <- unique(c(positive_patient_total, c(patient_id_pos, tcn_location_dataframe_whole_chrom_17)))
  }

# changes format of previous region of interest sample name to match format of COSMIC

patient_ID_positive <- paste("TCGA-", positive_patient_total, "-01", sep = "")
patient_ID_negative <- setdiff(aUPD_surv_analysis$patient_id, patient_ID_positive)
patient_ID_negative <- paste("TCGA-", patient_ID_negative, "-01", sep = "")

# changes to numeric to plot regions
chr1_location_df[, 2] <- as.numeric(chr1_location_df[, 2])
chr1_location_df[, 3] <- as.numeric(chr1_location_df[, 3])
patient_ID_positive
# finds pathogenic mutations in COSMIC for those that are positive in previous df
positive_patients <- chr1_location_df %>%
  filter(SAMPLE_NAME %in% patient_ID_positive) %>%
  filter(start > region_start & end < region_end)%>%
  filter(FATHMM_PREDICTION == "PATHOGENIC")

# finds pathogenic mutations in COSMIC for those that are negative in previous df
negative_patients <- chr1_location_df %>%
  filter(SAMPLE_NAME %in% patient_ID_negative) %>%
  filter(start > region_start & end < region_end)%>%
  filter(FATHMM_PREDICTION == "PATHOGENIC")
return(list(positive_patients, negative_patients))
}
