# this script generates another dataframe, tcn_location_dataframe, which shows the 
# genomic locations of all aUPD events, the patient they below to, and the percentage length of the aUPD event

# Extracts filenames of only files containing aUPD
total_copy_number_location <- function(filenames_clinical_ID, filenames_aUPD_results, dataset){
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_",
              dataset, "_filtered", sep = ""))
  filename <- c()
  for(file in reading_data_output$filenames_aUPD_results){
    if ("upd" %in% read.table(file, header = TRUE)$flagsign){
      filename <- c(filename, file)
    }
  }
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
  write.csv(data.frame(number_perc_patients_upd = paste(as.character(length(filename)), 
                                                        round(length(filename)/length(reading_data_output$filenames_aUPD_results)*100, digits = 1), "%"), 
                       total_number_aupd_patient_data = length(reading_data_output$filenames_aUPD_results)), "total_patients_with_aUPD.csv")
  filename
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_",
              dataset, "_filtered", sep = ""))
  # Extracts patient IDs from above filenames
  ID_aUPD <- c()
  for (file in filename){
    ID_aUPD <- c(ID_aUPD, str_sub(file, (reading_data_output$ID_end + 1), (reading_data_output$ID_end + 7)))
  }
  
  ### Identifying gene location of aUPD
  tcn_location <- NULL #number of upd events for each cancer patient
  for (i in 1:length(reading_data_output$filenames_aUPD_results)){
    aUPD <- read.table(reading_data_output$filenames_aUPD_results[[i]], header = TRUE)
    tcn_location_ind <- data.frame(rep(reading_data_output$ID_aUPD[i], nrow(aUPD)), aUPD$chromosome, 
                                   aUPD$tcnStart, aUPD$tcnEnd, aUPD$flagsign)
    tcn_location <- rbind(tcn_location, tcn_location_ind)
  }
  colnames(tcn_location) <- c("patient_id", "Chrom", "Start","End", "Flagsign")
  return(tcn_location)
}

