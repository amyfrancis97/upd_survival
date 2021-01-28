# function for extracting specific patient information from clinical data
extract_clinical_info <- function(event_interest){
  event <- c()
  for(file in reading_data_output$filenames_clinical_ID){
    clinical_dataframe <- xmlToDataFrame(file, collectNames = TRUE)
    event <- c(event, clinical_dataframe[, event_interest])
    event <- na.omit(event)[1:length(reading_data_output$filenames_clinical_ID)]
  }
  return(event)
}

# creates function
survival_dataframe <- function(filenames_clinical_ID, filenames_aUPD_results, dataset){

  # extracts all UPD data from filenames_aUPD_results
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_", dataset, "_filtered", sep = ""))
  upd_only <- NULL 
  for (file in reading_data_output$filenames_aUPD_results){
    aUPD <- read.table(file, header = TRUE)
    aUPD <- aUPD[aUPD$flagsign == "upd", ]
    upd_only[[file]] <- aUPD
}
  
  setwd("/Users/uw20204/Documents/mini_project_1_UPD/data/clinical")
  # uses previous extract_clinical_info function to extract the following information from patient files
  days_last_fu <- extract_clinical_info("days_to_last_followup")
  days_death <- extract_clinical_info("days_to_death")
  age <- extract_clinical_info("age_at_initial_pathologic_diagnosis")
  sex <- extract_clinical_info("gender")
  stage_event <- extract_clinical_info("stage_event")
  
  # transfer missing data from days until last fu from days until death
  for(indicie in which(days_last_fu == "")){
    days_last_fu[indicie] <- days_death[indicie]
}
  days_last_fu <- as.numeric(days_last_fu)
  
  # extracts censoring info
  # vital status = Alive: censored (1)
  # vital status = dead: not censored (2)
  censoring <- extract_clinical_info("vital_status")
  
  # new data frame for survival analysis 
  aUPD_surv_analysis <- data.frame(days_last_fu, censoring)
  aUPD_surv_analysis[aUPD_surv_analysis$censoring == "Alive", ]$censoring <- 1
  aUPD_surv_analysis[aUPD_surv_analysis$censoring == "Dead", ]$censoring <- 2
  aUPD_surv_analysis$censoring <- as.numeric(aUPD_surv_analysis$censoring)
  
  # change to age to groups <= 50 or >50 years as inn Tuna et al (2019)
  aUPD_surv_analysis$age <- age
  aUPD_surv_analysis[aUPD_surv_analysis$age <= 50, ]$age <- "<= 50"
  aUPD_surv_analysis[aUPD_surv_analysis$age > 50, ]$age <- "> 50"
  aUPD_surv_analysis$sex <- sex
  
  # extracts number of aUPD events for each patient
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_", dataset, "_filtered", sep = ""))
  no_upd_events <- c()
  for(i in reading_data_output$filenames_aUPD_results){
    aUPD_res <- read.table(i, header = TRUE)
    nrow_upd_events <- nrow(aUPD_res[aUPD_res$flagsign == "upd", ])
    no_upd_events <- c(no_upd_events, nrow_upd_events)
  }
  
  aUPD_surv_analysis$no_upd_events <- no_upd_events
  aUPD_surv_analysis$stage_event <- stage_event
  
  # renames stages to make more consistent, this would have to be altered for different cancer types but is not mandatory
  aUPD_surv_analysis[grepl("IV", aUPD_surv_analysis[["stage_event"]]), ]$stage_event <- "stage 4"
  aUPD_surv_analysis[grepl("III", aUPD_surv_analysis[["stage_event"]]) | grepl("7th", aUPD_surv_analysis[["stage_event"]]) 
                     | grepl("6th", aUPD_surv_analysis[["stage_event"]]), ]$stage_event <- "stage 3"
  aUPD_surv_analysis[grepl("II", aUPD_surv_analysis[["stage_event"]]), ]$stage_event <- "stage 2"
  aUPD_surv_analysis[grepl("I", aUPD_surv_analysis[["stage_event"]]), ]$stage_event <- "stage 1"
  aUPD_surv_analysis[grepl("T3", aUPD_surv_analysis[["stage_event"]]), ]$stage_event <- "stage 3"
  
  
  aUPD_surv_analysis <- cbind(as.data.frame(reading_data_output$ID_aUPD), aUPD_surv_analysis)
  colnames(aUPD_surv_analysis) <- c("patient_id", "days_last_fu", "censoring",
                                    "age", "sex", "no_upd_events", "stage_event")
  return(aUPD_surv_analysis)
}
