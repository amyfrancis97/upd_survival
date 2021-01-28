# this script reads the aUPD files from directory
# this script finds whether there are any missing clinical files for aUPD files
# then need to go back and delete the upd file from the documents
# use script version 1.2 to continue with file deleted
#if this script does not run then it is because there is no missing clinical data
# hence, this step can be skipped

dataset_details <- function(dataset){
# lists the files in each of the directories
filenames_aUPD_results <- list.files(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_", dataset, "_filtered", sep = ""))

# finds end location of "TCGA-" from filenames_aUPD_results string
ID_end <- as.vector(str_locate(filenames_aUPD_results[1], "TCGA-")[, "end"])

# extracts patient ID from filenames_aUPD_results based on location of TCGA- in the filename string
ID_aUPD <- c()
for (file in filenames_aUPD_results){ # ID_end from previous step
  ID_aUPD <- c(ID_aUPD, str_sub(file, (ID_end + 1), (ID_end + 7))) 
  #+1 & +7 is where ID starts and ends
}

setwd("/Users/uw20204/Documents/mini_project_1_UPD/data/clinical")
# extracts file names from clinical data
filenames_clinical <- list.files("/Users/uw20204/Documents/mini_project_1_UPD/data/clinical", 
                                 recursive=TRUE, pattern = "xml$") #filters to xml files
filenames_clinical <- Filter(function(x) grepl("org_clinical", x), filenames_clinical) #removes org_clinical files

# filters clinical file names to match patient ID of UPD results
# if it prints a patient ID then it does not have clinical data available
filenames_clinical_ID <- NULL
missing_clinical_patient <- NULL
for(i in 1:length(ID_aUPD)){
  if(length(filenames_clinical[grep(pattern = ID_aUPD[[i]], filenames_clinical)]) != 0){
  filenames_clinical_ID <- c(filenames_clinical_ID, 
                             filenames_clinical[grep(pattern = ID_aUPD[[i]], filenames_clinical)])
  }else{
    missing_clinical_patient <- ID_aUPD[[i]]
    setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/data/", dataset, "/aUPD_results_", dataset, "_filtered", sep = ""))
    unlink(filenames_aUPD_results[grep(pattern = missing_clinical_patient, filenames_aUPD_results)]) #deletes the aUPD file for which there is no clinical data
    }
}
return(list(filenames_aUPD_results = filenames_aUPD_results[grep(pattern = missing_clinical_patient, filenames_aUPD_results, invert = TRUE)], 
            ID_aUPD = ID_aUPD[ID_aUPD != missing_clinical_patient], # removes patient with missing clinical data from analysis 
            filenames_clinical = filenames_clinical, filenames_clinical_ID = filenames_clinical_ID, missing_clinical_patient = missing_clinical_patient, ID_end = ID_end))
}

