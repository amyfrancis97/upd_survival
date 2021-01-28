dataset_genes_COSMIC
saving_COSMIC <- function(dataset_genes_cosmic, dataset, region_description){
  #this onyl extracts info from segmental aUPD since the dataset use filters for this
  patients_per_mutated_gene <- function(dataset_genes_cosmic){
    sample_IDs_positive_unique_df <- NULL
    for(gene_positive in 1:length(unique(dataset_genes_cosmic$GENE_NAME))){
      sample_IDs_positive_gene <- dataset_genes_cosmic$SAMPLE_NAME[which(dataset_genes_cosmic$GENE_NAME ==
                                                                           unique(dataset_genes_cosmic$GENE_NAME)[gene_positive])]
      newdataset <- dataset_genes_cosmic[which(dataset_genes_cosmic$GENE_NAME ==
                                                             unique(dataset_genes_cosmic$GENE_NAME)[gene_positive]), ]
      number_of_pos_patients_with_mutation <- length(unique(sample_IDs_positive_gene))
      
      if(number_of_pos_patients_with_mutation == 1){
        patient_ID <- unique(sample_IDs_positive_gene)
        number_mutations_each_gene <- nrow(newdataset[newdataset$SAMPLE_NAME == sample_IDs_positive_gene, ])

      }else{
        number_mutations_each_gene <- NULL
        for(sample_name in unique(sample_IDs_positive_gene)){
          new_dataset <- dataset_genes_cosmic[which(dataset_genes_cosmic$GENE_NAME ==
                                                      unique(dataset_genes_cosmic$GENE_NAME)[gene_positive]), ]
          number_mutations_per_gene_per_patient <- nrow(new_dataset[new_dataset$SAMPLE_NAME == sample_name, ])
          number_mutations_each_gene[[sample_name]] <- number_mutations_per_gene_per_patient
        }
        number_mutations_each_gene <- paste(number_mutations_each_gene[1], number_mutations_each_gene[2], sep = ", ")
        
        patient_ID <- paste(unique(sample_IDs_positive_gene)[1], unique(sample_IDs_positive_gene)[2], sep = ",")
      }

      sample_IDs_positive_unique_df <- rbind(sample_IDs_positive_unique_df, data.frame(gene = unique(dataset_genes_cosmic$GENE_NAME)[gene_positive], 
                                                                                       number_of_pos_patients_with_mutation = number_of_pos_patients_with_mutation, 
                                                                                       patient_ID = patient_ID, number_mutations_per_gene = number_mutations_each_gene))
      } 
    return(sample_IDs_positive_unique_df)
  }
  
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
  write.csv(patients_per_mutated_gene(dataset_genes_cosmic[[1]]), file = paste(dataset, region_description, "_mutation_upd_positive_patients.csv"))
  write.csv(patients_per_mutated_gene(dataset_genes_cosmic[[2]]), file = paste(dataset, region_description, "_mutation_upd_negative_patients.csv"))
}

