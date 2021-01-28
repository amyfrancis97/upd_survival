# this script generates survival p-values from windows of 2M bases
# uses "sliding window" method for regional analysis, rather than coverage from GenomicRanges
window_analysis <- function(window_id, window_number){
  # tells you which regions of UPD overlap with the first window
  new_df <- subsetByOverlaps(granges_obj, GRanges(steps[window_number, ]))
  
  #extracts patients who have aUPD in this first window
  aUPD_pos <- new_df$patient_id
  #extracts patients who do not have aUPD in this window
  aUPD_neg <- setdiff(reading_data_output$ID_aUPD, aUPD_pos)
  #creates a new columns in surv_analysis for this first window
  aUPD_surv_analysis_2[, window_id] <- "x"
  length(aUPD_pos) #26 positive
  length(aUPD_neg) #69 negative
  
  for(i in aUPD_pos){
    aUPD_surv_analysis_2[aUPD_surv_analysis_2$patient_id == i, ][, window_id] <- "Positive"
  }
  for(i in aUPD_neg){
    aUPD_surv_analysis_2[aUPD_surv_analysis_2$patient_id == i, ][, window_id] <- "Negative"
  }
  return(aUPD_surv_analysis_2)
}

sliding_window_analysis <- function(dataset){
# filter out whole chromosomal aUPD- only looking at segmental aUPD for regional analysis
segmental_tcn <- tcn_location_dataframe_upd[tcn_location_dataframe_upd$aUPD_perc < 90, ]
granges_obj <- GRanges(segmental_tcn)
granges_df_for_window <- data.frame(granges_obj)
granges_obj
# set a couple of parameters we'll use to slice the series into chunks:
# window width (w) and the time step at which you want to end the first
# training set
w <- 2000000 
start <- 1

# create a function to add new columns to the aUPD_surv_analysis_2 df
# each column represents a new 2MB region specified by the above method
# the contents of the column is either "negative" for patients who dont have aUPD at that region
# and "positive" for those that do
# each row represents a different patient

# for loop carries out analysis for each chromosome individually
regional_df_list <- NULL
regions_patients_id <- NULL #saves the patient ID's that are positive for the regions for further mutation analysis

for(chrom_number in c(1:22)){
  aUPD_surv_analysis_2 <- aUPD_surv_analysis[, c(1:7)]
  # sources a new aUPD_surv_analysis_2 df in which to add new columns
  # uses above parameters to make a vector of the time steps at which each window will end
  steps <- seq(start + w, getChromInfoFromUCSC("hg19")[1:22, ]$size[chrom_number], by = w)
  steps <- data.frame(cbind(c(0, steps), c(steps-1, getChromInfoFromUCSC("hg19")[1:22, ]$size[chrom_number])))
  steps$chrom <- chrom_number
  colnames(steps) <- c("start", "end", "chrom")
  aUPD_surv_analysis_2
  # adds new regions to aUPD_surv_analysis_2 as columns and whether the patient is aUPD positive
  # or negative at that region using window_analysis function previously made
  window_numbers <- c(1:nrow(steps))
  window_names <- c(1:nrow(steps))
  window_names <- paste("window", window_names, sep = "")
  for(i in 1:nrow(steps)){
    aUPD_surv_analysis_2 <- window_analysis(window_names[i], window_numbers[i])
    }
  
  # calculates p-values for survival for each of the new regions
  statistical_vals <- NULL
  for(c in 8:length(colnames(aUPD_surv_analysis_2))){
    if(sum(aUPD_surv_analysis_2[, c] == "Positive") > 0){#this stops code from breaking if no patients have UPD in a particular region
      form <- as.formula(paste("surv_object", colnames(aUPD_surv_analysis_2)[c], sep = " ~ "))
      f1 <- survfit(form, data = aUPD_surv_analysis_2)
      coxph_test <- coxph(form, data = aUPD_surv_analysis_2)
      surv_log_rank <- survdiff(form, data =aUPD_surv_analysis_2)
      surv_log_rank <- pchisq(surv_log_rank$chisq, length(surv_log_rank$n)-1, lower.tail = FALSE)
      likelihood_ratio <- as.numeric(summary(coxph_test)$logtest[3])
      wald_test_stat <- summary(coxph_test)$coefficients[5]
      log_stat <- summary(coxph_test)$conf.int[2]
    }else{ # if no patients have aUPD at a region then it assigns p-values as zero
      surv_log_rank <- 0
      likelihood_ratio <- 0
      wald_test_stat <- 0
      log_stat <- 0
    }
    # generates new df comprising regions and p-values
    statistical_vals <- rbind(statistical_vals, data.frame(Region = colnames(aUPD_surv_analysis_2)[c],
                                                           surv_log_rank = surv_log_rank,
                                                           coxph_wald = wald_test_stat,
                                                           coxph_likelihood = likelihood_ratio, coxph_log_ratio = log_stat))
  }
  
  # adjusts p-values using "fdr"
  statistical_vals$surv_log_rank <-  p.adjust(statistical_vals$surv_log_rank, method = "fdr")
  statistical_vals$coxph_wald <-  p.adjust(statistical_vals$coxph_wald, method = "fdr")
  statistical_vals$coxph_likelihood <-  p.adjust(statistical_vals$coxph_likelihood, method = "fdr")
  statistical_vals$coxph_log_ratio <-  p.adjust(statistical_vals$coxph_log_ratio, method = "fdr")
  regional_df
  # generates region mid point to plot
  regional_df <- cbind(steps, statistical_vals[2:5])
  regional_df$chrom <- gsub(" ", "", paste("chr", as.numeric(regional_df$chrom, sep = NULL)))
  midpoint <- (regional_df$start + regional_df$end)/2
  regional_df <- cbind(regional_df[1:3], midpoint, regional_df[4:ncol(regional_df)])
  
  new_df_fin <- NULL
  for(region_name in colnames(aUPD_surv_analysis_2)[8: length(colnames(aUPD_surv_analysis_2))]){
    new_df <- data.frame(positive_patients = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Positive", ][, region_name]),
                         negative_patients = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Negative", ][, region_name]),
                         positive_alive = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Positive" & aUPD_surv_analysis_2$censoring == 1, ][, region_name]),
                         positive_dead = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Positive" & aUPD_surv_analysis_2$censoring == 2, ][, region_name]),
                         negative_alive = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Negative" & aUPD_surv_analysis_2$censoring == 1, ][, region_name]),
                         negative_dead = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Negative" & aUPD_surv_analysis_2$censoring == 2, ][, region_name]))
    new_df_fin <- rbind(new_df_fin, new_df)
  }
  regional_df <- cbind(regional_df, new_df_fin)
  regional_df
  # write csv for further analysis
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
  write.csv(regional_df, file = paste(chrom_number, "p_val_table_sliding_window.csv", sep = "_"))
  
  # convert p-values to -log10(p-values)
  regional_df[c(5:8)] <- -log10(regional_df[c(5:8)])
  
  # sets infinite and NA values to 0
  regional_df <- regional_df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  regional_df <- regional_df %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x))
  regional_df$pval_means <- rowMeans(regional_df[, c(5, 7)])
  
  # karyotype plotting windows and p-values
  jpeg(filename = paste(chrom_number, "plot_windows_90_threhsold.jpeg"), width=8,
       height=9,
       units="in",
       res=500)
  plot.new()
  kp <- plotKaryotype(chromosomes = paste("chr", as.character(c(chrom_number)), sep = ""))
  kpAddBaseNumbers(kp)
  kpAddCytobandLabels(kp)
  kpDataBackground(kp, data.panel = 1, r0=0.05, r1=0.4)
  kpAxis(kp, ymin=0, ymax=30, 
         r0=0.05, r1=0.4, col="gray50", cex=0.5)
  # plot p-values against chromosomal location
  kpPoints(kp, chr = regional_df$chrom, x=regional_df$midpoint,
           y = regional_df$surv_log_rank,ymin=0, ymax=30, r0=0.05, r1=0.4, 
           col="midnightblue", pch=16, cex=0.5)
  
  # adds genome wide significance threshold to plot
  kpLines(kp, chr = regional_df$chrom, x=c(0, regional_df$midpoint),
          y = rep(7.3, length(regional_df$coxph_log_ratio)),ymin=0, ymax=30, r0=0.05, r1=0.4, 
          col="red", pch=16, cex=0.5)
  # plots number of patients with UPD at that region
  kpAxis(kp, ymin=0, ymax=nrow(aUPD_surv_analysis), r0=0.45, r1=0.9, col="gray50", cex=0.5, numticks = 5)
  kpBars(kp, chr = regional_df$chrom, x0=regional_df$start, x1 = regional_df$end, 
         y1 = regional_df$positive_patients, ymin=0, ymax=nrow(aUPD_surv_analysis), r0=0.45, r1=0.9, col = "thistle1")
  kpBars(kp, chr = regional_df$chrom, x0=regional_df$start, x1 = regional_df$end, 
         y1 = regional_df$positive_dead, ymin=0, ymax=nrow(aUPD_surv_analysis), r0=0.45, r1=0.9, col = "orchid4")
  
  kpBars(kp, chr = regional_df$chrom, x0=regional_df$start, x1 = regional_df$end, y0 = regional_df$positive_patients,
         y1 = 116, ymin= 0, ymax=nrow(aUPD_surv_analysis), r0=0.45, r1=0.9, col = "lightsteelblue1")
  kpBars(kp, chr = regional_df$chrom, x0=regional_df$start, x1 = regional_df$end, y0 = regional_df$positive_patients,
         y1 = regional_df$positive_patients+regional_df$negative_dead, ymin= 0, ymax=nrow(aUPD_surv_analysis), r0=0.45, r1=0.9, col = "lightslategrey")
  
  legend(x = "topright", legend = 
           c("Negative Alive", "Negative Dead","Positive Alive","Positive Dead", 
             "Survival: Log Rank"),
         title =  "Number of Patients                                   P-Value Methods",
         fill = c("lightsteelblue1", "lightslategrey", "thistle1", "orchid4","midnightblue"),
         ncol = 3, cex = 0.6)
  dev.off()
  regional_df_list[[chrom_number]] <- regional_df
  regions_patients_id[[chrom_number]] <- aUPD_surv_analysis_2
}
return(list(regional_df_list, regions_patients_id))
}
