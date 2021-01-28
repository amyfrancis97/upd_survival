# this script uses genes discovered from sliding window analysis and previous script
# to plot onto karyotype plot

genes_w_ranges$midpoint <- genes_w_ranges$start + (0.5*genes_w_ranges$width)

# this next part copies script 12 but adds the genes into the plot
# this script generates survival p-values from windows of 2M bases


# filter out whole chromosomal aUPD- only looking at segmental aUPD for regional analysis
segmental_tcn <- tcn_location_dataframe_upd[tcn_location_dataframe_upd$aUPD_perc < 90, ]
granges_obj <- GRanges(segmental_tcn)
granges_df_for_window <- data.frame(granges_obj)
getChromInfoFromUCSC("hg19")[1:22, ]$size[5]

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

# function carries out analysis for each chromosome individually

plot_genes_karyotype <- function(chrom_number, genes_w_ranges_df){
  # sources a new aUPD_surv_analysis_2 df in which to add new columns
  aUPD_surv_analysis_2 <- aUPD_surv_analysis[, c(1:7)]
  # uses above parameters to make a vector of the time steps at which each window will end
  steps <- seq(start + w, getChromInfoFromUCSC("hg19")[1:22, ]$size[5], by = w)
  steps <- data.frame(cbind(c(0, steps), c(steps-1, getChromInfoFromUCSC("hg19")[1:22, ]$size[5])))
  steps$chrom <- 5
  colnames(steps) <- c("start", "end", "chrom")
  
  # adds new regions to aUPD_surv_analysis_2 as columns and whether the patient is aUPD positive
  # or negative at that region using window_analysis function previously made
  window_numbers <- c(1:nrow(steps))
  window_names <- c(1:nrow(steps))
  window_names <- paste("window", window_names, sep = "")
  for(i in 1:nrow(steps)){
    aUPD_surv_analysis_2 <- window_analysis(window_names[i], window_numbers[i])
  }
  aUPD_surv_analysis_2
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
  
  # generates region mid point to plot
  regional_df <- cbind(steps, statistical_vals[2:5])
  regional_df$chrom <- gsub(" ", "", paste("chr", as.numeric(regional_df$chrom, sep = NULL)))
  midpoint <- (regional_df$start + regional_df$end)/2
  regional_df <- cbind(regional_df[1:3], midpoint, regional_df[4:ncol(regional_df)])
  
  new_df_fin <- NULL
  for(region_name in colnames(aUPD_surv_analysis_2)[8: length(colnames(aUPD_surv_analysis_2))]){
    new_df <- data.frame(positive_patients = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Positive", ][, region_name]),
                         negative_patients = length(aUPD_surv_analysis_2[aUPD_surv_analysis_2[, region_name] == "Negative", ][, region_name]))
    new_df_fin <- rbind(new_df_fin, new_df)
  }
  regional_df <- cbind(regional_df, new_df_fin)
  regional_df
  # write csv for further analysis
  setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs")
  write.csv(regional_df, file = paste(5, "p_val_table_sliding_window.csv", sep = "_"))
  
  # convert p-values to -log10(p-values)
  regional_df[c(5:8)] <- -log10(regional_df[c(5:8)])
  regional_df
  # sets infinite and NA values to 0
  regional_df <- regional_df %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  regional_df <- regional_df %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x))
  regional_df$pval_means <- rowMeans(regional_df[, c(5, 7)])
  
  # karyotype plotting windows and p-values
  dev.off()
   setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs")
  jpeg(filename = paste(5, "plot_windows.jpeg"), width=12,
       height=6,
       units="in",
       res=500)
  plot.new()
  
  kp <- plotKaryotype(chromosomes = paste("chr", as.character(5), sep = ""), plot.type=2,
                      zoom = toGRanges(data.frame(paste("chr", as.character(5), sep = ""), mydata$start[1]-1000, mydata$end[length(mydata$end)]+1000)))
  kpAddBaseNumbers(kp)
  kpAddCytobandLabels(kp)
  kpPlotMarkers(kp, chr = oncogenes_ROI[, 1], x = oncogenes_ROI[, 2], 
                labels = oncogenes_ROI[, 8], data.panel=1, cex=0.8,
                ymax = 3.5, ignore.chromosome.ends = FALSE, label.color = "red", r0 = 0.05, r1 = 1.2)
  
  
  # plots genes in the region
  kpPlotMarkers(kp, chr = TSG_ROI[, 1], x = TSG_ROI[, 2], 
                labels = TSG_ROI[, 8], data.panel=2, cex=0.8,
                ymax = 3.5, ignore.chromosome.ends = FALSE, label.color = "purple")

   
  legend(x = "bottomleft", legend = 
           c("Oncogenes", "TSG's"),
         title =  "Gene type",
         fill = c("red", "purple"),
         ncol = 2, cex = 0.8, bty = "n")
  dev.off()
  regional_df_list[[chrom_number]] <- regional_df
}

# applies above funciton to different chromosome regions above genome wide significance threshold
plot_genes_karyotype(5, genes_w_ranges)
regional_df[which(regional_df$start == 0):which(regional_df$start == 820000), ]
regional_df
