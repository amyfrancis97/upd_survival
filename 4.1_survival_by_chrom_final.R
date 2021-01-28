# this script calculates survival based on aUPD at a chromosome level
aUPD_surv_analysis <- aUPD_surv_analysis[, 1:7]
# creates a function to find which patients are positive or negative for aUPD for each chromosome
adding_chrom <- function(z){
  new_tcn <- tcn_location_dataframe_upd %>% filter(Chrom == z)
  aUPD_pos <- unique(new_tcn[, 1])
  aUPD_neg <- setdiff(reading_data_output$ID_aUPD, aUPD_pos)
  aUPD_surv_analysis[, as.character(z)] <- "x"
  
  for(i in aUPD_pos){
    aUPD_surv_analysis[aUPD_surv_analysis$patient_id == i, ][, as.character(z)] <- "Positive"
  }
  
  for(i in aUPD_neg){
    aUPD_surv_analysis[aUPD_surv_analysis$patient_id == i, ][, as.character(z)] <- "Negative"
  }
  return(aUPD_surv_analysis)
}

# applies the function to each chromosome
for(i in 1:22){
  aUPD_surv_analysis <- adding_chrom(i)
}

# adds prefix to colnames (could be done in function)
colnames(aUPD_surv_analysis)[8:ncol(aUPD_surv_analysis)] <- 
  paste("Chr",colnames(aUPD_surv_analysis)[8:ncol(aUPD_surv_analysis)], sep = "")

# survival grouped by chromosome neg and pos
surv_object <- Surv(time = aUPD_surv_analysis$days_last_fu, event = aUPD_surv_analysis$censoring)
aUPD_surv_analysis
plot_chrom <- NULL
for(chrom in colnames(aUPD_surv_analysis)[8:length(colnames(aUPD_surv_analysis))]){
  f1 <- survfit(surv_object ~ aUPD_surv_analysis[, chrom], data = aUPD_surv_analysis)
  plot_chrom[[chrom]] <- ggsurvplot(f1, data = aUPD_surv_analysis, pval = TRUE, pval.method =TRUE ,title = chrom,
                                    legend.labs = c("Neg", "Pos"), legend.title = "aUPD Status",
                                    xlab = "Time (Days)", ggtheme = theme_cowplot())
}

surv_plot <- function(dataset){
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
  surv_plots_arrange_save <- function(a, b, c, d){
    plots_whole_chrom <- arrange_ggsurvplots(list(plot_chrom[[a]], plot_chrom[[b]], plot_chrom[[c]],
                                                  plot_chrom[[d]]),
                                             print = FALSE,
                                             title = NA,
                                             ncol = 2,
                                             nrow = 2,
                                             surv.plot.height = NULL,
                                             risk.table.height = NULL,
                                             ncensor.plot.height = NULL
    )
    ggsave(paste(a, b, c, d, sep = "", "survplot", ".pdf"), width = 10, height = 10, units = "in", plots_whole_chrom)
  }
  surv_plots_arrange_save(1, 2, 3, 4)
  surv_plots_arrange_save(5, 6, 7, 8)
  surv_plots_arrange_save(9, 10, 11, 12)
  surv_plots_arrange_save(13, 14, 15, 16)
  surv_plots_arrange_save(17, 18, 19, 20)
  
  plots_whole_chrom <- arrange_ggsurvplots(list(plot_chrom[[21]], plot_chrom[[22]]),
                                           print = FALSE,
                                           title = NA,
                                           ncol = 2,
                                           nrow = 2,
                                           surv.plot.height = NULL,
                                           risk.table.height = NULL,
                                           ncensor.plot.height = NULL
  )
  ggsave(paste(21, 22, sep = "", "survplot", ".pdf"),  width = 10, height = 10, units = "in", plots_whole_chrom)
  }
