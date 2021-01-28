# this script calculates survival based on aUPD at a chromosome level
tcn_location_dataframe_whole_chrom <- tcn_location_dataframe_upd %>%
  filter(aUPD_perc >= 90)
tcn_location_dataframe_segmental <- tcn_location_dataframe_upd %>%
  filter(aUPD_perc < 90)


segmental_whole_chrom_surv <- function(dataset, upd_type, upd_tcn_dataframe){
aUPD_surv_analysis <- aUPD_surv_analysis[, c(1:7)]
# creates a function to find which patients are positive or negative for aUPD for each chromosome
adding_chrom <- function(z){
  new_tcn <- tcn_location_dataframe_whole_chrom %>% filter(Chrom == z)
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
  paste("Chr", colnames(aUPD_surv_analysis)[8:ncol(aUPD_surv_analysis)], sep = "")

# plotting
surv_object <- Surv(time = aUPD_surv_analysis$days_last_fu, event = aUPD_surv_analysis$censoring)
plot_chrom <- NULL
for(chrom in colnames(aUPD_surv_analysis)[8:length(colnames(aUPD_surv_analysis))]){
  if(length(which(aUPD_surv_analysis[, chrom] == "Positive")) > 0){
    f1 <- survfit(surv_object ~ aUPD_surv_analysis[, chrom], data = aUPD_surv_analysis)
    plot_chrom[[chrom]] <- ggsurvplot(f1, data = aUPD_surv_analysis, pval = TRUE, pval.method =TRUE ,title = chrom,
                                      legend.labs = c("Neg", "Pos"), legend.title = "aUPD Status",
                                      xlab = "Time (Days)", ggtheme = theme_cowplot())
  }
}

setwd("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/READ")
ggsave("chr1_whole_chrom_surv_plot.jpeg", width = 5, height =5, units = "in", print(plot_todo))
setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep =""))
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
  ggsave(paste(a, b, c, d, sep = "", "survplot_", upd_type, ".pdf"), width = 10, height = 10, units = "in", plots_whole_chrom)
}

for(i in c(1, 5, 9, 13, 17)){
surv_plots_arrange_save(names(plot_chrom)[i], names(plot_chrom)[i+1], names(plot_chrom)[i+2], names(plot_chrom)[i+3])
}
  
plots_whole_chrom <- arrange_ggsurvplots(list(plot_chrom[[21]], plot_chrom[[22]]),
                                         print = FALSE,
                                         title = NA,
                                         ncol = 2,
                                         nrow = 1,
                                         surv.plot.height = NULL,
                                         risk.table.height = NULL,
                                         ncensor.plot.height = NULL
)
ggsave(paste("Chr21_22_", sep = "", "survplot_", upd_type, ".pdf"), width = 10, height = 5, units = "in", plots_whole_chrom)
##coxph
##

chrom_names <- names(plot_chrom)

fit.coxph <- purrr::map(chrom_names, ~coxph(as.formula(paste("surv_object ~ ", .x)), aUPD_surv_analysis))
dev.off()
forest <- NULL
for(i in 1:22){
  forest[[i]] <- ggforest(fit.coxph[[i]], data = aUPD_surv_analysis, main = paste("Hazard Ratio", "Chrom", i, upd_type))
}
ggforest(fit.coxph[[1]], data = aUPD_surv_analysis, main = paste("Hazard Ratio", "Chrom", "1", "segmental"))
fit.coxph[[1]]
aUPD_surv_analysis
for(i in c(1, 5, 9, 13, 17)){
ggsave(paste(upd_type, "coxph_plot", i, i+3, ".pdf", sep = "_"), ggarrange(forest[[i]], forest[[i+1]], forest[[i+2]], forest[[i+3]],
                                                                   ncol = 2, nrow = 2
),  width = 12, height = 7, units = "in")
}

ggsave(ggarrange(forest[[21]], forest[[22]],
          ncol = 2, nrow = 1
), width = 12, height = 3.5, units = "in")
  
}
