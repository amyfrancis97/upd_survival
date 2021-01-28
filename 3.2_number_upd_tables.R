# saves table with number aUPD events per chromosome
# saves table with number aUPD events per chromosome
tcn_location_dataframe_upd <- tcn_location_dataframe %>%
  filter(Flagsign == "upd")

number_upd_chrom_number_patients_upd <- function(dataset){
# % of UPD events for each chromosome
perc_chrom_df <- NULL
for(i in 1:22){
  perc_chrom <- nrow(tcn_location_dataframe_upd[tcn_location_dataframe_upd$Chrom == i, ])/ nrow(tcn_location_dataframe_upd) *100
  perc_chrom_df <- round(rbind(perc_chrom_df, data.frame(i, nrow(tcn_location_dataframe_upd[tcn_location_dataframe_upd$Chrom == i, ]), perc_chrom)), digits = 2)
}
colnames(perc_chrom_df) <- c("Chr_no", "no_upd_events", "perc_aUPD_events")

setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
write.csv(perc_chrom_df[order(perc_chrom_df$perc_aUPD_events), ], "number_upd_events_per_chrom.csv")
}
