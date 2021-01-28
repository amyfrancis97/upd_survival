# this script calculates the percentage coverage of the total length of the chromosome for each aUPD event
# adds a new column to tcn_location_dataframe which shows the percentage length of the event


chrom_lengths <- getChromInfoFromUCSC("hg19")

#new column for aUPD type (e.g. >80% is whole chromosome, <80% is segmental)
aUPD_type <- function(chrom_num, chrom_col_no_UCSC){
percentage_aUPD <- round((filter(tcn_location_dataframe_upd, Chrom == chrom_num)$End- 
                            filter(tcn_location_dataframe_upd, Chrom == chrom_num)$Start)/
                           chrom_lengths[chrom_lengths$chrom == chrom_col_no_UCSC, ]$size*100,
                         digits = 2)
return(percentage_aUPD)
  }

tcn_location_dataframe_upd$aUPD_perc <- NULL
for(i in 1:22){
aUPD_perc <- aUPD_type(i, chrom_lengths$chrom[i])
tcn_location_dataframe_upd[tcn_location_dataframe_upd$Chrom == i, "aUPD_perc"] <- aUPD_perc
}


#total number of aUPD regions:
total_aUPD <- sum(aUPD_surv_analysis$no_upd_events)
mean_aUPD <- mean(aUPD_surv_analysis$no_upd_events)
median_aUPD <- median(aUPD_surv_analysis$no_upd_events)
range_aUPD <- range(aUPD_surv_analysis$no_upd_events)
max_minus_min_aUPD <- max(aUPD_surv_analysis$no_upd_events) - min(aUPD_surv_analysis$no_upd_events)
df <- data.frame(total_aUPD, round(mean_aUPD, digits = 2), median_aUPD, range = "1-17", max_minus_min_aUPD)

perc_aupd_coverage <- function(dataset){
setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
write.csv(df, file = "stats_df.csv")

jpeg("perc_coverage_upd.jpeg", units = "in", width = 6, height = 5, res = 300)
plot.new()
plot(tcn_location_dataframe_upd$Chrom, tcn_location_dataframe_upd$aUPD_perc, xlab = "Chromosome",
     ylab = "aUPD % coverage of total chromosome")
abline(h = 90, col = "red", lwd = 2)
dev.off()
# percentage whole chrom & segmental
#90% threshold

segmental_num <- as.character(length(which(tcn_location_dataframe_upd$aUPD_perc < 90)))
segmental_perc <- paste("(", as.character(round(length(which(tcn_location_dataframe_upd$aUPD_perc < 90))/total_aUPD*100, digits = 2)), "%", ")", sep = "")
whole_chrom_num <- as.character(length(which(tcn_location_dataframe_upd$aUPD_perc >= 90)))
whole_chrom_perc <- paste("(", as.character(round(length(which(tcn_location_dataframe_upd$aUPD_perc >= 90))/total_aUPD*100, digits = 2)),"%", ")", sep = "")
paste(segmental_num, segmental_perc, sep = " ")
perc_seg_wc_table <- data.frame(total_aupd = total_aUPD,segmental = paste(segmental_num, segmental_perc, sep = " "), whole_chromosome = paste(whole_chrom_num, whole_chrom_perc, sep = " "))
write.csv(perc_seg_wc_table, file = "perc_seg_wc_df.csv")

}