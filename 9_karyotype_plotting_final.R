# visualising UPD regions across chromosome
karyoplot_distribution_upd <- function(dataset){
setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
tcn_location_dataframe_upd_2 <- tcn_location_dataframe_upd
tcn_location_dataframe_upd_2$Chrom <- paste("chr", tcn_location_dataframe_upd_2$Chrom, sep ="")
granges_obj_2 <- GRanges(tcn_location_dataframe_upd_2)
#looking at SORs of data

cov_gr_2 <- GRanges(coverage(granges_obj_2))
cov_2 <- coverage(cov_gr_2, weight="score")
 
jpeg("aupd_coverage_karyoplot_1_11.jpeg", units = "in", width = 15, height = 8, res = 300)
plot.new()
kp <- plotKaryotype(chromosomes = paste("chr", as.character(c(1:11)), sep = ""))
kpAddBaseNumbers(kp)
kpAddCytobandLabels(kp)
kpPlotRegions(kp, data=granges_obj_2, layer.margin = 0.05, lwd = 0.0, col = "darkgreen")
dev.off()

jpeg("aupd_coverage_karyoplot_12_22.jpeg", units = "in", width = 15, height = 8, res = 300)
plot.new()
kp <- plotKaryotype(chromosomes = paste("chr", as.character(c(12:22)), sep = ""))
kpAddBaseNumbers(kp)
kpAddCytobandLabels(kp)
kpPlotRegions(kp, data=granges_obj_2, layer.margin = 0.05, lwd = 0.0, col = "darkgreen")
dev.off()
}
getwd()
