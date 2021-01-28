# this script calculates cox proportional hazards model on a chromosomal level
coxph_chrom <- function(dataset){
  setwd(paste("/Users/uw20204/Documents/mini_project_1_UPD/UPD_survival/final_ouputs/", dataset, sep = ""))
chrom_names <- colnames(aUPD_surv_analysis)[8:ncol(aUPD_surv_analysis)]

dev.off()
forest <- NULL
for(i in 1:22){
  surv_minus_x <- aUPD_surv_analysis[aUPD_surv_analysis[, 7+1] != "x", ]
  fit.coxph <- purrr::map(chrom_names, ~coxph(as.formula(paste("surv_object ~ ", .x)), surv_minus_x))
  forest[[i]] <- ggforest(fit.coxph[[i]], data = surv_minus_x, main = paste("Hazard Ratio", "Chrom",i))
  ggsave(paste(i, "coxhaz.pdf"), forest[[i]], height = 4, width = 8, units = "in")
  if(isTRUE(class(fit.coxph)=="try-error")) { next }
}
}

