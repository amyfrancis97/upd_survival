# this script finds overlapping genes with the regions from sliding window analysis
# this again is a play-around script and not fully automated yet
# uses genome wide significance threshold to pull out regions of interest
# then finds regions of overlap with the txdb (hg19) genome



regional_df_list
# extracts p values that are greater than or equal to -log10(genome-wide significance threshold)
#chr5
mydata <-  regionaldf_list_patient_id[[1]][[5]][regionaldf_list_patient_id[[1]][[5]]$surv_log_rank >= -log10(5e-8), ]
plot(mydata$midpoint, mydata$surv_log_rank, cex = 1.5)
mydata$chrom <- "chr5"
chrom_p_val_table <- mydata
chrom_ranges <- chrom_p_val_table
chrom_ranges <- na.omit(chrom_ranges)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
seqlevels(txdb) <- seqlevels0(txdb)
seqlevels(txdb) <- paste("chr", 5, sep = "")
txdb_df <- data.frame(genes(txdb, single.strand.genes.only = FALSE)) #finds genes and their locations for txdb

##################################################################################
# for genomic regions that are too large to search ncbi all at once, must split df using code below
#splits larger dataframes as they cannot all be searched together on NCBI
# extracts p values that are greater than or equal to -log10(genome-wide significance threshold)

n <- 6
nr <- nrow(mydata)
mydata_split <- split(mydata, rep(1:ceiling(nr/n), each=n, length.out=nr))

genes_w_ranges <- NULL
for(splitdf in 1:length(mydata_split)){
  overlapping_genes <- subsetByOverlaps(GRanges(txdb_df), 
                                        GRanges(mydata_split[[splitdf]]), ignore.strand = TRUE,
                                        type = "any")
  id_return <- entrez_summary(db="gene", id=data.frame(overlapping_genes)$group_name) 
  genes_w_ranges <- rbind(genes_w_ranges, data.frame(overlapping_genes, t(data.frame(extract_from_esummary(
    id_return, c("name", "nomenclaturesymbol","description", "summary"))))))
}

#removes uncharacterised genes from analysis
genes_w_ranges <- genes_w_ranges[genes_w_ranges$nomenclaturesymbol != "", ]
genes_w_ranges <- genes_w_ranges[genes_w_ranges$summary != "", ]

#removed miRNA genes from analysis
genes_w_ranges <- genes_w_ranges[str_detect(genes_w_ranges$description, "microRNA") == FALSE, ]

for(i in colnames(genes_w_ranges)[8:11]){
  genes_w_ranges[, i] <- as.character(genes_w_ranges[, i])
}

setwd("/Users/uw20204/Documents/mini_project_1_UPD/data")
oncogenes <- read.csv(file = "Human_oncogenes_COSMIC.csv", header = TRUE)
tsg <- read.csv(file = "Human_TSGs_database.csv", header = TRUE, sep = "\t")

#tumour supressor genes and oncogenes within the region of interest
oncogenes_ROI <- genes_w_ranges[genes_w_ranges$name %in% oncogenes$Gene.Symbol, ]
TSG_ROI <- genes_w_ranges[genes_w_ranges$name %in% tsg$GeneSymbol, ]
genes_of_interest <- rbind(oncogenes_ROI, TSG_ROI)

