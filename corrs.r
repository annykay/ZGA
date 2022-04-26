path_ri <- "D:/study/OmicsData/project/totalRNA/"
path_li <- "D:/study/OmicsData/project/LiRibo/"
setwd(path_ri)
filenames_ri  <- list.files(pattern="*.txt", full.names=TRUE)
samples_ri <- sub(".*totalRNA_*", "", sub("*.quant.genes.*", "", filenames_ri))

setwd(path_li)
filenames_li  <- list.files(pattern="*.txt", full.names=TRUE)
samples_li <- sub(".*liRiboseq_*", "", sub("*.quant.genes.*", "", filenames_li))

columns <- c( "gene_id",
 "transcript_id(s)",
 "length",
 "effective_length",
 "expected_count",
 "TPM",
 "FPKM",
 "posterior_mean_count",
 "posterior_standard_deviation_of_count",
 "pme_TPM",
 "pme_FPKM",
 "TPM_ci_lower_bound",
 "TPM_ci_upper_bound",
 "TPM_coefficient_of_quartile_variation",
 "FPKM_ci_lower_bound",
 "FPKM_ci_upper_bound",
 "FPKM_coefficient_of_quartile_variation")
calculateCorrs <- function(path_li, path_ri, filenames_li, stages_li, filenames_ri, stages_ri, column_id){
  setwd(path_ri)
  for (i in c(1:length(filenames_ri))) {
    varname <- paste('ri', stages_ri[i], sep = "_")
    assign(varname, read.table(filenames_ri[i], skip = 1, col.names = columns))
  }

  setwd(path_li)
  for (i in c(1:length(filenames_li))){
    varname <- paste('li', stages_li[i], sep = "_")
    assign(varname, read.table(filenames_li[i], skip = 1, col.names = columns))
  }
  n <- length(stages_ri)
  m <- length(stages_li)
  corrs <- matrix(nrow=n, ncol=m)
  for (i in c(1:length(stages_ri))){
    varname_ri <- paste('ri', stages_ri[i], sep = "_")
    ri_data <- get(varname_ri)
    for (j in c(1:length(stages_li))){
      varname_li <- paste('li', stages_li[j], sep = "_")
      li_data <- get(varname_li)
      total <- merge(ri_data[, c(1,2,column_id)], li_data[, c(1,2,column_id)],
         by = c('gene_id', 'transcript_id.s.'), all.x = F, all.y = F)
      corr <-  cor(total[,3],total[,4], method = 'spearman' )

      corrs[i, j] <- corr

    }
  }
  return(corrs)
}

correlations <- calculateCorrs(path_li, path_ri, filenames_li, samples_li, filenames_ri, samples_ri, 6)
corrplot(correlations, method = 'color')

correlation_1 <- correlations[c(11,12,1,2,3,4,5,6,7,8,9,10), c(3,4,5,6,7,14,8,9,10,11,12,13)]
corrplot(correlation_1, method = 'color')
corrs_2 <- matrix(nrow=6, ncol=6)
for (i in c(1,3,5, 7, 9, 11)){
  for (j in c(1,3,5,7,9,11)){
    av <- 0.25*(correlation_1[i,j] + correlation_1[i+1, j] + correlation_1[i, j+1] + correlation_1[i+1, j+1])
    corrs_2[ i%/%2 +1, j%/%2 + 1] <- av
  }
}

colnames(corrs_2) <- c("MII", "1 Cell", "2 Cell", "4 Cell", "Mor", "BL")
row.names(corrs_2) <- c("MII", "1 Cell", "2 Cell", "4 Cell", "Mor", "BL")

pdf("D:/study/OmicsData/project/correlations.pdf")
corrplot(corrs_2, method = 'color',  addCoef.col = 'grey50', is.corr = FALSE, col.lim = c(min(corrs_2), max(corrs_2)), xlab = 'LiRibo', mar = c(2,2,2,2))

mtext(text = 'LiRibo-Seq', side = 3, line = 2, cex = 2.5, col = 'grey50')
mtext(text = 'totalRNA-Seq', side = 2, line = 2, cex = 2.5, col = 'grey50')
dev.off()
