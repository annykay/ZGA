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
      varname_li <- paste('li', stages_li[i], sep = "_")
      li_data <- get(varname_li)
      total <- merge(ri_data[, c(1,2,column_id)], li_data[, c(1,2,column_id)],
         by = c('gene_id', 'transcript_id.s.'), all.x = F, all.y = F)
      corr <-  cor(total[,3],total[,4])

      corrs[i, j] <- corr

    }
  }
  return(corrs)
}

calculateCorrs(path_li, path_ri, filenames_li, samples_li, filenames_ri, samples_ri, 6)
