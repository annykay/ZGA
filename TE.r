path_ri <- "D:/study/OmicsData/project/totalRNA/"
path_li <- "D:/study/OmicsData/project/LiRibo/"
setwd(path_ri)
filenames_ri  <- list.files(pattern="*.txt", full.names=TRUE)
samples_ri <- sub(".*totalRNA_*", "", sub("*.quant.genes.*", "", filenames_ri))

setwd(path_li)
filenames_li  <- list.files(pattern="*.txt", full.names=TRUE)
samples_li <- sub(".*liRiboseq_*", "", sub("*.quant.genes.*", "", filenames_li))

columns <- c( "gene_id",
 "transcript_id.s.",
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
AggregateCounts <- function(path, filenames, stages, column_id, columns){
  setwd(path)
  df_list <- list()
  for (i in c(1:length(filenames))) {
    varname <- paste('ri', stages[i], sep = "_")
    assign(varname, read.table(filenames[i], skip = 1, col.names = columns))

    print(colnames(get(varname)[c(1,2,column_id)]))
    df_list <- list(df_list, get(varname)[c(1,2,column_id)])
    df_list <- unlist(df_list)

  }

  total <- Reduce(function(x, y) merge(x, y, all=TRUE, by = c("gene_id",  "transcript_id.s.")), unlist(df_list)) 
  print('OK')
  colnames(total) <-  c("gene_id",  "transcript_id.s.", stages)
  return(total)
}

total_li <- AggregateCounts(path_li, filenames_li[c(3,4,5,6,7,14,8,9,10,11,12,13)],
samples_li[c(3,4,5,6,7,14,8,9,10,11,12,13)], 7, columns)

total_ri <- AggregateCounts(path_ri, filenames_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)],
samples_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)], 7, columns)
