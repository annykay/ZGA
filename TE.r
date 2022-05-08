library(reshape2)
library(ggplot2)
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
  N <- length(filenames)
  df_list <- vector("list", )
  j <- 1
  for (i in c(1:N)) {
    varname <- paste('ri', stages[i], sep = "_")
    assign(varname, read.table(filenames[i], skip = 1, col.names = columns))

    print(colnames(get(varname)[c(1,2,column_id)]))
    df_list[[j]] <- get(varname)[c(1,2,column_id)]
    j <- j+1

  }

  total <- Reduce(function(x, y) merge(x, y, all=TRUE, by = c("gene_id", "transcript_id.s.")), df_list)
  print('OK')
  colnames(total) <-  c("gene_id",  "transcript_id.s.", stages)
  return(total)
}


total_li <- AggregateCounts(path_li, filenames_li[c(3,4,5,6,7,14,8,9,10,11,12,13)],
samples_li[c(3,4,5,6,7,14,8,9,10,11,12,13)], 6, columns)

total_ri <- AggregateCounts(path_ri, filenames_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)],
samples_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)], 6, columns)

total <- merge(total_ri[apply(total_ri >1, 1,all),], total_li, by =  c("gene_id", "transcript_id.s."))
#total <- merge(total_ri, total_li, by =  c("gene_id", "transcript_id.s."))
all_stages <- unique(sub( "*_rep1", "", sub("*_rep2", "", samples_li[c(3,4,5,6,7,14,8,9,10,11,12,13)])))
all_stages_TE <- paste(all_stages, '_TE', sep="")
for (i in all_stages){
  rep1x <- paste(i, '_rep1', '.x', sep="")
  rep1y <- paste(i, '_rep1', '.y', sep="")
  rep2x <- paste(i, '_rep2', '.x', sep="")
  rep2y <- paste(i, '_rep2', '.y', sep="")
  te <- paste(i, '_TE', sep="")
  total[te] <- (total[rep1y] / total[rep1x] + total[rep2y] / total[rep2x]) /2

}
total_te <- total[, all_stages_TE]
#total.normalized <- t(apply(total_te, 1, function (x) (x - mean(x))/std(x)))
#total.log <- log10(total_te)
for_violin_plot <- melt(total_te)
#for_violin_plot <- for_violin_plot[,-1]

ggplot(for_violin_plot, aes(x=variable, y=value))+ geom_boxplot() + ylim(0,6)
