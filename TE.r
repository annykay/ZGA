library(reshape2)
library(ggplot2)
library(ggpubr)

#THIS PART IS PRESENTED IN CORRS.R. CAN BE REPALCED IN REPORT. START
#Setting path to files
path_ri <- "D:/study/OmicsData/project/totalRNA/"
path_li <- "D:/study/OmicsData/project/LiRibo/"

#Reading filenames
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
#THIS PART IS PRESENTED IN CORRS.R. CAN BE REPALCED IN REPORT. END

#Aggregating counts for different samples
AggregateCounts <- function(path, filenames, stages, column_id, columns){
  setwd(path)
  df_list <- list()
  N <- length(filenames)
  df_list <- vector("list", N)
  j <- 1
  for (i in c(1:N)) {
    varname <- paste('ri', stages[i], sep = "_")
    assign(varname, read.table(filenames[i], skip = 1, col.names = columns))

    print(colnames(get(varname)[c(1,column_id)]))
    df_list[[j]] <- get(varname)[c(1,column_id)]
    j <- j+1
  }
  total <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "gene_id"), df_list)
  print('OK')
  colnames(total) <-  c("gene_id", stages)
  return(total)
}


total_li <- AggregateCounts(path_li, filenames_li[c(3,4,5,6,7,14,8,9,10,11,12,13)],
samples_li[c(3,4,5,6,7,14,8,9,10,11,12,13)], 7, columns)

total_ri <- AggregateCounts(path_ri, filenames_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)],
samples_ri[c(11,12,1,2,3,4,5,6,7,8,9,10)], 7, columns)

#Calculating TE
total <- merge(total_ri[apply(total_ri >1, 1,all),], total_li, by =  "gene_id")
#total <- merge(total_ri, total_li, by =  "gene_id")
all_stages <- unique(sub( "*_rep1", "", sub("*_rep2", "", samples_li[c(3,4,5,6,7,14,8,9,10,11,12,13)])))
all_stages_TE <- paste(all_stages, '_TE', sep="")
for (i in all_stages){
  rep1x <- paste(i, '_rep1', '.x', sep="")
  rep1y <- paste(i, '_rep1', '.y', sep="")
  rep2x <- paste(i, '_rep2', '.x', sep="")
  rep2y <- paste(i, '_rep2', '.y', sep="")
  te <- paste(i, '_TE', sep="")
  total[te] <- (total[rep1y] / total[rep1x] + total[rep2y] / total[rep2x])/2

}
total_te <- total[, all_stages_TE]
for_violin_plot <- melt(total_te)
N <- ncol(total_te)-1
comparisons_1 <- vector("list", N)
for (i in c(1:N)){
  comparisons_1[[i]] <- c(colnames(total_te)[i], colnames(total_te)[i+1])
}

#Drawing boxplots and calculating statistics
#total_te1 <- total_te[apply(total_te < 6, 1,all),]
for_violin_plot1 <- melt(total_te)

stat.test <- compare_means(value~variable, data = for_violin_plot1, paired=F, method = "wilcox.test", p.adjust.method = "BH")
stat.test1 <- stat.test[c(1,6,10,13,15),]
ggboxplot(for_violin_plot, x = "variable", y = "value", outlier.shape = NA) +
 stat_pvalue_manual(stat.test1, label = "p.format", y.position = c(6,5.5,6,5.5,6)) +
 scale_x_discrete(name="", labels=all_stages) + theme(axis.text.x = element_text(size=14, angle=45)) +
 scale_y_continuous(name="TE", limits = quantile(for_violin_plot$value, c(0.1, 0.9)))

#Drawing plot of top-10 TE genes
top_10 <- round(nrow(total_te) /10 )
top_total <- total_te[c(1:top_10),]
N <- ncol(total_te)
for (i in c(1:N)){
  top_total[,i] <- sort(total_te[,i], decreasing = T)[c(1:top_10)]
}


for_violin_plot <- melt(top_total)
ggviolin(for_violin_plot, x = "variable", y = "value", outlier.shape = NA) +
 scale_x_discrete(name="", labels=all_stages) + theme(axis.text.x = element_text(size=14, angle=45)) +
 scale_y_continuous(name="TE", limits = quantile(for_violin_plot$value, c(0.1, 0.9)))+
 theme_classic()
ggsave('../top10.pdf')
head(aggregate(test_data$MII_rep1.x, by=list(Category=total$gene_id), FUN=sum))
