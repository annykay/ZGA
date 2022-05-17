library(DESeq2)
path_trt <- "D:/study/OmicsData/project/treatedRNA/"

#Reading filenames
setwd(path_trt)
filenames_trt  <- list.files(pattern="*.txt", full.names=TRUE)
samples_trt <- sub(".*RNA_Zygote_*", "", sub("*.quant.genes.*", "", filenames_trt))

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

minorZGA <- AggregateCounts(path_trt, filenames_trt, samples_trt, 7, columns)
Averaging <- function(x, c1){
  return(mean(x[c1]))
}
minorZGA$avAMA <- apply(minorZGA[,-1], 1, Averaging, c1 = c(1,2))
minorZGA$avControl <- apply(minorZGA[,-1], 1, Averaging, c1 = c(3,4))
#minorZGA[, c("avAMA", "avControl")] <- apply(minorZGA[,c(2,3,4,5)], 1, Averaging, c1 = c(1,2), c2  = c(3,4))
minorZGA <- minorZGA[minorZGA$avControl >=1, ]
minorZGA$FC <- minorZGA$avAMA / minorZGA$avControl

coldata <- data.frame(rownames = colnames(minorZGA)[c(2:5)], condition = c(1, 1, 0, 0))
rownames(coldata) <- coldata$rownames
coldata$replicate <- c(1,2,1,2) coldata <- coldata[, -1]
head(coldata)
coldata$condition <- factor(coldata$condition)
counts <- apply(minorZGA[,c(2,3,4,5)], 2, as.integer)
rownames(counts) <- minorZGA[,1]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

minorZGA$log2FC <- log2(minorZGA$FC)
minorZGA$FC_norm <- (minorZGA$FC - mean(minorZGA$FC))/ sd(minorZGA$FC )

minorZGA$p.val <- 1 - pnorm(abs(minorZGA$FC_norm))
minorZGA$p_adj <- p.adjust(minorZGA$p.val, method = 'BH')
res_minorZGA <- minorZGA[minorZGA$p_adj < 0.05, ]
total_miZGA <- total[total$gene_id %in% res_minorZGA$gene_id, ]

total_miZGA <- total[total$gene_id %in% res_minorZGA$gene_id, ]
total_te <- total[total$gene_id %in% res_minorZGA$gene_id, all_stages_TE]
for_violin_plot1 <- melt(total_te)
stat.test <- compare_means(value~variable, data = for_violin_plot1, paired=F, method = "wilcox.test", p.adjust.method = "BH")
stat.test1 <- stat.test[c(1,6,10,13,15),]
ggboxplot(for_violin_plot, x = "variable", y = "value", outlier.shape = NA) +
 stat_pvalue_manual(stat.test1, label = "p.format", y.position = c(6,5.5,6,5.5,6)) +
 scale_x_discrete(name="", labels=all_stages) + theme(axis.text.x = element_text(size=14, angle=45)) +
 scale_y_continuous(name="TE", limits = quantile(for_violin_plot$value, c(0.1, 0.9)))

plotBoxes <- function(data, name, all_stages_TE, total){
 coldata <- data.frame(rownames = colnames(data)[c(2,3,4,5)], condition = c(0, 0, 1, 1))
 rownames(coldata) <- coldata$rownames
 coldata$replicate <- c(1,2,1,2)
 coldata <- coldata[, -1]
 coldata$condition <- factor(coldata$condition)
 counts <- apply(data[,c(2,3,4,5)], 2, as.integer)
 rownames(counts) <- data[,1]
 counts <- counts[counts[, 3] >=1, ]
 counts <- counts[counts[, 4] >=1, ]

 dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)
 dds <- DESeq(dds)
 res <- results(dds)

 resZGA <- res[abs(res$log2FoldChange) > log2(5), ]
 resZGA <- resZGA[!is.na(resZGA$padj), ]
 resZGA <- resZGA[abs(resZGA$padj) <0.01, ]

 total_te <- total[total$gene_id %in% rownames(resZGA), all_stages_TE]
 for_violin_plot <- melt(total_te)
 stat.test <- compare_means(value~variable, data = for_violin_plot, paired=F, method = "wilcox.test", p.adjust.method = "BH")
 stat.test1 <- stat.test[c(1,6,10,13,15),]
 ggboxplot(for_violin_plot, x = "variable", y = "value", outlier.shape = NA) +
  stat_pvalue_manual(stat.test1, label = "p.format", y.position = c(4,3.5,4,3.5,4)) +
  scale_x_discrete(name="", labels=all_stages) + theme(axis.text.x = element_text(size=14, angle=45)) +
  scale_y_continuous(name="TE", limits = quantile(for_violin_plot$value, c(0.1, 0.9)))
 ggsave(name)
}
plotBoxes(total_ri[,c(1,4,5,6,7)], 'majorZGA.pdf', all_stages_TE, total)
plotBoxes(total_ri[,c(1,6,7,8,9)], 'twoCell.pdf', all_stages_TE, total)
plotBoxes(total_ri[,c(1,8,9,10,11)], 'fourCell.pdf', all_stages_TE, total)
total_ri$MII <- (total_ri[,2] +total_ri[,3])/2
m2genes <- total_ri[total_ri$MII>20, 1]
total_te <- total[total$gene_id %in% m2genes, all_stages_TE]
