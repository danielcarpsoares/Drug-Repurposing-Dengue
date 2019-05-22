library(DESeq2)

###Data upload
data <- read.delim("alinhamentoALL.txt", stringsAsFactors=FALSE, sep = "\t")

###Sets dataUnique has a matrix of intergers
dataUnique <- matrix(as.integer(unlist(data[,2:dim(data)[2]])), nrow = dim(data)[1])

###Column names start in the second column
rownames(dataUnique) <- data[,1]
colnames(dataUnique) <- colnames(data[,2:dim(data)[2]])

###Vectors with information for columns
alinhamento <- c("global", "local", "global", "local", "global", "local","global", "local",
                 "global", "local", "global", "local", "global", "local", "global", "local",
                 "global", "local", "global", "local", "global", "local", "global", "local",
                 "global", "local", "global", "local", "global", "local", "local", "global",
                 "local", "global", "local", "global", "local", "global", "local", "global",
                 "local", "local", "global", "local", "global", "local", "global", "global" )

tissue <- c("Ba", "Ba", "Hep", "Hep", "Ba", "Ba", "Enc", "Enc", "Hep", "Hep", "Ba", "Ba",
            "Enc", "Enc", "Hep", "Hep", "Ba", "Ba", "Enc", "Enc", "Hep", "Hep", "Ba", "Ba",
            "Hep", "Hep", "Ba", "Ba", "Hep", "Hep", "Ba", "Enc", "Enc", "Hep", "Hep", "Ba", 
            "Ba","Enc", "Enc", "Hep", "Hep", "Hep","Ba", "Ba", "Enc", "Enc", "Hep", "Ba")

###Pick only local aligned samples for Liver
heplocal <- c(4,10,16,22,35,41,42)

dataUniquehep <- dataUnique[,heplocal]

###Set conditions for DE
disease_state <- c("Infected","Infected","Infected","Infected",
                   "Control","Control","Control")

coldata <- data.frame(col.names = colnames(dataUniquehep),
                      group = disease_state)

DESeq2Table <- DESeqDataSetFromMatrix(dataUniquehep, colData = coldata, design = ~group)
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$group)

# Estimate size factors
DESeq2Table <- estimateSizeFactors(DESeq2Table)
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <- nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, contrast=c("group","Infected","Control"),
                     pAdjustMethod = "BH")

# remove missing pvalue data
DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]

# PLot results
plotMA(DESeq2Res)

# Ordering 
DESeq2Res_sig <- DESeq2Res[ DESeq2Res$pvalue < 0.05, ]
DESeq2Res_sig <- DESeq2Res_sig[ order(DESeq2Res_sig$pvalue), ]

up <- DESeq2Res_sig[DESeq2Res_sig$log2FoldChange>0,]
down <- DESeq2Res_sig[DESeq2Res_sig$log2FoldChange<0,]

write.table(file="denesDEHep.txt", DESeq2Res, col.names=T, row.names=T, sep="\t", quote=F)
write.table(file="upregHep.txt", row.names(up), sep = "\t")
write.table(file="downregHep.txt", row.names(down), sep="\t")