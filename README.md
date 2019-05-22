# Drug Repurposing Dengue

Repository with code for the analysis of Microarray and RNA-Seq data from Dengue cohorts.

**I - Introduction**

The workflow here described for the analysis of Microarray and RNA-Seq was performed in R. Many packages specific for Bioinformatic 
analysis were used such as "limma" or "DESeq2". Below are briefly summarized the workflows for the differential expression of Microarray
and RNA-Seq data as well as the data preparation to be used in CMap and all the plots created. 

Scripts with the code for each step here referred are in the correspondent folders.

**II - Differential Expression Evaluation in Array Data**

Analysis of Microarray data was done in three cohorts from Dengue cohorts: GSE18090, GSE38246 and GSE51808. Scripts used for each of the microarrays are in the correspondent folder. Packages used were "GEOquery", "limma", "genefilter" and the "hgu133plus2.db" database for the Affymetrix Human Genome U133 Plus 2.0 Array chip.

Below is described as an example the code for the analysis of GSE18090, scripts for GSE38246 and GSE51808 are similar although there are some slight changes due to the different chips used.

- Load of Necessary Packages
```R
library(GEOquery)
library(limma)
library(genefilter)
library(hgu133plus2.db)
```

- Load of Data
```
gset <- getGEO("GSE18090", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
```

- Remove of Duplicated Probes & Addition of Gene Annotations
```
gset@annotation = "hgu133plus2"
gset <- featureFilter(gset, require.entrez = T, remove.dupEntrez = T)
```

- Addition of Metadata and Selection of Samples
```
fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- "XXXXXXXX111111111100000000" #x-remove; 1-infected; 0-controls
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
```

- Differential Expression
```
exprs(gset) <- log2(exprs(gset))

sml <- paste("G", sml, sep="")
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=22000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="matrixDE18090.txt", row.names=F, sep="\t")
```

- Selection of Up and Down Regulated Genes
```
down <- tT[,"logFC"]<0
downreg <- tT[down,]
write.table(downreg$Gene.symbol, file="downreg18090.txt", row.names=F, sep="\t")

up <- tT[,"logFC"]>0
upreg <- tT[up,]
write.table(upreg$Gene.symbol, file="upreg18090.txt", row.names=F, sep="\t")
```
**III - Differential Expression Evaluation in RNASeq Data**

Analysis of RNA-Seq data was done in three samples from Dengue cohorts: Spleen, Hepatic and Encephalon. Scripts used for each of the RNA-Seq are in the correspondet folder. The package used for this analysis was "DESeq2".

Below is described as an example the code for the analysis of the Spleen data, scripts for Hepatic and Encephalic samples are similar.

- Load of Necessary Packages
```
library("DESeq2")
```

- Upload of Data and Selection of Samples
```
data <- read.delim("alinhamentoALL.txt", stringsAsFactors=FALSE, sep = "\t")

dataUnique <- matrix(as.integer(unlist(data[,2:dim(data)[2]])), nrow = dim(data)[1])

rownames(dataUnique) <- data[,1]

colnames(dataUnique) <- colnames(data[,2:dim(data)[2]])

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

baslocal <- c(2,6,12,18,24,28,31,37,44)

dataUnique <- dataUnique[,baslocal]

###remover individuos 875 e 900
dataUniquelocal <- dataUnique[,-c(5,6)]
```

- Defenition of Groups and Creation of DESeq Object
```
disease_state <- c("Infected","Infected","Infected","Infected",
                   "Control","Control","Control")

coldata <- data.frame(col.names = colnames(dataUniquelocal),
                      group = disease_state)

DESeq2Table <- DESeqDataSetFromMatrix(dataUniquelocal, colData = coldata, design = ~group)
```

- Differential Expression
``` 
colData(DESeq2Table)$condition <- factor(colData(DESeq2Table)$group)

DESeq2Table <- estimateSizeFactors(DESeq2Table)
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

DESeq2Table <- nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, contrast=c("group","Infected","Control"),
                     pAdjustMethod = "BH")

DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
plotMA(DESeq2Res)

DESeq2Res_sig <- DESeq2Res[ DESeq2Res$padj < 0.05, ]
DESeq2Res_sig <- DESeq2Res_sig[ order(DESeq2Res_sig$padj), ]
```

- Selection of Up and Down Regulated Genes
```
up <- DESeq2Res_sig[DESeq2Res_sig$log2FoldChange>0,]
down <- DESeq2Res_sig[DESeq2Res_sig$log2FoldChange<0,]

write.table(file="genesDEBA.txt", DESeq2Res, col.names=T, row.names=T, sep="\t", quote=F)
write.table(file="upregulatedBA.txt", row.names(up), sep = "\t")
write.table(file="downregulatedBA.txt", row.names(down), sep="\t")
```
**IV - Data Preparation for Drug Repurposing**

Selection of genes to be inputted in CMap. Scripts with code for the selection of genes with common patterns in two of the three analysed microarrays.

**V - Plot Creation**

Creation of plots for the demonstration of reesults. Volcano plots, barplots and Venn diagrams were created using the "ggplot2" package.
