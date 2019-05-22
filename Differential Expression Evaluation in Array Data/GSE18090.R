################################################################
###   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(genefilter)
library(hgu133plus2.db)

###Load series and platform data from GEO
gset <- getGEO("GSE18090", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

###Remove duplicated probes
gset@annotation = "hgu133plus2"
gset <- featureFilter(gset, require.entrez = T, remove.dupEntrez = T)

###Add metadata to dataset 
fvarLabels(gset) <- make.names(fvarLabels(gset))

###Group names for all samples: x-remove; 1-infected; 0-controls
gsms <- "XXXXXXXX111111111100000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

###Eliminate samples marked as "x"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

###Log2 transform
exprs(gset) <- log2(exprs(gset))

###Differential Expression
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
write.table(tT, file="genes18090deftest.txt", row.names=F, sep="\t")

#Select up and down regulated genes
down <- tT[,"logFC"]<0
downreg <- tT[down,]
write.table(downreg$Gene.symbol, file="downreg18090.txt", row.names=F, sep="\t")

up <- tT[,"logFC"]>0
upreg <- tT[up,]
write.table(upreg$Gene.symbol, file="upreg18090.txt", row.names=F, sep="\t")