#Load Packages
library(Biobase)
library(GEOquery)
library(genefilter)
library(limma)

#Load GSE
gset <- getGEO("GSE51808", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gset@annotation = "hgu133plus2"

featureNames(gset) <- gsub("_PM","",featureNames(gset))
controls <- grep("AFFX",featureNames(gset))
gset <- gset[-controls,]
gset <- featureFilter(gset, require.entrez = T, remove.dupEntrez = T)

# Add metadata 
fvarLabels(gset) <- make.names(fvarLabels(gset))

#Sample selection: x-remove; 1-infected; 0-controls
gsms <- "XXXX1X1XX1X111XXXX1X11XX1XXXXXXXXXXXXXXXXXXXXXX000000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=20000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="genesDE51808.txt", row.names=F, sep="\t")

###Select up & down regulated genes
down <- tT[,"logFC"]<0
downreg <- tT[down,]
write.table(downreg$Gene.symbol, file="downreg51808.txt", row.names=F, sep="\t")

up <- tT[,"logFC"]>0
upreg <- tT[up,]
write.table(upreg$Gene.symbol, file="upreg51808.txt", row.names=F, sep="\t")