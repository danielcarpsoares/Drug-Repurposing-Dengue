library(dplyr)

down38246 <- read.delim("genesdown38246.txt", header = F)
down51808 <- read.delim("genesdown51808.txt", header = F)
down18090 <- read.delim("genesdown18090.txt", header = F)

common38246_51808 <- intersect(down38246$V1, down51808$V1)
length(common38246_51808)

common18090_38246 <- intersect(down18090$V1, down38246$V1)
length(common18090_38246)

common18090_51808 <- intersect(down18090$V1, down51808$V1)
length(common18090_51808)

common18090_51808_38246 <- intersect(common18090_38246, common38246_51808)
length(common18090_51808_38246)

downgenes <- append(common18090_38246, common18090_51808)
downgenes <- append(downgenes, common38246_51808)
downgenes <- unique(downgenes)

write.table(downgenes, file="downgenes.txt", row.names=F, sep="\t")