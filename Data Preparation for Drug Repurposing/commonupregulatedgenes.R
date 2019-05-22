library(dplyr)

up38246 <- read.delim("upreg38246.txt", header = F)
up51808 <- read.delim("upreg51808.txt", header = F)
up18090 <- read.delim("upreg18090.txt", header = F)

common38246_51808 <- intersect(up38246$V1, up51808$V1)
length(common38246_51808)

common18090_38246 <- intersect(up18090$V1, up38246$V1)
length(common18090_38246)

common18090_51808 <- intersect(up18090$V1, up51808$V1)
length(common18090_51808)

common18090_51808_38246 <- intersect(common18090_38246, common38246_51808)
length(common18090_51808_38246)

upgenes <- append(common18090_38246, common18090_51808)
upgenes <- append(upgenes, common38246_51808)
upgenes <- unique(upgenes)