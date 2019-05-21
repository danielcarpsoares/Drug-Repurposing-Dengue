# Drug Repurposing Dengue

Repository with code for the analysis of Microarray and RNA-Seq data from Dengue cohorts

**I - Introduction**

The workflow here described for the analysis of Microarray and RNA-Seq was performed in R. Many packages specific for Bioinformatic 
analysis were used such as "limma" or "DESeq2". Below are briefly summarized the workflows for the differential expression of Microarray
and RNA-Seq data as well as the data preparation to be used in CMap and all the plots created. 

Scripts with the code for each step here referred are in the correspondent folders.

**II - Differential Expression Evaluation in Array Data**

Analysis of Microarray data was done in three cohorts from Dengue cohorts: GSE18090, GSE38246 and GSE51808. Scripts used for each of the microarrays are in the correspondent folder. Packages used were "GEOquery", "limma", "genefilter" and the "hgu133plus2.db" database for the Affymetrix Human Genome U133 Plus 2.0 Array chip.

**III - Differential Expression Evaluation in RNASeq Data**

Analysis of RNA-Seq data was done in three samples from Dengue cohorts: Spleen, Hepatic and Encephalon. Scripts used for each of the RNA-Seq are in the correspondet folder. The package used for this analysis was "DESeq2".

**IV - Data Preparation for Drug Repurposing**

**V - Plot Creation**
