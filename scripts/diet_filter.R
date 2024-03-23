
### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Filter dietary peptides based on blastp alignment
library(dplyr)

## loading files
args = commandArgs(trailingOnly = TRUE)
genus_database <- read.table(args[1], header=TRUE, sep="\t", quote = "", check.names=FALSE, comment.char="") 
genus_database <- subset(genus_database, select=c("Protein", "Genus_TaxID"))

genus_blastout <- read.table(args[2], header=TRUE, sep="\t", quote = "", check.names=FALSE, comment.char="")
genus_blastout$Protein <- gsub("_[0-9]*$", "", genus_blastout$Query)

# combine files
genus_blastout <- merge( genus_blastout, genus_database, by = "Protein", all.x = TRUE)

# filter
count_peptides <- subset(genus_blastout, evalue<0.01) %>% group_by(Query) %>% arrange(evalue) %>% slice_max(with_ties = TRUE, order_by=bitscore) %>% summarize(Count=length(GenusTaxID[Genus_TaxID==GenusTaxID]))
count_peptides <- as.data.frame(count_peptides)
count_peptides <- subset(count_peptides, Count>0)
count_peptides$Protein <- gsub("_[0-9]*$", "", count_peptides$Query)

write.table(count_peptides, args[3], sep="\t", row.names=FALSE, quote=FALSE)