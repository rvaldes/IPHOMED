### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

### Selecting species with a callability > 70%

#reading files
args = commandArgs(trailingOnly=TRUE)
coverage <- read.table(args[1], header=FALSE, sep="\t")
colnames(coverage) <- c("Assembly", "Total", "Count")

# processing callability
coverage$Per <- 100 * coverage$Count / coverage$Total

# plotting distribution
pdf("Species_coverage.pdf")
hist(coverage$Per, breaks=50)
dev.off()

# filtering
coverage <- subset(coverage, Per >= 70)
coverage$Assembly <- gsub("[.]1_.*", ".1", coverage$Assembly)
coverage$Assembly <- gsub("[.]2_.*", ".2", coverage$Assembly)
coverage$Assembly <- gsub("[.]3_.*", ".3", coverage$Assembly)
coverage$Assembly <- gsub("[.]4_.*", ".4", coverage$Assembly)
write.table(coverage$Assembly, args[2], row.names=FALSE, sep="\t", col.names=FALSE, quote=FALSE) 