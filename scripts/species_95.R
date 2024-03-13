
### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas

### Bacterial species abundances < 95%

## loading libraries
library(dplyr)
library(tidyr)
library(purrr)
library(readr)


## loading bacterial species taxIDs
species_taxIDs <- read.table("automatic/bacteria-species.detailed.txt", header=FALSE, sep="\t", quote = "", check.names=FALSE, comment.char="")
colnames(species_taxIDs) <- c("TaxID", "Species", "Level")
genbank_genomes <- read.table("assembly_summary.genbank.txt", header=TRUE, sep="\t", quote = "", check.names=FALSE, comment.char="", skip=1)
refseq_genomes <- read.table("assembly_summary.refseq.txt", header=TRUE, sep="\t", quote = "", check.names=FALSE, comment.char="", skip=1)

## select representative genomes
taxIDs <- unique(refseq_genomes$taxid)
for (i in 1:length(taxIDs)){
	refseq_taxID <- subset( refseq_genomes, taxid == taxIDs[i] )
	refseq_taxID %>% arrange(refseq_category, total_gene_count, genome_rep) %>% 
}

## reading bracken outputs
customized_read_tsv <- function(file){
	id <- gsub('bracken/(.*)/bracken', '\\1', file)
    read_tsv(file) %>%
        mutate(sampleName = id)
}

bracken <- list.files(path="bracken",pattern='bracken$',recursive=TRUE, full.names=TRUE) %>%
    lapply(customized_read_tsv) %>% 
    reduce(bind_rows) %>% 
    select(sampleName, taxonomy_id, new_est_reads) %>% 
    pivot_wider(names_from = sampleName, values_from = new_est_reads) 

bracken <- as.data.frame(bracken)
bracken[is.na(bracken)] <- 0

## selecting bacterial species
bracken <- subset(bracken, taxonomy_id %in% species_taxIDs$TaxID)
row.names(bracken) <- bracken$taxonomy_id
bracken <- subset(bracken, select=-c(taxonomy_id))

counts <- as.data.frame(rowSums(bracken) )
colnames(counts) <- c("reads")

## select <95
counts <- merge(counts, species_taxIDs, by.x="row.names", by.y="TaxID")
counts$Percentage <- 100 * counts$reads / sum(counts$reads)
counts <- counts %>% dplyr::arrange(-Percentage)
counts[,"Cum_Percentage"] <- cumsum(counts$Percentage)
counts$number <- row.names(counts)
counts$number <- as.numeric(counts$number)
colnames(counts)[1] <- "TaxID"

candidates <- subset(counts, Cum_Percentage<=95)
write.table(candidates, "bacteria.95.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## selecting assembly_accession
genomes <- subset(species_genomes, as.character(taxid) %in% as.character(candidates$TaxID) )
final <- genomes %>% dplyr::distinct(ncbi_species_taxid, .keep_all = TRUE)
write.table(final$ncbi_genbank_assembly_accession, "genomes.all.path", row.names=FALSE, sep="\t", quote=FALSE)

