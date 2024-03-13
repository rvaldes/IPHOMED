
### Integrated Proteomics of HOst-MicrobiomE-Diet (IPHOMED)
### Rafael Valdés-Mas (Elinav Lab)

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
#refseq_genomes <- read.table("assembly_summary.refseq.txt", header=TRUE, sep="\t", quote = "", check.names=FALSE, comment.char="", skip=1)

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
genomes <- subset(genbank_genomes, as.character(species_taxid) %in% as.character(candidates$TaxID) )
taxIDs <- unique(genomes$species_taxid)
reference_genomes <- genomes[FALSE, ]
for (i in 1:length(taxIDs)){
    message(i)
	genbank_taxID <- subset(genomes, species_taxid == taxIDs[i] )
	if (length(subset(genbank_taxID, refseq_category=="reference genome" )$species_taxid) > 0){
        genbank_taxID <- subset(genbank_taxID, refseq_category=="reference genome" ) %>% arrange(-as.numeric(protein_coding_gene_count))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }else if (length(subset(genbank_taxID, refseq_category=="representative genome" )$species_taxid) > 0){
        genbank_taxID <- subset(genbank_taxID, refseq_category=="representative genome" ) %>% arrange(-as.numeric(protein_coding_gene_count))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }else if (length( subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0 & assembly_level=="Complete Genome")$species_taxid ) > 0){
        genbank_taxID <- subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0 & assembly_level=="Complete Genome") %>% arrange(-as.numeric(protein_coding_gene_count))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }else if (length( subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0 & genome_rep=="Full")$species_taxid ) > 0){
        genbank_taxID <-  subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0 & genome_rep=="Full") %>% arrange(-as.numeric(protein_coding_gene_count))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }else if (length( subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0)$species_taxid ) > 0){
        genbank_taxID <- subset(genbank_taxID, as.numeric(protein_coding_gene_count) > 0) %>% arrange(-as.numeric(protein_coding_gene_count))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }else{
        genbank_taxID <- genbank_taxID %>% arrange(-as.numeric(genome_size))
        reference_genomes <- rbind(reference_genomes, genbank_taxID[1,])
    }
}

write.table(reference_genomes, "genomes.info.tsv", row.names=FALSE, sep="\t", quote=FALSE)
write.table(reference_genomes$assembly_accession, "genome_assembly.tsv", row.names=FALSE, col.names = FALSE, sep="\t", quote=FALSE)