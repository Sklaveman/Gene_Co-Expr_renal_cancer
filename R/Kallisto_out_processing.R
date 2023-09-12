library(AnnotationHub)
library(ensembldb)
library(tximport)
library(dplyr)
library(readr)
library(stringr)
library(dplyr)
library(tibble)


#Connect to Connect to AnnotationHub
#Return the Ensembl EnsDb information for Homo sapiens
ah <- AnnotationHub()
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

#Extract the latest Ensembl release using the AnnotationHub ID to subset the object
human_ens <- human_ens[["AH109606"]]

#Extract info about genes 
genedb <- genes(human_ens, return.type = "data.frame")
genedb$symbol <- ifelse(genedb$symbol == "", genedb$gene_id, genedb$symbol)
genedb2 <- dplyr::select(genedb, gene_id, symbol)

#Extract info about transcripts 
txdb <- transcripts(human_ens, return.type = "data.frame")
txdb2 <- dplyr::select(txdb, tx_id_version, gene_id)

#Join transcripts, genes IDs and symbols into one table
#Transcript ID must be the first column!
tx2gene <- right_join(txdb2, genedb2, by='gene_id')
tx2gene$symbol <- ifelse(tx2gene$symbol == "", tx2gene$gene_id, tx2gene$symbol)

#Optional conversion of normalized TPM counts into non-normalized counts
##Read in metadata file 
metadata_mod2 <- read_csv("Immune_project/data/metadata/metadata_mod2.csv",
                          col_types = cols(subject_id = col_character()))

#Obtain absolute paths to kallisto output abundance tables
samples <- list.dirs(path = "kal_quant_out_new", full.names = T)[2:264]
files <- file.path(samples, "abundance.tsv")
names(files) <- str_replace(samples, "kal_quant_out_new/", "")

# Estimate counts from abundances
txi <- tximport(files, 
                type="kallisto", 
                tx2gene=tx2gene[,c("tx_id_version", "gene_id")],
                countsFromAbundance="lengthScaledTPM")

#Filter out pseudogene IDs
genedb_no_pseudo <- genedb %>% filter(!grepl('pseudogene', gene_biotype, ignore.case ="True") )

#Construct the data frame with abundances
abud <- data.frame(txi$abundance)
colnames(abud) <- str_replace(colnames(abud), "kal_quant_out_new.", "")
abud <- abud %>% rownames_to_column(., var = "gene_id")
abud <- right_join(genedb_no_pseudo[,c("gene_id", "symbol")], abud, by="gene_id")
abud <- abud %>% filter(!is.na(abud$symbol))
abud <- abud[,2:265]

#Construct the data frame with counts
raw_counts <- data.frame(round(txi$counts))
colnames(raw_counts) <- str_replace(colnames(raw_counts), "kal_quant_out_new.", "")
raw_counts[(rownames(raw_counts) == "ENSG00000171163" | rownames(raw_counts) == "ENSG00000185220"), 1:5]
raw_counts <- raw_counts %>% rownames_to_column(., var = "gene_id")
raw_counts <- right_join(genedb_no_pseudo[,c("gene_id", "symbol")], raw_counts, by="gene_id")
raw_counts <- raw_counts %>% filter(!is.na(raw_counts$symbol))
raw_counts <- raw_counts[,2:265]

#Save data frames of counts and abundances
write.table(abud, "Immune_project/data/patients_abundances.csv", quote=F, sep=',', na="NaN", row.names=F)
write.table(raw_counts, "Immune_project/data/patients_counts.csv", quote=F, sep=',', na="NaN", row.names=F)

#Save gene annotation data frame in workspace file
save(genedb_no_pseudo, file = "genedb_no_pseudo.RData")
