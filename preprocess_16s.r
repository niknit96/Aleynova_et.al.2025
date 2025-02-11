library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

dir = as.character(read.table("dir.txt"))
setwd(dir)



asv_meta = read_qza("./16s/dada2/FeatureTable[Frequency]_16s.qza") # Table of metagenome data reads for each ASV by samples

asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)

# Obtain ASVs (16s)
seq_16s = read_qza("./16s/dada2/FeatureData[Sequence]_16s.qza")
seq_16s = as.data.frame(seq_16s$data)
seq_16s$ASV = row.names(seq_16s)
colnames(seq_16s) = c("Sequence", "ASV")
seq_16s = as_tibble(seq_16s)
#


tax_meta = read_qza("./16s/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_16s.qza") # Taxonomy table for metagenome data
tax_meta = parse_taxonomy(tax_meta$data, trim_extra=FALSE)
tax_meta[is.na(tax_meta)] <- "uncultured"

# Swap values ​​from Swap_unidentified in tax_meta to previous taxon level
Swap_unidentified = c("uncultured", "unidentified", "metagenome", "bacteriap25", "Unknown", "Incertae_Sedis")
for(Swap in Swap_unidentified) {
    for(j in 1:7) {
        for (i in 1:nrow(tax_meta)) {
            if(grepl(Swap, tax_meta[i,j])) {
                tax_meta[i,j] = tax_meta[i,j-1] }
                else { tax_meta[i,j] = tax_meta[i,j]
            }       
        }
    }
}


#

# Filtering metagenome data from non-significant taxa
tax_meta = filter(tax_meta, 
Phylum != "d__Bacteria" & 
Phylum != "Unassigned" &
Kingdom != "d__Eukaryota" & 
Phylum != "k__Fungi" & 
Kingdom != "d__Archaea" &
Kingdom != "k__Viridiplantae" &
Genus != "g__Mitochondria" &
Genus != "g__Chloroplast")
#

tax_suppl = tax_meta
tax_suppl$ASV = row.names(tax_suppl)
tax_suppl = left_join(tax_suppl, seq_16s, by = "ASV")
tax_suppl = tax_suppl[,c("Sequence", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")]
colnames(tax_suppl) = c("Amplicon sequence variant", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")



tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")



sampledata_16s = read.table(file="./Metagenome metadata/sampledata_16s_meta.txt", header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome
row.names(sampledata_16s) = sampledata_16s$SRA


# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_16s)

tax_16s = tax_table(as.matrix(tax_meta))
asv_16s = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq_16s = phyloseq(asv_16s, tax_16s, sampledata_meta)


# phyloseq object for sow (16S)
sampledata_16s_sow = read.table(file="./16s/Culture-dependent data/sampledata_16s_sow.txt", header = TRUE, row.names = 1, colClasses = "character", sep="\t")
asv_16s_sow = read.table(file="./16s/Culture-dependent data/asv_16s_sow.txt", header = TRUE, row.names = 1, sep="\t")
tax_16s_sow = read.table(file="./16s/Culture-dependent data/tax_16s_sow.txt", header = TRUE, row.names = 1, sep="\t")
sampledata_16s_sow = sample_data(sampledata_16s_sow)
tax_16s_sow = tax_table(as.matrix(tax_16s_sow))
asv_16s_sow = otu_table(asv_16s_sow, taxa_are_rows = TRUE)
physeq_16s_sow = phyloseq(asv_16s_sow, tax_16s_sow, sampledata_16s_sow)


filtered_data = as.data.frame(colSums(otu_table(physeq_16s)))
filtered_data = as.data.frame(t(t(filtered_data)))
colnames(filtered_data) = "Sequences after filtration used in analysis"
filtered_data$SRA = row.names(filtered_data)
raw_data = read_qza("./16s/dada2/16s-stats-dada2.qza")$data
raw_data$SRA = row.names(raw_data)


data_16s = left_join(raw_data,filtered_data, by="SRA")

data_16s = data_16s[,c("SRA", "input", "Sequences after filtration used in analysis")]
colnames(data_16s) = c("SRA", "Raw paired-end reads", "Sequences after filtration used in analysis")
data_16s$Sum_raw = sum(data_16s[,"Raw paired-end reads"])
data_16s$Mean_raw = mean(data_16s[,"Raw paired-end reads"])
data_16s$Median_raw = median(data_16s[,"Raw paired-end reads"])
data_16s$Sum_filtered = sum(data_16s[,"Sequences after filtration used in analysis"])
data_16s$Mean_filtered = mean(data_16s[,"Sequences after filtration used in analysis"])
data_16s$Median_filtered = median(data_16s[,"Sequences after filtration used in analysis"])

sampledata_16s = left_join(sampledata_16s, data_16s, by="SRA")

write.table(sampledata_16s, file="Sampledata_16s_with_summary.txt", quote=FALSE, sep="\t")

save.image(file='./16s/physeq_16s.RData')
