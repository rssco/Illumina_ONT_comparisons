# NovaSeq read processing 

## Table of contents
[1. Bacteria](#novaseq_bacteria)  
- [Cutadapt](#novaseq_bacteria_cutadapt)  
- [Split 16S/18S](#novaseq_split)  
- [DADA2 16S](#novaseq_16S_dada2)  
- [DADA2 18S](#novaseq_18S_dada2)  
- [Phyloseq with 16S and 18S](#ps_16S_18S) 
- [Remove plastid reads](#novaseq_remove_plastid) 
- [Decontamination with Microdecon](#novaseq_microdecon) 
- [Remove Homopolymer](#novaseq_homopolymer) 
- [Remove chimeras vsearch](#novaseq_vsearch) 
- [Identity Ascophyllum ASVs](#novaseq_asco) 
- [Final novaseq ps object for sequencing comparison](#novaseq_ps_final) 
  
[2. Fungi](#novaseq_fungi)  
- [Cutadapt](#novaseq_fungi_cutadapt)  
- [DADA2](#novaseq_fungi_dada2)  


## 1. Bacteria <a name="novaseq_bacteria"></a>

### Cutadapt <a name="novaseq_bacteria_cutadapt" ></a>
```bash
#!/usr/bin/env bash
#SBATCH --job-name=cutadapt
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1 
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load cutadapt

PRIMER_F=GTGYCAGCMGCCGCGGTAA
PRIMER_R=CCGYCAATTYMTTTRAGTTT
MIN_LENGTH=220
output=03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_FILTERED/02_CUTADAPT

for sample in $(cat sample.txt); \
do \
   cutadapt  -m ${MIN_LENGTH} --discard-untrimmed --report=minimal -g ${PRIMER_F} -G ${PRIMER_R} ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz -o ${output}/${s
ample}_R1_cutadapt.fastq.gz -p ${output}/${sample}_R2_cutadapt.fastq.gz; done
```

### Split 16S/18S <a name="novaseq_split"></a>
https://astrobiomike.github.io/amplicon/16S_and_18S_mixed

```bash
#!/usr/bin/env bash
#SBATCH --job-name=magicblast
#SBATCH -p fast
#SBATCH --mem=7G
#SBATCH --cpus-per-task 4
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL
#SBATCH --array 1-129


module load magicblast/1.5.0

query=$(ls 03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/01_CUTADAPT/*R1* | awk "NR==$SLURM_ARRAY_TASK_ID")
database=03_NOVASEQ_METAB/02_DATABASES
out=03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/01_MBLAST_OUT


magicblast -db ${database}/pr2-magicblast-db -query ${query} -query_mate  ${query%_R1_cutadapt.fastq.gz}_R2_cutadapt.fastq.gz -infmt fastq -out ${out}/$(basename ${qu
ery%_R1_cutadapt.fastq.gz}_mblast_out.txt) -outfmt tabular -num_threads $SLURM_CPUS_PER_TASK -splice F -no_unaligned
```
Filtering magic-Blast output

```bash
for i in *txt; do cut -f 1,2,3,7,8,16 $i | sed '1d' | sed '1d' | sed 's/# Fields: //' | tr " " "_" | awk -F $'\t' ' BEGIN { OFS=FS } NR==1 { $7="%_query_aln"; print $0 } NR>1 { print $0, ($5-$4)/$6*100 } ' > ${i%_mblast_out.txt}_mblast_out_mod.txt; done

for i in *; awk ' $3 > 90 && $7 > 35 ' $i | cut -f 1 | uniq -d > ${i%_mblast_out_mod.txt)_18S_headers.txt; done
```

Splitting fastq files intro 16S/18S
```bash
#!/usr/bin/env bash
#SBATCH --job-name=split
#SBATCH -p fast
#SBATCH --mem=5G
#SBATCH --cpus-per-task 2
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL
#SBATCH --array 2-2

input=$(ls /shared/projects/seabioz/finalresult/03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/01_CUTADAPT/*R1_cutadapt* | awk "NR==$SLURM_ARRAY_TASK_ID")

python split_16S_18S_reads.py -f ${input} -r ${input%_R1_cutadapt.fastq.gz}_R2_cutadapt.fastq.gz -E ${input%_R1_cutadapt.fastq.gz}_18S_headers.txt
```

### DADA2 16S <a name="novaseq_16S_dada2"></a>
```r
library(dada2); packageVersion("dada2")


path="03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/03_SPLIT_16S/"
files <- list.files(path)

fnFs_16S <- sort(list.files(path, pattern="_R1_cutadapt.fastq.gz", full.names = TRUE))
fnRs_16S <- sort(list.files(path, pattern="_R2_cutadapt.fastq.gz", full.names = TRUE))

sample.names_16S <- sapply(strsplit(basename(fnFs_16S),  "_cutadapt"), `[`, 1)

filtFs_16S <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_16S, "_R1_filt_dada2.fastq.gz"))
filtRs_16S <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_16S, "_R2_filt_dada2.fastq.gz"))
names(filtFs_16S) <- sample.names_16S
names(filtRs_16S) <- sample.names_16S

out_16S <- filterAndTrim(fnFs_16S, filtFs_16S, fnRs_16S, filtRs_16S, maxN=0, maxEE=c(2,2), truncLen = c(231,230), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=
TRUE) 
out_16S


# SAVE 
saveRDS(out_16S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/01_out_filt_split_16
S.rds")


errF_16S <- learnErrors(filtFs_16S, multithread=TRUE)
errR_16S <- learnErrors(filtRs_16S, multithread=TRUE)
dadaFs_16S <- dada(filtFs_16S, err=errF_16S, multithread=TRUE)
dadaRs_16S <- dada(filtRs_16S, err=errR_16S, multithread=TRUE)
mergers_16S <- mergePairs(dadaFs_16S, filtFs_16S, dadaRs_16S, filtRs_16S, verbose=TRUE)
seqtab_16S <- makeSequenceTable(mergers_16S)
seqtab.nochim_16S <- removeBimeraDenovo(seqtab_16S, method="consensus", multithread=TRUE, verbose=TRUE)
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/06_seqtabnochim_16S.RData"
)
getN_16S <- function(x) sum(getUniques(x))
track_16S <- cbind(out_16S, sapply(dadaFs_16S, getN_16S), sapply(dadaRs_16S, getN_16S), sapply(mergers_16S, getN_16S), rowSums(seqtab.nochim_16S))
colnames(track_16S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_16S) <- sample.names_16S
track_16S

save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/07_tracks_16S.RData")

taxa_16S <- assignTaxonomy(seqtab.nochim_16S, "03_NOVASEQ_METAB/02_DATABASES/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa.print_16S <- taxa_16S # Removing sequence rownames for display only
rownames(taxa.print_16S) <- NULL

# SAVE 
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/08_assigntaxo_16S.RData")
saveRDS(taxa_16S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/09_taxa_16S.rds")
saveRDS(seqtab.nochim_16S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/10_seqtab.n
ochim_16S.rds")
```

### DADA2 18S <a name="novaseq_18S_dada2"></a>

```r
library(dada2); packageVersion("dada2")


path="03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/04_SPLIT_18S/"
files <- list.files(path)

fnFs_18S <- sort(list.files(path, pattern="_R1_cutadapt.fastq.gz", full.names = TRUE))
fnRs_18S <- sort(list.files(path, pattern="_R2_cutadapt.fastq.gz", full.names = TRUE))

sample.names_18S <- sapply(strsplit(basename(fnFs_18S), "_cutadapt"), `[`, 1)

filtFs_18S <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_18S, "_R1_filt_dada2.fastq.gz"))
filtRs_18S <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_18S, "_R2_filt_dada2.fastq.gz"))
names(filtFs_18S) <- sample.names_18S
names(filtRs_18S) <- sample.names_18S

out_18S <- filterAndTrim(fnFs_18S, filtFs_18S, fnRs_18S, filtRs_18S, maxN=0, maxEE=c(2,2), truncLen = c(231,230), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=T
RUE) 
out_18S

save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/12_out_filt_18S.RData")

errF_18S <- learnErrors(filtFs_18S, multithread=TRUE)
errR_18S <- learnErrors(filtRs_18S, multithread=TRUE)
dadaFs_18S <- dada(filtFs_18S, err=errF_18S, multithread=TRUE)
dadaRs_18S<- dada(filtRs_18S, err=errR_18S, multithread=TRUE)

mergers_18S <- mergePairs(dadaFs_18S, filtFs_18S, dadaRs_18S, filtRs_18S, verbose=TRUE, justConcatenate = TRUE)
seqtab_18S <- makeSequenceTable(mergers_18S)
seqtab.nochim_18S <- removeBimeraDenovo(seqtab_18S, method="consensus", multithread=TRUE, verbose=TRUE)

save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/08_seqtabnochim_18S.RData")

getN_18S <- function(x) sum(getUniques(x))
track_18S <- cbind(out_18S, sapply(dadaFs_18S, getN_18S), sapply(dadaRs_18S, getN_18S), sapply(mergers_18S, getN_18S), rowSums(seqtab.nochim_18S))
colnames(track_18S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_18S) <- sample.names_18S
track_18S

taxa_18S <- assignTaxonomy(seqtab.nochim_18S, "03_NOVASEQ_METAB/02_DATABASES/pr2_version_5.0.0_SSU_dada2.fasta.gz",
                       taxLevels = c("Kingdom","Supergroup","Division", "Subdivision", "Class","Order","Family","Genus","Species"), multithread=TRUE)

taxa.print_18S <- taxa_18S 
rownames(taxa.print_18S) <- NULL

# SAVE 
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/14_dada_pipeline_18S.RData"
)
saveRDS(taxa_18S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/15_taxa_18S.rds")
saveRDS(seqtab.nochim_18S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/16_seqtab.n
ochim_18S.rds")
```

### Phyloseq with 16S and 18S <a name="ps_16S_18S"></a>
```r
taxa_16S <- readRDS("00_PHYLOSEQ_OBJECTS/00_taxa_16S.rds")
taxa_18S <- readRDS("00_PHYLOSEQ_OBJECTS/00_taxa_18S.rds")
seqtab.nochim_16S <- readRDS("00_PHYLOSEQ_OBJECTS/00_seqtab.nochim_16S.rds")
seqtab.nochim_18S <- readRDS("00_PHYLOSEQ_OBJECTS/00_seqtab.nochim_18S.rds")

seqtab.nochim_16S <- seqtab.nochim_16S[-c(112,113,114),]

seqtab.nochim_all <- cbind(seqtab.nochim_16S, seqtab.nochim_18S)
dim(seqtab.nochim_all)[2] #74370

taxa_18S %<>% as.data.frame() %>% dplyr::select(Kingdom, Division, Class, Order, Family, Genus) %>% as.matrix
colnames(taxa_18S)[2] <- "Phylum"
taxa_all <- rbind(taxa_16S, taxa_18S)


sample <- sample_data(sample)
ps <- phyloseq(otu_table(seqtab.nochim_all, taxa_are_rows = FALSE),
               tax_table(taxa_all),sample )

dna <- Biostrings::DNAStringSet(taxa_names(ps)) 
names(dna) <- taxa_names(ps) 
ps <- merge_phyloseq(dna, ps) 
taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps))) 


saveRDS(ps, "00_PHYLOSEQ_OBJECTS/01_ps.rds")
```

### Remove plastid reads <a name="novaseq_remove_plastid"></a>
```r
tax <- ps_sup@tax_table %>% as.data.frame()

mito <- tax %>% as.data.frame() %>% filter(Family=="Mitochondria") #1219
otu_mito <- rownames(mito[1])

chloro <- tax %>% as.data.frame() %>% filter(Order=="Chloroplast") #1225
otu_chloroplast <- rownames(chloro[1])

plas <- tax %>% as.data.frame() %>% filter(Kingdom=="Eukaryota:plas") #809
otu_plas <- rownames(plas[1])

all_ASV = taxa_names(ps_sup)
all_ASV <- all_ASV[!(all_ASV %in% otu_chloroplast)]
all_ASV <- all_ASV[!(all_ASV %in% otu_mito)]
all_ASV <- all_ASV[!(all_ASV %in% otu_plas)]

ps_sup_clean <- ps_sup
ps_sup_clean=prune_taxa(all_ASV, ps_sup)
ps_sup_clean

ps_sup_clean = merge_phyloseq(ps_sup_clean, sample_data(sample))

saveRDS(ps_sup_clean, "00_PHYLOSEQ_OBJECTS/03_ps_sup_clean.rds")
```


### Decontamination with Microdecon <a name="novaseq_microdecon"></a>

```r
## change first column 
tax <- ps_sup_clean@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_sup_clean@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"


otu %<>% as.data.frame() %>%  relocate("16S_T_neg_MARS_BACT_S275_R1", .before="16S_CECILE_1_AVRIL_BACT_S281_R1") 
otu %<>% as.data.frame() %>%  relocate("16S_T_neg_pool_2_BACT_S274_R1", .before="16S_CECILE_1_AVRIL_BACT_S281_R1") 


otu_tax <- merge(otu, tax)
otu_tax %<>% unite("taxonomy", Kingdom:Genus, sep = ";")

# Transform column in numeric
otu_tax[, 2:124] <- sapply(otu_tax[, 2:124], as.numeric)

clean_data <-  decon(data = otu_tax, numb.blanks = 3, numb.ind = c(120), taxa = TRUE, thresh = 1)

#ps
otu <- clean_data$decon.table %>% as.data.frame()
otu %<>% select(-c(Mean.blank, taxonomy))

otu2 <- otu[,-1]
rownames(otu2) <- otu[,1]

rm(otu)
otu <- otu2
rm(otu2)

otu <- as.matrix(otu)


OTU = otu_table(otu, taxa_are_rows = TRUE)
 
ps_decon <- phyloseq(OTU, ps_sup_clean@sam_data, ps_sup_clean@tax_table, ps_sup_clean@refseq)
ps_decon

saveRDS(ps_decon, "00_PHYLOSEQ_OBJECTS/04_ps_decon.rds")
```

### Remove Homopolymer <a name="novaseq_homopolymer"></a>
```r
ref <- ps_decon@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps_decon@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_decon@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref, "01_Tables/05_otu_tax_ref_after_microdecon.csv", sep=";", quote=FALSE)
# in this table, remove by hand all sequences which starts with more than 5 C. 


otu_tax_ref <- read_tsv("01_Tables/07_otu_tax_ref_without_homopolymer.tsv")
otu_tax_ref %<>% tibble::column_to_rownames("ASV")

otu <- otu_tax_ref %>% dplyr::select(`16S_CECILE_1_AVRIL_BACT_S281_R1`:`16S_R5_male_BACT_S270_R1`) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Genus) %>% as.matrix()

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)

ps_decon_wt_cc <- phyloseq(OTU, TAX, ps_decon@refseq, ps_decon@sam_data)

saveRDS(ps_decon_wt_cc, "00_PHYLOSEQ_OBJECTS/05_ps_wt_homopoly.rds")
```

### Remove chimeras vsearch <a name="novaseq_vsearch"></a>
```bash
#!/usr/bin/env bash
#SBATCH --job-name=chimera
#SBATCH --partition fast
#SBATCH --mem 5G
#SBATCH --cpus-per-task 3 
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type END

module load vsearch


input=03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/08_CHIMERA_VSEARCH/bact_euk.fasta
output_chimera=03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/08_CHIMERA_VSEARCH/chimera.fasta
output_no_chimera=03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/08_CHIMERA_VSEARCH/no_chimera.fasta
db=databases/SILVA_138.1_SSURef_NR99_tax_silva.fasta

vsearch --threads 3 --uchime_ref ${input} --chimera ${output_chimera} --nonchimeras ${output_no_chimera} --db ${db}
```

```r
asv_no_chimera <- read.table("01_Tables/08_asv_no_chimera.tsv", sep="\t", header=FALSE, row.names = 1)
asv_no_chimera <- rownames(asv_no_chimera[2])

all_ASV = taxa_names(ps_decon_wt_cc)

all_ASV <- all_ASV[(all_ASV %in% asv_no_chimera)]

ps_decon_wt_cc_no_chimera <- ps_decon_wt_cc
ps_decon_wt_cc_no_chimera=prune_taxa(all_ASV, ps_decon_wt_cc)
ps_decon_wt_cc_no_chimera

saveRDS(ps_decon_wt_cc_no_chimera, "00_PHYLOSEQ_OBJECTS/06_ps_no_chimeras.rds")
```

### Identity Ascophyllum ASVs <a name="novaseq_asco"></a>

```r
ref <- ps_decon_wt_cc_no_chimera@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps_decon_wt_cc_no_chimera@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_decon_wt_cc_no_chimera@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref, "01_Tables/09_otu_tax_ref_after_vsearch.csv", sep=";", quote=FALSE)

#all Silvetia + Fucus has been merged + ASV_73217, ASV_70703,ASV_71087, ASV_71721, ASV_70539, ASV_70858,  ASV_71415, ASV_70581, ASV_71374, ASV_73644 as there 99% identity with fucus or silvetia. (clustering did with vsearch)

asco <- read_tsv("01_Tables/10_otu_tax_ref_after_vsearch_1asv_asco.tsv")
  
otu <- asco %>% dplyr::select(ASV,`16S_CECILE_1_AVRIL_BACT_S281_R1`:`16S_R5_male_BACT_S270_R1`)
otu %<>% tibble::column_to_rownames("ASV") %>% as.matrix()
tax <- asco %>% dplyr::select(ASV,Kingdom:Genus) 
tax %<>% tibble::column_to_rownames("ASV") %>% as.matrix()

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)

ps_decon_asco <- phyloseq(OTU, TAX, ps_decon_wt_cc_no_chimera@refseq, ps_decon_wt_cc_no_chimera@sam_data)

saveRDS(ps_decon_asco, "00_PHYLOSEQ_OBJECTS/07_ps_decon_asco.rds")
```

### Final novaseq ps object for sequencing comparison <a name="novaseq_ps_final"></a>

```r
ps_decon_asco <- readRDS("../01_Microdecon_all/00_PHYLOSEQ_OBJECTS/07_ps_decon_asco.rds")
sample <- read.table("01_Tables/01_samples_sites.csv", sep=";", header=TRUE, dec=".", row.names = 1)

ps_sites_asco = phyloseq(ps_decon_asco@otu_table, ps_decon_asco@tax_table, ps_decon_asco@refseq, sample_data(sample))

ps_sites_asco = filter_taxa(ps_sites_asco, function(x) sum(x) > 0, TRUE)
ps_sites_asco

saveRDS(ps_sites_asco,"00_PHYLOSEQ_OBJECTS/04_ps_sites_asco.rds")
```

## 2. Fungi <a name="novaseq_fungi"></a>
### Cutadapt <a name="novaseq_fungi_cutadapt"></a>
```bash
#!/usr/bin/env bash
#SBATCH --job-name=cutadapt
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load cutadapt

PRIMER_F=AACTTTYRRCAAYGGATCWCT
PRIMER_R=AGCCTCCCGCTTATTGATATGCTTAART
MIN_LENGTH=220
output=03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/01_CUTADAPT

for sample in $(cat sample.txt); \
do \
   cutadapt  -m ${MIN_LENGTH} --discard-untrimmed --report=minimal -g ${PRIMER_F} -G ${PRIMER_R} ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz -o ${output}/${s
ample}_R1_cutadapt.fastq.gz -p ${output}/${sample}_R2_cutadapt.fastq.gz; done
```

### DADA2 <a name="novaseq_fungi_dada2"></a>

```r
library(dada2); packageVersion("dada2")


path="03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/01_CUTADAPT"
files <- list.files(path)

fnFs_ITS <- sort(list.files(path, pattern="_R1_cutadapt.fastq.gz", full.names = TRUE))
fnRs_ITS <- sort(list.files(path, pattern="_R2_cutadapt.fastq.gz", full.names = TRUE))

sample.names_ITS <- sapply(strsplit(basename(fnFs_ITS), "_R"), `[`, 1)

filtFs_ITS <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_ITS, "_R1_filt_dada2.fastq.gz"))
filtRs_ITS <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names_ITS, "_R2_filt_dada2.fastq.gz"))
names(filtFs_ITS) <- sample.names_ITS
names(filtRs_ITS) <- sample.names_ITS

out_ITS <- filterAndTrim(fnFs_ITS, filtFs_ITS, fnRs_ITS, filtRs_ITS, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
out_ITS

save.image("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/01_out_filt_split_ITS.rds")

#load("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/01_out_filt_split_ITS.rds")

errF_ITS <- learnErrors(filtFs_ITS, multithread=TRUE)
errR_ITS <- learnErrors(filtRs_ITS, multithread=TRUE)
dadaFs_ITS <- dada(filtFs_ITS, err=errF_ITS, multithread=TRUE)
dadaRs_ITS <- dada(filtRs_ITS, err=errR_ITS, multithread=TRUE)

mergers_ITS <- mergePairs(dadaFs_ITS, filtFs_ITS, dadaRs_ITS, filtRs_ITS, verbose=TRUE)
seqtab_ITS <- makeSequenceTable(mergers_ITS)
seqtab.nochim_ITS <- removeBimeraDenovo(seqtab_ITS, method="consensus", multithread=TRUE, verbose=TRUE)



getN_ITS <- function(x) sum(getUniques(x))
track_ITS <- cbind(out_ITS, sapply(dadaFs_ITS, getN_ITS), sapply(dadaRs_ITS, getN_ITS), sapply(mergers_ITS, getN_ITS), rowSums(seqtab.nochim_ITS))
colnames(track_ITS) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_ITS) <- sample.names_ITS
track_ITS

taxa_ITS <- assignTaxonomy(seqtab.nochim_ITS, "03_NOVASEQ_METAB/02_DATABASES/sh_general_UNITE_release_dynamic_25.07.2023.fasta",t
ryRC=TRUE, multithread=TRUE)

taxa.print_ITS <- taxa_ITS
rownames(taxa.print_ITS) <- NULL

# SAVE 
save.image("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/04_dada_pipeline_ITS.RData")
saveRDS(taxa_ITS, "03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/05_taxa_ITS.rds")
saveRDS(seqtab.nochim_ITS, "03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/06_seqtab.nochim_ITS.rds")
```
