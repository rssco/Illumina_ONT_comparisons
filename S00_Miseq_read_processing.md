# MiSeq read processing

## Table of content  
[1. Bacteria](#miseq_bacteria)  
- [DADA2](#miseq_bacteria_dada2)   
- [Phyloseq](#miseq_phyloseq)  
- [Decontamination with microdecon](#miseq_decontamination)  
- [Remove mock](#miseq_remove_mock)  
- [Remove plastid reads](#miseq_remove_plastid_reads)  
[2. Fungi](#miseq_fungi)  
- [Cutadapt](#miseq_fungi_cutadapt)  
- [DADA2](#miseq_fungi_dada2)  
- [Phyloseq](#miseq_fungi_ps)
- [Decontamination with Micrdecon](#miseq_fungi_microdecon)
- [Remove mock](#miseq_fungi_remove_mock)


## 1. Bacteria <a name="miseq_bacteria"></a>
### DADA2 <a name="miseq_bacteria_dada2"></a>

```r
# Library
library(tidyverse)
library(dada2)
library(DECIPHER)

# Data
path="300_200/01_BACTERIA/"
files <- list.files(path)

path_filt="300_200/01_BACTERIA/01_CLEANED_DATA/"
filtRs <- list.files(path_filt, pattern = "_R_filt.fastq.gz", full.names = TRUE)
filtFs <- list.files(path_filt, pattern = "_F_filt.fastq.gz", full.names = TRUE)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "01_CLEANED_DATA/", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Read filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=c(17,19),
              compress=TRUE, multithread=TRUE) 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# Assign taxonomy 
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

### Phyloseq <a name="miseq_phyloseq"></a>

Ps construction
```r
#Library
library(phyloseq)
library(microdecon)

sample <- read.table("01_Tables/sample.csv", sep=";", header=TRUE, row.names = 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(sample),
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps)) 
names(dna) <- taxa_names(ps) 
ps <- merge_phyloseq(ps, dna) 
taxa_names(ps) <- paste0("OTU", seq(ntaxa(ps))) 

saveRDS(ps, "00_Phyloseq_objects/01_ps.rds")
```

### Decontamination with microdecon <a name="miseq_decontamination"></a>
```r
tax <- ps@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "otu_id"

otu <- ps@otu_table %>% as.data.frame() %>% t()
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "otu_id"


otu %<>% as.data.frame() %>%  relocate("T-neg-NOV-BACT", .before="L1I1-MARS-BACT") 
otu %<>% as.data.frame() %>%  relocate("T-neg-MARS-BACT", .before="L1I1-MARS-BACT")


otu_tax <- merge(otu, tax)
otu_tax %<>% unite("taxonomy", Kingdom:Species, sep = ";")

# Transform column in numeric
otu_tax[, 2:29] <- sapply(otu_tax[, 2:29], as.numeric)

#run microdecon
clean_data <-  decon(data = otu_tax, numb.blanks = 2, numb.ind = c(24,2), taxa = TRUE, thresh = 1)
```

ps decontaminated
```r
tax <- ps@tax_table %>% as.data.frame()
otu <- clean_data$decon.table %>% as.data.frame()
sample <- ps@sam_data %>% as.data.frame()
otu %<>% select(-c(Mean.blank, taxonomy))

otu2 <- otu[,-1]
rownames(otu2) <- otu[,1]

rm(otu)
otu <- otu2


otu <- as.matrix(otu)
tax <- as.matrix(tax)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
samples = sample_data(sample)
 
ps_decon <- phyloseq(OTU, TAX, samples)
```

### Remove mock <a name="miseq_remove_mock"></a>
```r
ps_decon_wt_mock <- prune_samples(sample_names(ps_decon) != "MOCK-BACT-mars", ps_decon)
ps_decon_wt_mock <- prune_samples(sample_names(ps_decon) != "MOCK-BACT-nov", ps_decon)
```

### Remove plastid reads <a name="miseq_remove_plastid_reads"></a>
```r
tax <- ps_decon_wt_mock@tax_table %>% as.data.frame()

noaff <- tax %>% as.data.frame() %>% filter(is.na(Kingdom)) #175
otu_noaff <- rownames(noaff[1]) 

mito <- tax %>% as.data.frame() %>% filter(Family=="Mitochondria") #162
otu_mito <- rownames(mito[1])

chloro <- tax %>% as.data.frame() %>% filter(Order=="Chloroplast") #82
otu_chloroplast <- rownames(chloro[1])

all_ASV = taxa_names(ps_decon_wt_mock)
all_ASV <- all_ASV[!(all_ASV %in% otu_chloroplast)]
all_ASV <- all_ASV[!(all_ASV %in% otu_noaff)]
all_ASV <- all_ASV[!(all_ASV %in% otu_mito)]

ps_decon_wt_mock_clean <- ps_decon_wt_mock
ps_decon_wt_mock_clean=prune_taxa(all_ASV, ps_decon_wt_mock)
ps_decon_wt_mock_clean
```

## 2. Fungi <a name="miseq_fungi"></a>
### Cutadapt <a name="miseq_fungi_cutadapt"></a>
```r
#Library
library(dada2)
library(ShortRead)
library(Biostrings)
library(parallel)

# Data
path <- "02_Sequences/"
path_filt <- "02_Sequences/filtN/"
path.cut <- "02_Sequences/cutadapt/"
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Identify primers
FWD <- "AACTTTYRRCAAYGGATCWCT" 
REV <- "AGCCTCCCGCTTATTGATATGCTTAART"  

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

source("https://raw.githubusercontent.com/benjjneb/dada2/master/R/filter.R")
source("https://raw.githubusercontent.com/benjjneb/dada2/master/R/RcppExports.R")

path_filt="02_Sequences/filtN/"

fnFs.filtN <- list.files(path_filt, pattern = "_R1_001.fastq.gz", full.names=TRUE) 
fnRs.filtN <- list.files(path_filt, pattern = "_R2_001.fastq.gz", full.names=TRUE) 

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) 


# Cutadapt
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#path to cutadapt
cutadapt <- "../Library/Python/3.9/bin/cutadapt"

# Remove primers with cutadapt

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Verify no more primers}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
```

### DADA2 <a name="miseq_fungi_dada2"></a>

```r
# data
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
sample.names

filtFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

path_cutadapt_filt="02_Sequences/cutadapt/filtered/"

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Read filtering
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(260,230), maxN=0, maxEE=c(2,5), truncQ=2, minLen = 50, rm.phix=TRUE,compress=TRUE, multithread=TRUE) 

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
    getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
    "nonchim")
rownames(track) <- sample.names
track

# Assign taxonomy
unite.ref <- "03_Databases/sh_general_release_dynamic_25.07.2023.fasta "  
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# SAVE
saveRDS(seqtab.nochim,"05_Envi_Phyloseq_object/03_seqtab.nochim.RData")
```

### Phyloseq object construction <a name="miseq_fungi_ps"></a>

```r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(sample),
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps)) 
names(dna) <- taxa_names(ps) 
ps <- merge_phyloseq(ps, dna) 
taxa_names(ps) <- paste0("OTU", seq(ntaxa(ps))) 

saveRDS(ps, "05_Envi_Phyloseq_object/04_ps.rds" )
```


### Decontamination with Microdecon <a name="miseq_fungi_microdecon"></a>

```r
## change first column 
tax <- ps@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "otu_id"

otu <- ps@otu_table %>% as.data.frame() %>% t()
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "otu_id"


otu %<>% as.data.frame() %>%  relocate("T-neg-NOV-FUNGI", .before="L1I1-MARS-FUNGI") 
otu %<>% as.data.frame() %>%  relocate("T-neg-MARS-FUNGI", .before="L1I1-MARS-FUNGI")


otu_tax <- merge(otu, tax)
otu_tax %<>% unite("taxonomy", Kingdom:Species, sep = ";")

# Transform column in numeric
otu_tax[, 2:29] <- sapply(otu_tax[, 2:29], as.numeric)


clean_data <-  decon(data = otu_tax, numb.blanks = 2, numb.ind = c(24,2), taxa = TRUE, thresh = 1)
```

New ps object, decontaminated + saving
```r
#recup tables from ps 
tax <- ps@tax_table %>% as.data.frame()
otu <- clean_data$decon.table %>% as.data.frame()
sample <- ps@sam_data %>% as.data.frame()
otu %<>% select(-c(Mean.blank, taxonomy))

#First column as rownames
otu2 <- otu[,-1]
rownames(otu2) <- otu[,1]

rm(otu)
otu <- otu2
rm(otu2)

#transform tables
otu <- as.matrix(otu)
tax <- as.matrix(tax)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
samples = sample_data(sample)
 
#Make new ps object  
ps_decon <- phyloseq(OTU, TAX, samples, ps@refseq)

saveRDS(ps_decon, "05_Envi_Phyloseq_object/05_ps_decon.rds")
```

### Remove mock <a name="miseq_fungi_remove_mock"></a>

```r
ps_decon_wt_mock <- prune_samples(sample_names(ps_decon) != "MOCK-FUNGI-mars", ps_decon)
ps_decon_wt_mock <- prune_samples(sample_names(ps_decon_wt_mock) != "MOCK-FUNGI-nov", ps_decon_wt_mock)
ps_miseq <- ps_decon_wt_mock

saveRDS(ps_decon_wt_mock, "05_Envi_Phyloseq_object/06_ps_decon_wtout_mock.rds")
```