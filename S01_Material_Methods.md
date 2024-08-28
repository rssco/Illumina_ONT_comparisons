## A. MiSeq read processing






## B. NovaSeq read processing
### 1. Bacteria

#### Cutadapt
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

#### Split 16S/18S
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

#### DADA2 16S
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
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/02_err_16S.RData")
dadaFs_16S <- dada(filtFs_16S, err=errF_16S, multithread=TRUE)
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/03_dadaFs_16S.RData")
dadaRs_16S <- dada(filtRs_16S, err=errR_16S, multithread=TRUE)
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/04_dadaRs_16S.RData")
mergers_16S <- mergePairs(dadaFs_16S, filtFs_16S, dadaRs_16S, filtRs_16S, verbose=TRUE)
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/05_mergers_16S.RData")
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


load("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/07_tracks_16S.RData")
taxa_16S <- assignTaxonomy(seqtab.nochim_16S, "03_NOVASEQ_METAB/02_DATABASES/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa.print_16S <- taxa_16S # Removing sequence rownames for display only
rownames(taxa.print_16S) <- NULL

# SAVE 
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/08_assigntaxo_16S.RData")
saveRDS(taxa_16S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/09_taxa_16S.rds")
saveRDS(seqtab.nochim_16S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/10_seqtab.n
ochim_16S.rds")
```

#### DADA2 18S

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


saveRDS(out_18S, "03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/11_out_filt_split_18S
.rds")
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/12_out_filt_18S.RData")

errF_18S <- learnErrors(filtFs_18S, multithread=TRUE)
errR_18S <- learnErrors(filtRs_18S, multithread=TRUE)
dadaFs_18S <- dada(filtFs_18S, err=errF_18S, multithread=TRUE)
dadaRs_18S<- dada(filtRs_18S, err=errR_18S, multithread=TRUE)
save.image("03_NOVASEQ_METAB/03_BACTERIA_ANALYSIS/02_ALL/02_SPLIT_16S_18S_PIPELINE/05_RSTUDIO_OUTPUTS/13_err_dada_18S.RData")

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


### 2. Fungi
#### Cutadapt
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

#### DADA2

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
save.image("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/02_err_ITS.RData")
dadaFs_ITS <- dada(filtFs_ITS, err=errF_ITS, multithread=TRUE)
dadaRs_ITS <- dada(filtRs_ITS, err=errR_ITS, multithread=TRUE)
save.image("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/03_dada_ITS.RData")

mergers_ITS <- mergePairs(dadaFs_ITS, filtFs_ITS, dadaRs_ITS, filtRs_ITS, verbose=TRUE)
save.image("03_NOVASEQ_METAB/04_FUNGI_ANALYSIS/02_RSTUDIO_OUTPUTS/04_mergers_ITS.RData")
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

## C. ONT read processing
### 1. Bacteria
#### Cutadapt

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

PRIMER_F=AGAGTTTGATCMTGGCTCAG
#PRIMER_R=AAGTCGTAACAAGGTAACC #REVERSE COMPLEMENT
#MAE=0.04
#E=0.2
#MIN_LENGTH=800

#for i in 01_BACTERIA/*gz;
#do cutadapt -O ${#PRIMER_F} -e ${E} --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --revcomp --report=minimal -g ${PRIMER_F} ${i} | \
#   cutadapt - -O ${#PRIMER_R} -e ${E}  --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --report=minimal -a ${PRIMER_R} -o ${i%.fastq.gz}_cutadap
t.fastq.gz; done
```


### 2. Fungi
#### Cutadapt

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
##SBATCH --array 1-1

module load cutadapt

PRIMER_F=TGTACACACCGCCCGTCG
PRIMER_R=GTCTTGAAACACGGACCA #REVERSE COMPLEMENT
E=0.2

for i in /09_NGSPECIESID/01_QUALITY_FILTERING/*_chopper.fastq;
do cutadapt -O ${#PRIMER_F} -e ${E} --discard-untrimmed --revcomp --report=minimal -g ${PRIMER_F} ${i} | \
   cutadapt - -O ${#PRIMER_R} -e ${E} --discard-untrimmed --revcomp --report=minimal -a ${PRIMER_R} -o ${i%.fastq}_cutadapt.fastq; done
```











## A. Map sample collection

### 0. Libraries
```r include=FALSE
library(ggmap)
library(ggplot2)
```

### 1. Tables
```r
lieux <- data.frame(Lieux=c('Lieu 1', 'Lieu 2', 'Lieu 3'),
                             Latitude=c(48.84941, 48.85253, 48.85346),
                             Longitude=c(-3.078291, -3.075012, -3.077533))

villes <- data.frame(Lieux=c('Brest', 'Pleubian'),
                             Latitude=c(48.40, 48.81),
                             Longitude=c(-4.47798, -3.15305))
```

### 2. Maps
```{r}
register_stadiamaps("YOUR_API_KEY", write = TRUE)

Pleubian <- get_stadiamap(bbox = c(left = -3.095, bottom = 48.846, right = -3.060, top = 48.8641), zoom = 16, maptype = "stamen_watercolor") %>% ggmap() + 
  geom_point(data=lieux, aes(x = Longitude, y=Latitude),size=1.5) +
  geom_text(aes(x=Longitude, y=Latitude, label=Lieux), data=lieux, hjust=1.2, size=2.5) +
  xlab(expression(paste("Longitude"))) +
  ylab(expression(paste("Latitude"))) +
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )

Brittany <- get_stadiamap(bbox = c(left = -5, bottom = 47.5, right = -3, top = 49), zoom = 9, maptype = "stamen_watercolor") %>% ggmap() + 
  geom_point(data=villes, aes(x = Longitude, y=Latitude),size=1.5) +
  geom_text(aes(x=Longitude, y=Latitude, label=Lieux), data=villes, hjust=1.2, size=2.5) +
  xlab(expression(paste("Longitude"))) +
  ylab(expression(paste("Latitude"))) +
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )
```

![Figure 1 | Sampling locations.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/01_Methods_sampling.png)<!-- -->
