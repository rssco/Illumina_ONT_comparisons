# Performed assignation with SAMBA on Slurm cluster 

Same parameters with and without adaptative sampling

https://gitlab.ifremer.fr/bioinfo/workflows/samba
 
Changes in conf/nanopore.config
- Minimum length of raw nanopore reads to keep: nanopore_read_minlength = "0"  
- Maximale length of raw nanopore reads to keep: nanopore_read_maxlength = "100000"  
- minimap2_db = "../finalresult/02_ONT/06_SAMBA/samba/database/silva_v138.1_16S_NR99_SEQ_k15.mmi"  
- ref_tax = "../finalresult/02_ONT/06_SAMBA/samba/database/silva_v138.1_16S_NR99_TAX.txt"  
- Taxonomic level to analyse data. Can be: Kingdom;Phylum;Class;Order;Genus;Species: tax_rank = "Species"  
- Define the Kingdom studied e.g for 16S nanopore reads: "Bacteria" or for 18S nanopore reads: "Eukaryota": kingdom = "Bacteria"  

```bash
#!/usr/bin/env bash
#SBATCH --job-name=samba
#SBATCH -p fast
#SBATCH --mem=150G
#SBATCH --cpus-per-task 28
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL


module load nextflow slurm-drmaa graphviz


nextflow run ../finalresult/02_ONT/00_BASECALLING_GUPPY_ADAPTATIVE_SAMPLING_v6.1.5/03_SAMBA/samba/main.nf -profile singularity,nanopore -c ../finalresult/02_ONT/00_BASECALLING_GUPPY_ADAPTATIVE_SAMPLING_v6.1.5/03_SAMBA/samba/abims.config
```


# On Rstudio
## 1. Libraries
```r
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(ggh4x)
library(ggforce)
```

## 2. Without adaptative sampling
```r
#load ps + read length table
ps_without_adapt <- readRDS("02_Phyloseq_objects/03_ps_witout_adapt_SAMBA.rds")
read_length <- read.table("01_Tables/01_read_id_lenght_without_adapt.tsv", sep="\t", header=TRUE)

# Add read length
otu <- ps_with_adapt@otu_table %>% as.data.frame()
otu <- cbind(rownames(otu), otu)
colnames(otu)[1] <- "read_id"
tax <- ps_with_adapt@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
colnames(tax)[1] <- "read_id"

merged <- merge(otu, read_length, by="read_id")
merged <- merge(merged, tax, by="read_id")

# Plot 
class <- c("Cyanobacteriia", "Gammaproteobacteria","Alphaproteobacteria")

merged$Class[!merged$Class %in% class] <- "Other bacteria"

palette_class <- c(
                "Other bacteria"="grey",
                  "Cyanobacteriia"="#00441B",
                  "Gammaproteobacteria"="#B03F00",
                "Alphaproteobacteria"= "#08306B")

merged$Class <- factor(merged$Class, levels = c('Other bacteria',
                                          "Cyanobacteriia",
                                          "Gammaproteobacteria",
                                          "Alphaproteobacteria"))

without <- ggplot(merged, aes(x=pb, fill=Class)) +  geom_histogram(position = "stack") +
  scale_x_continuous(limits=c(0, 2000))+
      scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    aspect.ratio = 1) +
     ylab("Number of reads") +
  xlab("Read length (pb)")+
  theme_bw() +
  ggtitle("A. Without adaptative sampling")

# Plot Zoom
zoom_without <- ggplot(merged, aes(x = pb, fill=Class)) + geom_histogram(position = "stack") +
    scale_x_continuous(limits=c(0, 2000))+
      scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    aspect.ratio = 1) +
     ylab("Number of reads") +
  xlab("Read length (pb)")+
  theme_bw() +
  ggtitle("C. Without adaptative sampling")+
facet_zoom(ylim = c(0,100000), split = TRUE, zoom.size = 4)

```

## 2. With adaptative sampling
```r
#load ps + read length table
ps_with_adapt <- readRDS("02_Phyloseq_objects/04_ps_with_adapt_SAMBA.rds")
read_length <- read.table("01_Tables/03_read_id_length_with_adapt.tsv", sep="\t", header=TRUE)

#Add read length
otu <- ps_with_adapt@otu_table %>% as.data.frame()
otu <- cbind(rownames(otu), otu)
colnames(otu)[1] <- "read_id"
tax <- ps_with_adapt@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
colnames(tax)[1] <- "read_id"

merged <- merge(otu, read_length, by="read_id")
merged <- merge(merged, tax, by="read_id")

# Plot
class <- c("Cyanobacteriia", "Gammaproteobacteria","Alphaproteobacteria")

merged$Class[!merged$Class %in% class] <- "Others bacteria"

palette_class <- c(
                "Others bacteria"="grey",
                  "Cyanobacteriia"="#00441B",
                  "Gammaproteobacteria"="#B03F00",
                "Alphaproteobacteria"= "#08306B")

merged$Class <- factor(merged$Class, levels = c('Others bacteria',
                                          "Cyanobacteriia",
                                          "Gammaproteobacteria",
                                          "Alphaproteobacteria"))

with <- ggplot(merged, aes(x=pb, fill=Class)) +  geom_histogram(position = "stack") +
  scale_x_continuous(limits=c(0, 2000))+
      scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    aspect.ratio = 1) +
     ylab("Number of reads") +
  xlab("Read length (pb)")+
  theme_bw() +
  ggtitle("B. With adaptative sampling")

#Plot Zoom
zoom_with <- ggplot(merged, aes(x = pb, fill=Class)) + geom_histogram(position = "stack") +
    scale_x_continuous(limits=c(0, 2000))+
      scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    aspect.ratio = 1) +
     ylab("Number of reads") +
  xlab("Read length (pb)")+
  theme_bw() +
  ggtitle("D. With adaptative sampling")+
facet_zoom(ylim = c(0,100000), split = TRUE, zoom.size = 4)
```

## 3. Combined figures

```r
ab <- gridExtra::grid.arrange(without, with, ncol=2)
cd <- gridExtra::grid.arrange(zoom_wihtout, zoom_with, ncol=2)
z <- gridExtra::grid.arrange(ab, cd, nrow=2)
ggsave("04_Figures/03_combined_normal_zoom_with_without.pdf", z ,width = 20, height = 10)
```

![Figure 4 | Adaptative sampling.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/03_combined_normal_zoom_with_without.png)<!-- -->
