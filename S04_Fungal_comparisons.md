
# 1. Libraries
```r
library(phyloseq)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(magrittr)
library(microbiome)
library(RColorBrewer)
library(microDecon)
library(paletteer)
library(gridExtra)
library(vegan)
library(rstatix)
library(ggpubr)
library(ape)
library(ggrepel)
library(reshape2)
library(ggh4x)
library(colorspace)
```

# 2. Tables
```r
ps_novaseq <- readRDS("~/Documents/ownCloud/11_Rstudio/08_Novaseq/05_Fungi_analysis/02_Sites/00_Phyloseq_objects/01_ps_sites_fungi.rds") #after microdecon and removing CCC
ps_miseq <- readRDS("~/Documents/ownCloud/11_Rstudio/05_Fungi_miseq_metaB/05_Envi_Phyloseq_object/06_ps_decon_wtout_mock2.rds") #after microdecon


ps_ms_ns_ont <- readRDS("00_Phyloseq_objects/01_ps_ms_ns_ont.rds") # ps with miseq,  novaseq and ONT

ps_ms_ns_ont_trans <- readRDS("00_Phyloseq_objects/02_ps_ms_ns_ont_trans.rds") # ps with miseq,  novaseq and ONT
```


# 3. Generate one ps per all sequencing methods
## Novaseq 
```r
otu_ns <- ps_novaseq@otu_table %>% as.data.frame() 
tax_ns <- ps_novaseq@tax_table %>% as.data.frame()
ref_ns <- ps_novaseq@refseq %>% as.data.frame()

write.table(otu_ns, "01_Tables/otu_ns.csv", sep=";", quote=FALSE)
write.table(tax_ns, "01_Tables/tax_ns.csv", sep=";", quote=FALSE)
write.table(ref_ns, "01_Tables/ref_ns.csv", sep=";", quote=FALSE)
```

## Miseq

```r
otu_ms <- ps_miseq@otu_table %>% as.data.frame() #434 ASVs, 24 samples
tax_ms <- ps_miseq@tax_table %>% as.data.frame()
ref_ms <- ps_miseq@refseq %>% as.data.frame()

write.table(otu_ms, "01_Tables/otu_ms.csv", sep=";", quote=FALSE)
write.table(tax_ms, "01_Tables/tax_ms.csv", sep=";", quote=FALSE)
write.table(ref_ms, "01_Tables/ref_ms.csv", sep=";", quote=FALSE)

```

## Combine Miseq, NovaSeq, ONT together by hand, then phyloseq object
```r
sample_ms_ns_ont <- read.table("01_Tables/02_sample_ms_ns_ont.csv", header=TRUE, sep=";", row.names = 1)
otu_ms_ns_ont <- read.table("01_Tables/01_otu_ms_ns_ont.csv", header=TRUE, sep=";", row.names = 1)
tax_ms_ns_ont <- read.table("01_Tables/03_tax_ms_ns_ont.csv", header=TRUE, sep=";", row.names = 1)

otu_ms_ns_ont %<>% as.matrix()
tax_ms_ns_ont %<>% as.matrix()

OTU = otu_table(otu_ms_ns_ont, taxa_are_rows = TRUE)
TAX = tax_table(tax_ms_ns_ont)
SAMPLE = sample_data(sample_ms_ns_ont)

ps_ms_ns_ont <- phyloseq(OTU, TAX, SAMPLE)
ps_ms_ns_ont = filter_taxa(ps_ms_ns_ont, function(x) sum(x) > 0, TRUE)

saveRDS(ps_ms_ns_ont, "00_Phyloseq_objects/01_ps_ms_ns_ont.rds")
```

# 4. Normalization

```r
total = median(sample_sums(ps_ms_ns_ont)) # Median= 144989
standf = function(x, t=total) round(t * (x / sum(x))) # Standardize abundances to the median sequencing depth
ps_filt = transform_sample_counts(ps_ms_ns_ont, standf)

ps_ms_ns_ont_trans = filter_taxa(ps_filt, function(x) sum(x > total* 0.0003) > 0, TRUE)


saveRDS(ps_ms_ns_ont_trans, "00_Phyloseq_objects/02_ps_ms_ns_ont_trans.rds")
```

# 5. Barplot at genus level 

```r
ps_agglomerate <- tax_glom(ps_ms_ns_ont_trans, taxrank = 'Species', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Species <- as.character(taxa$Species)



species <- c("Mycophycias_ascophylli", "Moheitospora_sp", "Metschnikowia_bicuspidata", "Malassezia_restricta")

taxa$Species[!taxa$Species %in% species] <- "Other fungi"

taxa$Sequencing_type <- factor(taxa$Sequencing_type, levels = c("MISEQ", "NOVASEQ", "ONT"))
taxa$Month <- factor(taxa$Month, levels = c("November", "March"))
taxa %<>% filter(!Sample=="L1I5_NOV21_NS")
taxa %<>% filter(!Sample=="L1I3_NOV_MS")
taxa %<>% filter(!Sample=="L1I3_NOV_ONT")

taxa$Species <- factor(taxa$Species, levels = c('Other fungi',
                                            "Malassezia_restricta",
                                            "Metschnikowia_bicuspidata",
                                            "Moheitospora_sp",
                                            "Mycophycias_ascophylli"))
                                            



HowManyPhyla<-5
getPalette = colorRampPalette(c("#2A5676" , "#CB8215","#00441B","#645A9F","#767676")) 
PhylaPalette = rev(getPalette(HowManyPhyla))

ggplot(data=taxa, aes(x=Alga, y=Abundance, fill=Species))  +
  facet_nested(Sequencing_type ~ Month+Site, space = "free", scales = "free") +
  geom_bar(stat="identity",  position = "fill") +
  scale_color_manual(values = PhylaPalette ) + scale_fill_manual(values = PhylaPalette) +
  theme(axis.line = element_line(colour = "black"),
        strip.text=element_text(size=15),
        legend.title = element_text(size=15),
        legend.key.size = unit(0.9, 'cm'),
        legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    legend.text = element_text(size=15),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_blank(),
    axis.title=element_text(size=20)) +
 guides(fill=guide_legend(ncol = 2)) +
   ylab("Relative abundance") +
  pdf("02_Figures/01_barplots_fungi_species.pdf", width = 10, height = 14 )
```

![Figure 6 | Fungal comparisons.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/05_combined_mock_consensus_barplot.png)<!-- -->
