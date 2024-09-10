# Mock analysis with IDTAXA, SAMBA and KRAKEN2 
## 0. Library
```r
library(tidyverse)
library(magrittr)
library(phyloseq)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(paletteer)
library(ggh4x)
library(colorspace)
library(gghighlight)
```

## 1. Tables 

```r
kraken <- read.table("02_Tables/01_abundance_table_kraken2.csv", sep=";", header=TRUE, dec = ".", row.names = 1)
ps_samba <- readRDS("02_Tables/02_phyloseq_samba_all_assignation_Species.rds")
idtaxa <- read.table("02_Tables/03_idtaxa_abundance_table.csv", sep=";", header=TRUE, dec=".", row.names = 1)
```

## 2. Transform tables
### Transform kraken2 table
```r
taxonomy <- kraken %>% select(taxonomy) %>%  separate(taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";",convert=TRUE)

otu <- kraken %>% dplyr::select(-c(taxonomy,observation_name,observation_sum))

write.table(taxonomy, "02_Tables/01_tax_kraken2.csv", sep=";", quote=FALSE)
write.table(otu, "02_Tables/01_otu_kraken2.csv", sep=";", quote=FALSE)
```

### Transform samba ps

```r
otu_samba <- ps_samba@otu_table %>% as.data.frame()
tax_samba <- ps_samba@tax_table %>% as.data.frame()

write.table(otu_samba, "02_Tables/02_otu_samba.csv", sep=";", quote=FALSE)
write.table(tax_samba, "02_Tables/02_tax_samba.csv", sep=";", quote=FALSE)
```

### Transform idtaxa table

```r
taxonomy <- idtaxa %>% select(consensus_taxo) %>%  separate(consensus_taxo,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";",convert=TRUE)

otu <- idtaxa %>% dplyr::select(-c(consensus_taxo))

write.table(taxonomy, "02_Tables/03_tax_idtaxa.csv", sep=";", quote=FALSE)
write.table(otu, "02_Tables/03_otu_idtaxa.csv", sep=";", quote=FALSE)
```


## 3. Combined phyloseq object
Concatenate by hand otu/tax/sample table PR2+SAMBA+KRAKEN2

```r
otu <- read.table("02_Tables/04_otu_all.csv", header=TRUE, sep=";", row.names = 1)
sample <- read.table("02_Tables/04_sample_all.csv", header=TRUE, sep=";", row.names = 1)
tax <- read.table("02_Tables/04_tax_all.csv", header=TRUE, sep=";", row.names = 1)

tax <- as.matrix(tax)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
SAMPLE = sample_data(sample)


ps_idtaxa_samba_kraken2 <- phyloseq(OTU, TAX, SAMPLE)
ps_idtaxa_samba_kraken2
```
## SAVE 
```r
saveRDS(ps_idtaxa_samba_kraken2, "00_Phyloseq_objects/01_ps_idtaxa_samba_kraken2.rds")

# Remove NOV MOCK + L2I4_MARS
ps_idtaxa_samba_kraken2 <- subset_samples(ps_idtaxa_samba_kraken2, ID!="MOCK_NOV" & ID!="L2I4_MARS")
saveRDS(ps_idtaxa_samba_kraken2, "00_Phyloseq_objects/02_ps_idtaxa_samba_kraken2_ssNOV_L2I4.rds")
```

## 3. Plot mock comparisons 
### Plot Order
```r
ps_mock <- subset_samples(ps_idtaxa_samba_kraken2, Individual=="MOCK")

ps_mock = filter_taxa(ps_mock, function(x) sum(x) > 0, TRUE)
ps_mock

ps_agglomerate <- tax_glom(ps_mock, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)

taxa[is.na(taxa)] <- "Non assigned"

taxa %<>% mutate(Kingdom = ifelse(grepl("Unknown ", Kingdom), "unclassified_", Kingdom)) %>% 
   mutate(Phylum = ifelse(grepl("Unknown ", Phylum), sub("Unknown ", "unclassified_", Phylum), Phylum)) %>% 
  mutate(Class = ifelse(grepl("Unknown ", Class), sub("Unknown ", "unclassified_", Class), Class)) %>% 
  mutate(Order = ifelse(grepl("Unknown ", Order), sub("Unknown ", "unclassified_", Order), Order)) %>%  
  mutate(Family = ifelse(grepl("Unknown ", Family),sub("Unknown ", "unclassified_", Family), Family)) %>% 
  mutate(Genus = ifelse(grepl("Unknown ", Genus), sub("Unknown ", "unclassified_", Genus), Genus))

class <- c("Bacilli","Gammaproteobacteria","Alphaproteobacteria","Non assigned", "unclassified_Bacteria", "unclassified_Firmicutes", "unclassified_Proteobacteria")
taxa$Class[!taxa$Class %in% class] <- "Not present"

taxa$Sample <- factor(taxa$Sample, levels = c("MOCK_MARS_BACT_ONT_IDTAXA", "MOCK_MARS_BACT_ONT_SAMBA","MOCK_MARS_BACT_ONT_KRAKEN2"))
taxa$Class <- factor(taxa$Class, levels = c("Bacilli","Gammaproteobacteria","Alphaproteobacteria", "unclassified_Bacteria", "unclassified_Firmicutes", "unclassified_Proteobacteria","Not present", "Non assigned"))


order <- ggplot(data=taxa, aes(x=Sample, y=Abundance, fill=Class))  +
  ggtitle("A. Mock at order level")+
  geom_bar(stat="identity", alpha=0.9,  position = "fill") +
  scale_x_discrete(labels=c('IDTAXA', 'SAMBA', 'KRAKEN2'))+
  scale_color_manual(values = c("Alphaproteobacteria"= "#08306B", 
                              "Gammaproteobacteria"="#B03F00",
                              "Bacilli" = "#984EA3",
                              "unclassified_Bacteria" = "#CCCCCC",
                              "unclassified_Firmicutes" = "#CCCCCC",
                              "unclassified_Proteobacteria" = "#CCCCCC",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999")) +
  scale_fill_manual(values = c("Alphaproteobacteria"= "#08306B", 
                              "Gammaproteobacteria"="#B03F00",
                              "Bacilli" = "#984EA3",
                              "unclassified_Bacteria" = "#CCCCCC",
                              "unclassified_Firmicutes" = "#CCCCCC",
                              "unclassified_Proteobacteria" = "#CCCCCC",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
    legend.text = element_text(size=12),
    axis.title=element_text(size=12),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    legend.title = element_text(size=13), 
    aspect.ratio = 1) +
    ylab("Relative abundance")
```


### Plot Genus 
```r
ps_mock <- subset_samples(ps_idtaxa_samba_kraken2, Individual=="MOCK")

ps_mock = filter_taxa(ps_mock, function(x) sum(x) > 0, TRUE)
ps_mock

ps_agglomerate <- tax_glom(ps_mock, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)
#write.table(taxa, "02_Tables/07_taxa_mock_agglomerate_genus.csv", quote=FALSE, sep=";")

taxa[is.na(taxa)] <- "Non assigned"

taxa %<>% mutate(Kingdom = ifelse(grepl("Unknown ", Kingdom), "unclassified_", Kingdom)) %>% 
   mutate(Phylum = ifelse(grepl("Unknown ", Phylum), sub("Unknown ", "unclassified_", Phylum), Phylum)) %>% 
  mutate(Class = ifelse(grepl("Unknown ", Class), sub("Unknown ", "unclassified_", Class), Class)) %>% 
  mutate(Order = ifelse(grepl("Unknown ", Order), sub("Unknown ", "unclassified_", Order), Order)) %>%  
  mutate(Family = ifelse(grepl("Unknown ", Family),sub("Unknown ", "unclassified_", Family), Family)) %>% 
  mutate(Genus = ifelse(grepl("Unknown ", Genus), sub("Unknown ", "unclassified_", Genus), Genus))

genus <- c("Sulfitobacter", "Cobetia", "Vibrio", "Bacillus", "Ruegeria", "Paraglaciecola", "Leucothrix", "Marinomonas", "Roseovarius", "Pseudoalteromonas",
           "unclassified_Rhodobacteraceae", "unclassified_Alteromonadaceae",
           "unclassified_Vibrionaceae","unclassified_Thiotrichaceae",
           "unclassified_Alteromonadaceae","unclassified_Thiotrichaceae",  
           "unclassified_Bacillaceae",
  "Non assigned")

taxa$Genus[!taxa$Genus %in% genus] <- "Not present"


taxa$Sample <- factor(taxa$Sample, levels = c("MOCK_MARS_BACT_ONT_IDTAXA", "MOCK_MARS_BACT_ONT_SAMBA","MOCK_MARS_BACT_ONT_KRAKEN2"))
taxa$Genus <- factor(taxa$Genus, levels = c("Sulfitobacter", "Cobetia", "Vibrio", "Bacillus", "Ruegeria", "Paraglaciecola", "Leucothrix", "Marinomonas", "Roseovarius", "Pseudoalteromonas",           "unclassified_Rhodobacteraceae",
           "unclassified_Vibrionaceae","unclassified_Thiotrichaceae",
           "unclassified_Alteromonadaceae", "unclassified_Bacillaceae",
            "Not present", "Non assigned"))

genus <- ggplot(data=taxa, aes(x=Sample, y=Abundance, fill=Genus))  +
  ggtitle("B. Mock at genus level")+
  geom_bar(stat="identity", alpha=0.9,  position = "fill") +
  scale_x_discrete(labels=c('IDTAXA', 'SAMBA', 'KRAKEN2'))+
scale_color_manual(values = c("Sulfitobacter" = "#E41A1C", 
                              "Cobetia" = "#377EB8",
                              "Vibrio" = "#FFFF33",
                              "Bacillus" = "#984EA3",
                              "Ruegeria" = "#FF7F00",
                              "Paraglaciecola" = "#4DAF4A",
                              "Leucothrix" = "#A65628", 
                              "Marinomonas" = "#F781BF",
                              "Roseovarius" = "#7570B3",
                              "Pseudoalteromonas" = "#ABDDA4",
                              "unclassified_Rhodobacteraceae" = "#CCCCCC",
                              "unclassified_Alteromonadaceae" = "#CCCCCC",
                              "unclassified_Vibrionaceae" = "#CCCCCC",
                              "unclassified_Thiotrichaceae"= "#CCCCCC",
                              "unclassified_Bacillaceae" = "#CCCCCC",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999")) +
  scale_fill_manual(values = c("Sulfitobacter" = "#E41A1C", 
                              "Cobetia" = "#377EB8",
                              "Vibrio" = "#FFFF33",
                              "Bacillus" = "#984EA3",
                              "Ruegeria" = "#FF7F00",
                              "Paraglaciecola" = "#4DAF4A",
                              "Leucothrix" = "#A65628", 
                              "Marinomonas" = "#F781BF",
                              "Roseovarius" ="#7570B3",
                              "Pseudoalteromonas" = "#ABDDA4",
                              "unclassified_Rhodobacteraceae" = "#CCCCCC",
                              "unclassified_Alteromonadaceae" = "#CCCCCC",
                              "unclassified_Vibrionaceae" = "#CCCCCC",
                              "unclassified_Thiotrichaceae"= "#CCCCCC",
                              "unclassified_Bacillaceae" = "#CCCCCC",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
    legend.text = element_text(size=12),
    axis.title=element_text(size=12),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    legend.title = element_text(size=13), 
    aspect.ratio = 1) +
    ylab("Relative abundance")# +
     # pdf("03_Figures/02_barplot_mock_present_notpresent.pdf", width = 10, height = 10)
```

### Plot Species
```r
ps_mock <- subset_samples(ps_idtaxa_samba_kraken2, Individual=="MOCK")

ps_mock = filter_taxa(ps_mock, function(x) sum(x) > 0, TRUE)
ps_mock

ps_agglomerate <- tax_glom(ps_mock, taxrank = 'Species', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Species <- as.character(taxa$Species)

taxa[is.na(taxa)] <- "Non assigned"

taxa %<>% mutate(Kingdom = ifelse(grepl("Unknown ", Kingdom), "unclassified_", Kingdom)) %>% 
   mutate(Phylum = ifelse(grepl("Unknown ", Phylum), sub("Unknown ", "unclassified_", Phylum), Phylum)) %>% 
  mutate(Class = ifelse(grepl("Unknown ", Class), sub("Unknown ", "unclassified_", Class), Class)) %>% 
  mutate(Order = ifelse(grepl("Unknown ", Order), sub("Unknown ", "unclassified_", Order), Order)) %>%  
  mutate(Family = ifelse(grepl("Unknown ", Family),sub("Unknown ", "unclassified_", Family), Family)) %>% 
  mutate(Genus = ifelse(grepl("Unknown ", Genus), sub("Unknown ", "unclassified_", Genus), Genus)) %>% 
    mutate(Species = ifelse(grepl("Unknown ", Species), sub("Unknown ", "unclassified_", Species), Species))

species <- c("Sulfitobacter_undaria", "Cobetia_litoralis", "Vibrio_splendidus", "Bacillus_safensis", "Ruegeria_meonggei", "Paraglaciecola_mesophila", "Leucothrix_pacifica", "Marinomonas_ushuaiensis", "Roseovarius_nanhaiticus", "Pseudoalteromonas_marina", "unclassified_Paraglaciecola", "unclassified_Bacillus
","unclassified_Sulfitobacter","unclassified_Cobetia","unclassified_Ruegeria","unclassified_Leucothrix","unclassified_Vibrio","unclassified_Pseudoalteromonas","unclassified_Marinomonas","unclassified_Roseovarius","Non assigned")

taxa$Species[!taxa$Species %in% species] <- "Not present"

taxa$Sample <- factor(taxa$Sample, levels = c("MOCK_MARS_BACT_ONT_IDTAXA", "MOCK_MARS_BACT_ONT_SAMBA","MOCK_MARS_BACT_ONT_KRAKEN2"))

taxa$Species <- factor(taxa$Species, levels = c("Sulfitobacter_undaria", "Cobetia_litoralis", "Vibrio_splendidus", "Bacillus_safensis", "Ruegeria_meonggei", "Paraglaciecola_mesophila", "Leucothrix_pacifica", "Marinomonas_ushuaiensis", "Roseovarius_nanhaiticus", "Pseudoalteromonas_marina", "unclassified_Sulfitobacter","unclassified_Cobetia", "unclassified_Vibrio","unclassified_Bacillus","unclassified_Ruegeria","unclassified_Paraglaciecola","unclassified_Leucothrix","unclassified_Marinomonas","unclassified_Roseovarius","unclassified_Pseudoalteromonas","Not present","Non assigned"))
          
species <- ggplot(data=taxa, aes(x=Sample, y=Abundance, fill=Species))  +
  ggtitle("C. Mock at species level")+
  geom_bar(stat="identity", alpha=0.9, position = "fill") +
  scale_x_discrete(labels=c('IDTAXA', 'SAMBA', 'KRAKEN2'))+
scale_color_manual(values = c("Sulfitobacter_undaria" = "#E41A1C", 
                              "Cobetia_litoralis" = "#377EB8",
                              "Vibrio_splendidus" = "#FFFF33",
                              "Bacillus_safensis" = "#984EA3",
                              "Ruegeria_meonggei" = "#FF7F00",
                              "Paraglaciecola_mesophila" = "#4DAF4A",
                              "Leucothrix_pacifica" = "#A65628", 
                              "Marinomonas_ushuaiensis" = "#F781BF",
                              "Roseovarius_nanhaiticus" = "#7570B3",
                              "Pseudoalteromonas_marina" = "#ABDDA4",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999",
                              "unclassified_Sulfitobacter" = "#CCCCCC", 
                              "unclassified_Cobetia" = "#CCCCCC",
                              "unclassified_Vibrio" = "#CCCCCC",
                              "unclassified_Bacillus" = "#CCCCCC",
                              "unclassified_Ruegeria" = "#CCCCCC",
                              "unclassified_Paraglaciecola" = "#CCCCCC",
                              "unclassified_Leucothrix" = "#CCCCCC", 
                              "unclassified_Marinomonas" = "#CCCCCC",
                              "unclassified_Roseovarius" = "#CCCCCC",
                              "unclassified_Pseudoalteromonas" = "#CCCCCC")) +
  scale_fill_manual(values = c("Sulfitobacter_undaria" = "#E41A1C", 
                              "Cobetia_litoralis" = "#377EB8",
                              "Vibrio_splendidus" = "#FFFF33",
                              "Bacillus_safensis" = "#984EA3",
                              "Ruegeria_meonggei" = "#FF7F00",
                              "Paraglaciecola" = "#4DAF4A",
                              "Leucothrix_pacifica" = "#A65628", 
                              "Marinomonas_ushuaiensis" = "#F781BF",
                              "Roseovarius_nanhaiticus" ="#7570B3",
                              "Pseudoalteromonas_marina" = "#ABDDA4",
                              "unclassified_Sulfitobacter" = "#CCCCCC", 
                              "unclassified_Cobetia" = "#CCCCCC",
                              "unclassified_Vibrio" = "#CCCCCC",
                              "unclassified_Bacillus" = "#CCCCCC",
                              "unclassified_Ruegeria" = "#CCCCCC",
                              "unclassified_Paraglaciecola" = "#CCCCCC",
                              "unclassified_Leucothrix" = "#CCCCCC", 
                              "unclassified_Marinomonas" = "#CCCCCC",
                              "unclassified_Roseovarius" = "#CCCCCC",
                              "unclassified_Pseudoalteromonas" = "#CCCCCC",
                              "Not present" = "#B10026",
                              "Non assigned" = "#999999"))+
   theme_bw()+
  theme(axis.title.x = element_blank(),
    legend.text = element_text(size=12),
    axis.title=element_text(size=12),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    legend.title = element_text(size=13), 
    aspect.ratio = 1 ) +
    ylab("Relative abundance")#+
    #  pdf("03_Figures/03_barplot_mock_species.pdf", width = 10, height = 10)

```


# Bacterial comparisons on entire dataset
## 1. Libraries
```r
library(phyloseq)
library(tidyverse)
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
library(ranacapa)
```

## 2. Generate one ps per all sequencing methods
### Novaseq 
```r
otu_ns <- ps_novaseq@otu_table %>% as.data.frame() #6996 ASVS, 24 samples
tax_ns <- ps_novaseq@tax_table %>% as.data.frame()

write.table(otu_ns, "01_Tables/otu_ns.csv", sep=";", quote=FALSE)
write.table(tax_ns, "01_Tables/tax_ns.csv", sep=";", quote=FALSE)
```

### Miseq
```r
otu_ms <- ps_miseq@otu_table %>% as.data.frame() #4252 ASVs, 24 samples
tax_ms <- ps_miseq@tax_table %>% as.data.frame()

write.table(otu_ms, "01_Tables/otu_ms.csv", sep=";", quote=FALSE)
write.table(tax_ms, "01_Tables/tax_ms.csv", sep=";", quote=FALSE)
```

### ONT
```r
otu_ont <- ps_ont@otu_table %>% as.data.frame() #739 ASVs
tax_ont <- ps_ont@tax_table %>% as.data.frame()

write.table(otu_ont, "01_Tables/otu_ont.csv", sep=";", quote=FALSE)
write.table(tax_ont, "01_Tables/tax_ont.csv", sep=";", quote=FALSE)
```

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
ps_ms_ns_ont

ps_ms_ns_ont = filter_taxa(ps_ms_ns_ont, function(x) sum(x) > 0, TRUE)
ps_ms_ns_ont
```

## SAVE
```r
ps_ms_ns_ont <- subset_samples(ps_ms_ns_ont, !ID=="L2I4_MARS")
ps_ms_ns_ont = filter_taxa(ps_ms_ns_ont, function(x) sum(x) > 0, TRUE)

saveRDS(ps_ms_ns_ont, "00_Phyloseq_objects/01_ps_ms_ns_ont_without_L1I4_MARS.rds")
```

## 2. Normalization
```r
total = median(sample_sums(ps_ms_ns_ont)) # Median= 36288
standf = function(x, t=total) round(t * (x / sum(x))) # Standardize abundances to the median sequencing depth
ps_filt = transform_sample_counts(ps_ms_ns_ont, standf)

ps_ms_ns_ont_trans = filter_taxa(ps_filt, function(x) sum(x > total* 0.0003) > 0, TRUE)
ps_ms_ns_ont_trans
```

## SAVE
```r
saveRDS(ps_ms_ns_ont_trans, "00_Phyloseq_objects/02_ps_ms_ns_ont_trans_without_L2I4_MARS.rds")
```

```r
ps_agglomerate <- tax_glom(ps_ms_ns_ont_trans, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)
write.table(taxa, "01_Tables/04_taxa_agglomerate_genus.csv", sep=";", quote=FALSE)

#Clean table by hand to concatenate with ;NA when unassigned rank
taxa <- read.table("01_Tables/04_taxa_agglomerate_genus.csv", sep=";", header=TRUE, dec=".")

taxa$Sequencing_type <- factor(taxa$Sequencing_type, levels = c("MISEQ", "NOVASEQ", "ONT"))
taxa$Month <- factor(taxa$Month, levels = c("November", "March"))

taxa %<>% filter(!Genus=="Ascophyllum")
class <- c("Alphaproteobacteria","Gammaproteobacteria", "Cyanobacteriia","Planctomycetes","Verrucomicrobiae","Bacteroidia", "Phycisphaerae")

taxa$Class[!taxa$Class %in% class] <- "Others microorganisms"

palette_class <- c(
                "Others microorganisms"="grey",
                "Phycisphaerae"="#FFADC1",
                "Verrucomicrobiae"="#EAAF29",
                "Planctomycetes"="#645A9F",
                  "Bacteroidia"="#853250",
                  "Cyanobacteriia"="#00441B",
                  "Gammaproteobacteria"="#B03F00",
                "Alphaproteobacteria"= "#08306B")

taxa$Class <- factor(taxa$Class, levels = c('Others microorganisms',
                                          "Phycisphaerae",
                                          "Bacteroidia",
                                          "Verrucomicrobiae",
                                          "Planctomycetes",
                                          "Cyanobacteriia",
                                          "Gammaproteobacteria",
                                          "Alphaproteobacteria"))


barplot <- ggplot(data=taxa, aes(x=Algae, y=Abundance, fill=Class))  +
  ggtitle("D. Sequencing comparison at class level")+
  facet_nested(Sequencing_type ~ Month+Site, space = "free", scales = "free") +
      scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity", position = "fill") +
  theme(axis.line = element_line(colour = "black"),
        strip.text=element_text(size=13),
        legend.title = element_text(size=13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
     legend.text = element_text(size=13),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"),
    axis.text.y = element_text(size = 10, color="black"),
        axis.title.x = element_blank()) +
   ylab("Relative abundance") 
```



# 3. Combined plots
```r
x <- gridExtra::grid.arrange(order, genus, species, nrow=3)
z <- gridExtra::grid.arrange(x, barplot, nrow=1)
ggsave("04_combined_mock_barplot_analysis.pdf", z, height = 15, width = 25)
```

![Figure 5 | Mock and bacterial comparisons.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/04_combined_mock_barplot_analysis.png)<!-- -->
