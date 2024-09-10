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
```

### 2. Maps
```{r}
register_stadiamaps("YOUR_API_KEY", write = TRUE)

France <- get_stadiamap( bbox = c(left = -6, bottom = 41, right = 10, top = 52), zoom = 6, maptype = "stamen_watercolor")
ggmap(France) + 
  theme_void()+
  theme(
    plot.title = element_text(colour = "orange"), 
    panel.border = element_rect(colour = "grey", fill=NA, size=2)
  )+ pdf("Figures/01_France.pdf", height = 5, width = 5)

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

![Figure 1 | Sampling locations.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/01_Sampling_map.png)<!-- -->



## B. Read sequencing characteristics 

```r
library(ggplot2)
library(tidyverse)
library(magrittr)
library(ggh4x)
```

```r
char <- read.table("03_tables/02_sequencing_characteristics_keep_only_ss_IDTAXA.csv", sep=";", dec =".", header=TRUE)


char$Sequencing_type <- factor(char$Sequencing_type, levels = c("MISEQ", "NOVASEQ", "SAMBA", "IDTAXA", "NGSpeciesID", "WITHOUT ADAPTATIVE SAMPLING","ADAPTATIVE SAMPLING")) 
char$Quality_filtering_steps <- factor(char$Quality_filtering_steps, levels = c("Input", "Keep after first filtering", "Keep for final analysis after second filtering")) 
char$Target <- factor(char$Target, levels = c("Evaluation of adaptative sampling","Sequencing comparisons 16S","Sequencing comparisons ITS")) 
 
 ggplot(char, aes(y = Percent, x= Sequencing_type, fill=Quality_filtering_steps)) +
        geom_col(position = position_dodge()) +
   facet_nested(~Target,  space = "free", scales = "free")+
     theme_bw()+
   scale_color_manual(values=c("Input"="#1F78B4","Keep after first filtering"="#FDAE61", "Keep for final analysis after second filtering"= "#E31A1C"))+   
   scale_fill_manual(values=c("Input"="#1F78B4","Keep after first filtering"="#FDAE61", "Keep for final analysis after second filtering"= "#E31A1C"))+
        xlab("Samples")   +
   theme(strip.text.y = element_text(angle = 0))+
   pdf("04_Figures/01_Sequencing_characteristics2.pdf", height = 5, width = 10)
```



