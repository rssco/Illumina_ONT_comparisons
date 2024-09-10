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

![Figure 1 | Sampling locations.](https://github.com/rssco/Illumina_ONT_comparisons/blob/main/01_Figures/01_Methods_sampling.png)<!-- -->
