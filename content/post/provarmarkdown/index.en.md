---
title: "Pre-Processing Seurat"
author: "Roger Casals"
date: 2020-12-01T21:13:14-05:00
categories: ["R"]
tags: ["R Markdown", "single cell RNA", "Seurat"]
---



# R Markdown



 Primer de tot, carreguem les dades. Del SCC

# Squamos cell carcionma 

```{r}

#setwd("C:/Users/Roger Casals/OneDrive/Escriptori/UOC ROGER/2n semestre/TFM/Final 5_5")
setwd("C:/Users/Natàlia/Desktop/UOC ROGER/2n semestre/TFM/Primera prova")
metadata2 <- data.table::fread("GSE123813_scc_metadata.txt.gz")
expression_data <- read.table("GSE123813_scc_scRNA_counts.txt.gz", header=TRUE, row.names=1, sep="\t", check.names=FALSE)

```


You can also embed plots. See Figure <a href="#fig:pie">1</a> for example:


```r
par(mar = c(0, 1, 0, 1))
pie(
  c(280, 60, 20),
  c('Sky', 'Sunny side of pyramid', 'Shady side of pyramid'),
  col = c('#0292D8', '#F7EA39', '#C4B632'),
  init.angle = -50, border = NA
)
```

<div class="figure">
<img src="{{< blogdown/postref >}}index.en_files/figure-html/pie-1.png" alt="A fancy pie chart." width="672" />
<p class="caption"><span id="fig:pie"></span>Figure 1: A fancy pie chart.</p>
</div>
