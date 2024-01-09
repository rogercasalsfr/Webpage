---
title: "Pre-Processing Seurat"
author: "Roger Casals"
date: 2020-12-01T21:13:14-05:00
categories: ["R"]
tags: ["R Markdown", "single cell RNA", "Seurat"]
---

{{< icon name="download" pack="fas" >}} Download the HTML {{< staticref "uploads/index.en.html" "newtab" >}}CV{{< /staticref >}}.

# R Markdown



 Primer de tot, carreguem les dades. Del SCC

# Squamos cell carcionma 

```r

#setwd("C:/Users/Roger Casals/OneDrive/Escriptori/UOC ROGER/2n semestre/TFM/Final 5_5")
setwd("C:/Users/Natàlia/Desktop/UOC ROGER/2n semestre/TFM/Primera prova")
metadata2 <- data.table::fread("GSE123813_scc_metadata.txt.gz")
expression_data <- read.table("GSE123813_scc_scRNA_counts.txt.gz", header=TRUE, row.names=1, sep="\t", check.names=FALSE)

```




## Càrrega de llibreries

Carreguem les diferents llibreries

```r
library(Seurat)
library(SeuratWrappers)
library(monocle3)
#library(slingshot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(dyno)
library(tidyverse)
library(dynwrap)
#library(celldex)
library(SingleCellExperiment)
#library(SingleR)
library(dyneval)

```




## Preprocessament amb Seurat

### Objecte Seurat

```r


aa <- CreateSeuratObject(counts=expression_data, project="GSE123813", metadata=metadata2)

bbb <- aa



#Ordenem les dades per tal d'obtenir les identificacions correctes.

metadata2 <- metadata2[order(metadata2$cell.id),]
aa <- aa[, order(colnames(aa))]
aa@meta.data <- aa@meta.data[order(row.names(aa@meta.data)), ]



```




### Afegim metadata a Seurat

```r

aa@meta.data$cluster <- metadata2$cluster
aa@meta.data$treatment <- metadata2$treatment


table(aa@meta.data$cluster)

```




ENS CENTRAREM AMB LES CD4 perquè hem vist a la teoria que les "naive" es poden diferenciar en Treg o Th17 o Tfh


### Passos preprocessament amb Seurat


#### Control de qualitat

```r

#Obtenim els percentatges de gens mitocontrials i ribosòmics, per saber si les cèl·lules estan estressades.
aa <- PercentageFeatureSet(aa, "^MT-", col.name="percent_mito")
aa <- PercentageFeatureSet(aa, "^RP[SL]", col.name = "percent_ribo")
aa <- PercentageFeatureSet(aa, "^HB[^(P)]", col.name = "percent_hb")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

VlnPlot(aa, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()

FeatureScatter(aa, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)


CD4tots <- aa

pre.markers <- FindAllMarkers(CD8tots ,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


pre.markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) -> top15
DoHeatmap(CD8tots, features = top15$gene, size=2)



unique(aa@meta.data$cluster)

CD8tots <- subset(aa, cluster %in% c("CD8_ex", "CD8_mem", "CD8_eff", "CD8_naive", "CD8_ex_act", "CD8_act"))

CD4tots <- CD8tots

CD4tots <- subset(CD4tots, subset=nFeature_RNA < 4000 & percent_mito < 8 )

table(CD4tots@meta.data$cluster)


#Normalitzem les dades
CD4tots <- NormalizeData(CD4tots)

#Busquem els gens que contribueixen més significativament a la variabilitat d'una matriu d'expressió gènica.
#En seleccionem 2000 per defecte.
CD4tots <- FindVariableFeatures(CD4tots, selection.method = "vst", nfeatures = 2000)


#Escalem les dades mitjançant les funció que ens proporciona Seurat.
CD4tots <- ScaleData(CD4tots)

#Reduïm dimensió amb el RunPCA  (Utilitzem el número 16 ja que l'utilitzen al experiment.)
CD4tots <- RunPCA(CD4tots, npcs = 16)

                  #features = VariableFeatures(object = CD4tots))
#Utilitzem aquesta funció per trobar veïns propers mitjançant l'expressió gènica
CD4tots <- FindNeighbors(CD4tots, dims = 1:16)



#Ens quedem amb el cluster 0.3, ja que és el que millor s'adapta als subgrups que tenimm.


CD4tots <- FindClusters(CD4tots, resolution = c(0.7))


CD4tots <- RunUMAP(CD4tots, dims=1:16)



#Realitzem un DimPlot mitjançant per tal de que sigui etiquetat de diverses formes.
DimPlot(CD4tots, reduction = "umap", group.by = "treatment")

#Realitzem el plot per veure com estan ubicades les cèl·lules anotades en l'espai.
DimPlot(CD4tots, reduction="umap", group.by="cluster", label=T)


#Escollim quin serà el millor cluster.
#DimPlot(CD4tots, reduction = "umap", group.by = "RNA_snn_res.0.9", label = TRUE)
#DimPlot(CD4tots, reduction = "umap", group.by = "RNA_snn_res.0.7", label = TRUE)
#DimPlot(CD4tots, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(CD4tots, reduction = "umap", group.by = "RNA_snn_res.0.7", label = TRUE)
#DimPlot(CD4tots, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)



```




S'ha escollit el cluster 0.3, ja que és el mateix que utilitzen a l'article, i el que engloba millor els diferents subtipus que tenim.

## Rename idents

```r

CD4tots <- RenameIdents(CD4tots, "2" = "CD8_Naive")
CD4tots <- RenameIdents(CD4tots, "7" = "CD8_Naive")


CD4tots <- RenameIdents(CD4tots, "0" = "CD8_mem")
CD4tots <- RenameIdents(CD4tots, "3" = "CD8_mem")


CD4tots <- RenameIdents(CD4tots, "1" = "CD8_ex")
CD4tots <- RenameIdents(CD4tots, "8" = "CD8_ex")
CD4tots <- RenameIdents(CD4tots, "4" = "CD8_ex")

CD4tots <- RenameIdents(CD4tots, "5" = "CD8_act")
CD4tots <- RenameIdents(CD4tots, "11" = "CD8_act")

CD4tots <- RenameIdents(CD4tots, "9" = "CD8_ex_act")

CD4tots <- RenameIdents(CD4tots, "6" = "CD8_eff")

#El 10 no el posem

```
```r

CD4tots <- RenameIdents(CD4tots, "5" = "CD8_Naive")
CD4tots <- RenameIdents(CD4tots, "14" = "CD8_Naive")


CD4tots <- RenameIdents(CD4tots, "0" = "CD8_mem")
CD4tots <- RenameIdents(CD4tots, "8" = "CD8_mem")


CD4tots <- RenameIdents(CD4tots, "1" = "Naive")
CD4tots <- RenameIdents(CD4tots, "3" = "CD8_ex")
CD4tots <- RenameIdents(CD4tots, "9" = "CD8_ex")
CD4tots <- RenameIdents(CD4tots, "13" = "CD8_ex")

CD4tots <- RenameIdents(CD4tots, "7" = "CD8_act")
#CD4tots <- RenameIdents(CD4tots, "11" = "CD8_act")

CD4tots <- RenameIdents(CD4tots, "15" = "CD8_ex_act")

CD4tots <- RenameIdents(CD4tots, "12" = "CD8_eff")

CD4tots <- RenameIdents(CD4tots, "2" = "Tfh")
CD4tots <- RenameIdents(CD4tots, "16" = "Tfh")


CD4tots <- RenameIdents(CD4tots, "4" = "Treg")
CD4tots <- RenameIdents(CD4tots, "10" = "Treg")
CD4tots <- RenameIdents(CD4tots, "17" = "Treg")

CD4tots <- RenameIdents(CD4tots, "6" = "Th17")


```



Creem una nova columna a la metadata amb els "idents" d'aquests clusters

```r
CD4tots@meta.data$idents <- Idents(CD4tots)

```




### Només grups interessants
Realitzem un subset, identificant només els grups que hem pogut nombrar amb els clusters, i els representem per a veure'ls.

```r
CD8tots <- subset(CD4tots, idents %in% c("CD8_ex", "CD8_mem", "CD8_eff", "CD8_Naive", "CD8_ex_act", "CD8_act", "Naive", "Treg", "Th17", "Tfh"))
```


```r
CD8tots <- subset(CD4tots, cluster %in% c("CD8_ex", "CD8_mem", "CD8_eff", "CD8_naive", "CD8_ex_act", "CD8_act"))

CD4totssubset <- CD8tots
CD4totssubset <- subset(CD4totssubset, idents %in% c("CD8_ex", "CD8_mem", "CD8_eff", "CD8_Naive", "CD8_ex_act", "CD8_act"))

plot1 <- DimPlot(CD4totssubset, reduction = "umap", group.by="treatment")

plot2 <- DimPlot(CD4totssubset, reduction = "umap", group.by="idents", label=T)


plot1 | plot2






```


### Identificats segons condició

Identifiquem els nous grups segons la seva condició de tractament, és a dir, si són "pre" o "post", els separem mitjançant un guió baix.


```r

CD4totssubset$celltype.cnd <- paste0(CD4totssubset$treatment, "_", CD4totssubset$idents)

DimPlot(CD4totssubset, reduction = "umap", group.by="celltype.cnd", label=T)

```



Identifiquem els Idents del cluster, com al que hem definit.

```r

Idents(CD4totssubset) <- CD4totssubset$celltype.cnd


```


### Proporcions cel·lulars

Calculem el percentatge de cada cèl·lula

```r

prop.table(table(CD4totssubset@meta.data$celltype.cnd))*100


table(CD4totssubset@meta.data$celltype.cnd)

```
```r


CD4totssubset@meta.data$cluster_redefined <- CD4totssubset@meta.data$cluster



```



```r

CD8totssubset <- CD4totssubset


CD8totssubsetpre <- subset(CD8totssubset, celltype.cnd %in% c("pre_CD8_act", "pre_CD8_eff", "pre_CD8_ex", "pre_CD8_ex_act", "pre_CD8_mem", "pre_CD8_Naive"))
CD8totssubsetpost <- subset(CD8totssubset, celltype.cnd %in% c("post_CD8_act", "post_CD8_eff", "post_CD8_ex", "post_CD8_ex_act", "post_CD8_mem", "post_CD8_Naive"))




library(Seurat)


library(SeuratDisk)
library(SeuratData)

#SaveH5Seurat(CD8totssubset, "totscd8", overwrite=FALSE, verbose=TRUE)
#SaveH5Seurat(CD8totssubsetpre, "precd8", overwrite=FALSE, verbose=TRUE)
#SaveH5Seurat(CD8totssubsetpost, "postcd8", overwrite=FALSE, verbose=TRUE)


```




