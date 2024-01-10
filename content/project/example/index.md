---
date: "2016-04-27T00:00:00Z"
external_link: ""
image:
  caption: Photo by rawpixel on Unsplash
  focal_point: Smart
links:
- icon: twitter
  icon_pack: fab
  name: Follow
  url: 
#slides: {{< icon name="download" pack="fas" >}} Download the HTML {{< staticref "uploads/index.en.html" "newtab" >}}Pre-Processing Seurat{{< /staticref >}}. 
summary: An example of using the in-built project page.
tags:
- Master thesis
title: Comparison of Single Cell trajectory inference methods in cancer immunotherapy
url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""
---

The main goal of this project is to discern immune cell trajectories and gene regulatory networks within a tumor microenvironment. The project's efficacy is under assessment. Another aim is to identify differentially expressed genes to delineate the trajectories and conduct enrichment analysis.

The data utilized was sourced from the GSE123813 repository, exclusively employing samples from SCC patients. The data is from "pre" and "post" immunotherapy treatment

The code is housed in a folder and encompasses the preprocessing steps integral to the Seurat workflow. A comparative analysis is conducted on the trajectory inference methods of Dynverse and Monocle3. Enrichment analysis is carried out using ClusterProfiler, integrated into the Monocle3 files.
