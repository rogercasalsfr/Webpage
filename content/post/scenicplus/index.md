---
authors:
- admin
categories: []
date: "2019-02-05T00:00:00Z"
image:
  caption: ""
  focal_point: ""
lastMod: "2019-09-05T00:00:00Z"
projects: []
subtitle: inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
summary: inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
tags: []
title: SCENIC+
---

```python

#Install and import libraries

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pyranges
import sys
# Set stderr to null to avoid strange messages from ray
_stderr = sys.stderr

import pickle
from scenicplus.wrappers.run_scenicplus import run_scenicplus
null = open(os.devnull,'wb')


#Change directories
#!wget -O pbmc_tutorial/data/utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
#!wget -O pbmc_tutorial/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
#!chmod +x pbmc_tutorial/bedToBigBed


biomart_host = "http://sep2019.archive.ensembl.org/"


#with open('/home/roger/Github/OpenProblems/Conda/scenicplus/scplus_obj.pkl', 'rb') as f:
   # scplus_obj = pickle.load(f)

scplus_obj = pd.read_pickle("/home/roger/Github/OpenProblems/Conda/scenicplus/scplus_obj.pkl")


```


