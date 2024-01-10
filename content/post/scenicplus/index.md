---
authors: Roger Casals Franch
- admin
categories: []
date: "2019-02-05T00:00:00Z"
image:
  caption: ""
  focal_point: ""
lastMod: "2019-09-05T00:00:00Z"
projects: []
subtitle: Last step of scenic+
summary: We find the regulons
tags: []
title: Scenic+ find regulons
---


```python
from IPython.core.display import Image
Image('https://www.python.org/static/community_logos/python-logo-master-v3-TM-flattened.png')
```


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



```python

run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_celltype'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '/home/roger/Github/OpenProblems/Conda/scenicplus/TF_names_v_1.01.txt',
        save_path = os.path.join('/home/roger/Github/OpenProblems/Conda/scenicplus/Tutorial'),
        #biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '/home/roger/Github/OpenProblems/Conda/scenicplus/',
        n_cpu = 10,
        _temp_dir = '/home/roger/tempdir_rayspill')
        
```

