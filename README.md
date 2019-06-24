[![Build Status](https://travis-ci.org/NCBI-Hackathons/HLAClustRView.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/HLAClustRView)
[![codecov](https://codecov.io/gh/NCBI-Hackathons/HLAClustRView/branch/master/graph/badge.svg)](https://codecov.io/gh/NCBI-Hackathons/HLAClustRView)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
<img src="/vignettes/HLAClustRViewLogo.png" alt=""/>
</p>

# HLA typing clustering and visualization based on specific similarity metrics #


The **HLAClustRView** package implements specialized similarity metrics that
quantify the similarity between human leukocyte antigen (HLA) typing from 
multiple samples. Using these metrics, the package enables hierarchical 
cluster analysis and can bind RNA-seq expression (when available) to the 
clusters.

## Goal ##

The HLA genes plays a primary biological role in 
the regulation of the immune system and in the outcomes of human organ 
transplantation. However, the official HLA typing format is tedious to process.
So far, not much analysis has been done from the clustering perspective for
this type of data.
The goal of the **HLAClustRView** package is to implement appropriate 
similarity metrics on HLA typing to enable sample clustering and more complex
analysis.

## Description of Package Functionality ##

The vignette of the **HLAClustRView** package, which is a document that 
provides a task-oriented description of the package functionality, contains the 
most up-to-date information.

## Package Details ##

* **Version:** 0.99.0
* **Depends:** R (>= 3.3)
* **Imports:** purrr, dplyr, tidyr, tibble, utils, methods, stringr, data.table, stats, ComplexHeatmap, graphics, rlang
* **Suggests:** testthat, knitr, rmarkdown, circlize, BiocStyle
    
## Package Workflow ##


<p align="center">
<img src="/vignettes/HLAdesign.jpg" alt=""/>
</>

## Citing ##

If you use the *HLAClustRView* package 
for a publication, we would ask you to cite the following:

> Pascal Belleau, Adewunmi Adelaja, Astrid Deschênes, Santiago Medina, Nissim Ranade and Allissa Dillman (2018). HLAClustRView: HLA typing clustering and visualization based on specific similarity metrics. R package version 0.99.0.

## Authors ##

[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau"),
[Adewunmi Adelaja](https://www.linkedin.com/in/adewunmi-adelaja-2b3b1635/ "Adewunmi Adelaja"), 
[Astrid Deschênes](https://www.linkedin.com/in/astriddeschenes "Astrid Deschênes"),
[Santiago Medina](https://github.com/santiago1234), [Nissim Ranade](https://www.linkedin.com/in/nissim-ranade-4029b3b5 "Nissim Ranade") and Allissa Dillman

## License ##

This package and the underlying *HLAClustRView* code are distributed under 
the MIT license. You are free to use and redistribute this software. 

For more information on MIT license see: [https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT)

## Bugs/Feature requests ##

If you have any bugs or feature requests, 
[let us know](https://github.com/NCBI-Hackathons/Integrating-HLA-typing-methods-and-RNA-seq/issues). 

Thanks!
