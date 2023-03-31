# MiMed

Title: MiMed: Comprehensive microbiome causal mediation analysis on user-friendly web interfaces

Version: 1.0.0

Maintainer: Hyo Jung Jang <hyojung.jang@stonybrook.edu>

Description: MiPair is an integrative web cloud service for design-based comparative analysis with paired microbiome data. Pairing (or blocking) is a design technique that is widely used in comparative microbiome studies to efficiently control for the effects of potential confounders (e.g., genetic, environmental, or behavioral factors). Some typical paired (block) designs for human microbiome studies are repeated measures designs that profile each subject's microbiome twice (or more than twice) 1) for pre and post treatments to see the effects of a treatment on microbiome, or 2) for different organs of the body (e.g., gut, mouth, skin) to see the disparity in microbiome between (or across) organs. MiPair enables comprehensive comparative analysis in sequence for such paired microbiome studies on user-friendly web environments. Detailed features are as follows.

* A variety of data uploading, quality controlling, analytic and graphical procedures that produce publishable data, tables, and plots
* Comparative analysis between (or across) groups
* Comparative analysis between baseline (or reference) and other groups
* Parametric or non-parametric tests for incomplete or complete block designs
* Both ecological (alpha- and beta-diversity) and taxonomic (e.g., phylum, class, order, family, genus, species) analysis

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'gridExtra', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip'); GitHub ('LDM', 'volcano3D')

License: GPL 1, GPL 2 

## URLs

* Web application (online implementation): http://mipair.micloud.kr
* GitHub repository (local implementation): https://github.com/yj7599/MiPairGit

## References

* Jang H, Koh H, Gu W, Kang B. (2022) Integrative web cloud computing and analytics using MiPair for design-based comparative analysis with paired microbiome data. Scientific Reports 12:20465 https://doi.org/10.1038/s41598-022-25093-6.

## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiPairGit', 'yj7599', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiPair, please report in Issues (https://github.com/YJ7599/MiPairGit/issues) or email Hyo Jung Jang (hyojung.jang@stonybrook.edu).

* Tip 1. For the local implementation, to ensure no fail builds, you may need to remove all the pre-installed R packages on your computer or start it all over again after re-installing R and R Studio completely.
