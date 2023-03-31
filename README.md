# MiMed

Title: Comprehensive microbiome causal mediation analysis using MiMed on user-friendly web interfaces

Version: 1.0.0

Maintainer: Hyo Jung Jang <hyojung.jang@stonybrook.edu>

Description: MiMed is a unified web cloud platform for microbiome mediation analysis that probes for the roles of the microbiome as a mediator that transmits the effect of a treatment (e.g., environmental, behavioral or medical exposures) to an outcome (e.g., health or disease status). The main features of MiMed are as follows. First, MiMed can survey the microbiome in various spheres (1) as a whole microbial ecosystem using different ecological measures (e.g., alpha- and beta-diversity indices) (referred as Community-level Analysis) or (2) as individual microbial taxa (e.g., phyla, classes, orders, families, genera, species) using different data normalization methods (referred as Taxonomy-level Analysis). Second, MiMed enables covariate-adjusted analysis to control for potential confounding factors (e.g., age, gender), which is essential to enhance the causality of the results especially for observational studies. Third, MiMed enables a breadth of statistical inferences in both mediation effect estimation and significance testing. Finally, MiMed provides flexible and easy-to-use data processing and analytic modules and makes nice graphical representations.

NeedsCompilation: No

Depends: R(≥ 4.1.0)

Imports: Bioconductor ('BiocParallel', 'biomformat', 'phyloseq'); CRAN ('betareg', 'BiasedUrn', 'BiocManager', 'bios2mds', 'CompQuadForm', 'dashboardthemes', 'devtools', 'DiagrammeR', 'dirmult', 'dplyr', 'DT', 'ecodist', 'entropart', 'erer', 'fBasics', 'forestplot', 'fossil', 'ggplot2', 'ggthemes', 'googleVis', 'gridExtra', 'gridGraphics', 'gridExtra', 'compositions', 'GUniFrac', 'htmltools', 'ICSNP', 'lme4', 'lmerTest', 'MiRKAT', 'nlme', 'patchwork', 'phangorn', 'picante', 'plotly', 'PMCMRplus', 'quantreg', 'remotes', 'reticulate', 'rgl', 'rmarkdown', 'robCompositions', 'robustbase', 'seqinr', 'shiny', 'shinydashboard', 'shinyjs', 'shinyWidgets', 'stringr', 'tidyverse', 'vegan', 'xtable', 'zCompositions', 'zip', ); GitHub ('LDM', 'volcano3D')

License: GPL 1, GPL 2 

## URLs

* Web application (online implementation): 
* GitHub repository (local implementation): https://github.com/yj7599/MiMedGit

## References

* Jang H, Park S, Koh H. (2023) Comprehensive microbiome causal mediation analysis using MiMed on user-friendly web interfaces (*In review*). 


## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 

shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiMedGit', 'yj7599', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiPair, please report in Issues (https://github.com/YJ7599/MiMedGit/issues) or email Hyo Jung Jang (hyojung.jang@stonybrook.edu).

* Tip 1. For the local implementation, to ensure no fail builds, you may need to remove all the pre-installed R packages on your computer or start it all over again after re-installing R and R Studio completely.
