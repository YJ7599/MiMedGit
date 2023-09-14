library(biomformat)
library(phangorn)
#library(bios2mds)


otu.tab.clean.no.tree <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  #tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  biom <- merge_phyloseq(otu.tab, tax.tab, sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}

otu.tab.clean.no.tax.tab <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  biom <- merge_phyloseq(otu.tab, tree, sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}



otu.tab.clean.no.both <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  biom <- merge_phyloseq(otu.tab,  sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}


preprocess.tax.tab = function(tax.tab){
  trans.tax.tab <- matrix(NA, nrow(tax.tab), 7)
  tax.list <- strsplit(as.character(tax.tab$Taxon), ";")
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_0__", "", tax.list[[i]][grepl("D_0__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,1] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_1__", "", tax.list[[i]][grepl("D_1__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,2] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_2__", "", tax.list[[i]][grepl("D_2__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,3] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_3__", "", tax.list[[i]][grepl("D_3__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,4] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_4__", "", tax.list[[i]][grepl("D_4__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,5] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_5__", "", tax.list[[i]][grepl("D_5__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,6] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_6__", "", tax.list[[i]][grepl("D_6__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,7] <- taxon
    }
  }
  rownames(trans.tax.tab) <- tax.tab$Feature.ID
  colnames(trans.tax.tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  return(trans.tax.tab)
}

biom.check.samples <- function(otu.table, sample.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  sam.dat <- sample_data(sample.data)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  
  if (length(intersect(colnames(otu.tab), rownames(sam.dat))) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#h4("Error: There is no common samples among OTU table and Sample Data")
#h4("Error: There is no common OTUs among OTU table, taxonomic table and tree tip labels")
biom.check.otu <- function(otu.table, tax.table, tre.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  tax.tab <- tax_table(tax.table)
  tree <- phy_tree(tre.data)
  
  if (length(intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))) == 0) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}



biom.clean.no.tree <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  #tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(rownames(otu.tab), rownames(tax.tab)) 
  
  
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  #tree <- prune_taxa(ind.com.otu, tree)
  # if(!is.rooted(tree)) {
  #   tree <- phangorn::midpoint(tree)
  # }
  
  if (tax.tab.c) {
    tax.tab <- tax.tab.clean(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  
  biom <- merge_phyloseq(otu.tab, tax.tab,  sam.dat)
  
  biom <- otu.tab.clean.no.tree(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

biom.clean.no.both <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  
  otu.tab <- otu_table(biom)
  #tax.tab <- tax_table(biom)
  #tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  #ind.com.otu <- intersect(rownames(otu.tab), rownames(tax.tab))

  #if (length(ind.com.otu) == 0) {
  #  stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  #}
  
  #ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  #otu.tab <- otu.tab[ind.com.1,]
  #ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  #tax.tab <- tax.tab[ind.com.2,]
  #tree <- prune_taxa(ind.com.otu, tree)
  
  biom <- merge_phyloseq(otu.tab,  sam.dat)
  
  biom <- otu.tab.clean.no.both(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

biom.clean.no.tax.tab <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  #tax.tab <- tax_table(biom)
  
  # if (kingdom != "all") {
  #   
  #   ind <- is.element(tax.tab[,1], kingdom)
  #   rownames(tax.tab)[ind]
  #   biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  # }
  
  otu.tab <- otu_table(biom)
  #tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(rownames(otu.tab), tree$tip.label)
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  #ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  #tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  
  
  biom <- merge_phyloseq(otu.tab, tree, sam.dat)
  
  biom <- otu.tab.clean.no.tax.tab(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

biom.check.otu.without.tree <- function(otu.table, tax.table) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  tax.tab <- tax_table(tax.table)
  
  if (length(intersect(rownames(otu.tab), rownames(tax.tab))) == 0) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

biom.check.otu.without.tax.tab <- function(otu.table, tre.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  tree <- phy_tree(tre.data)
  
  if (length(intersect(rownames(otu.tab), tree$tip.label)) == 0) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}


#methods
REFERENCE_CHECK_M <- function(data_transform = "", method_name = "", FDR = ""){
  reference_lists <- NULL
  
  ## Data_transform
  if(data_transform == "CLR (Default)"){
    reference_lists = c(reference_lists, "Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-60.")
  } else if (data_transform == "CLR"){
    reference_lists = c(reference_lists, "Aitchison J. The statistical analysis of compositional data. J R Stat Soc Series B Stat Methodol. 1982;44(2):139-60.")
  }
  
  ## Method_name
  if (method_name == "Imai Method (Default)"){
    reference_lists = c(reference_lists, "Imai K., Keele L. & Tingley D. Identification, inference and sensitivity analysis for causal mediation effects. Stat Sci. 1, 51-71 (2010).", "Imai K., Keele L. & Tingley D. A general approach to causal mediation analysis. Psychol Methods. 15, 309-334 (2010)." )
  }
  else if(method_name == "Sobel Test"){
    reference_lists = c(reference_lists, "Sobel M. E. Asymptotic confidence intervals for indirect effects in structural equation models. Sociological Methodology. 
                        13, 290-312 (1982).")
  }
  else if (method_name == "DACT"){
    reference_lists = c(reference_lists, "Liu Z. et al. Large-scale hypothesis testing for causal mediation effects with applications in genome-wide epigenetic studies. 
                        J Am Stat Assoc. 117, 67-81 (2022).")
  }
  else if (method_name == "Preacher and Hayes") {
    reference_lists = c(reference_lists,  "Preacher K. J. & Hayes A. F. SPSS and SAS procedures for estimating indirect effects in simple mediation models. Behav Res 
                        Methods. 36, 717-731 (2004).",
                        "Preacher K. J. & Hayes A. F. Asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models.
                        Behav Res Methods. 40, 879-91 (2008).")
  }
  else if (method_name == "MedTest") {
    reference_lists = c(reference_lists, "Zhang J., Wei Z. & Chen J. A distance-based approach for testing the mediation effect of the human microbiome. Bioinformatics.
                        34, 1875-1883 (2018). ")
  }
  
  ## FDR
  if(FDR == "Yes"){
    reference_lists = c(reference_lists, "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300.")
  }
  
  print(length(reference_lists))
  
  if(length(reference_lists) == 0) {
    return(NULL)
  } 
  else if (length(reference_lists) == 1) {

    return(reference_lists)
  } 
  else {
    for (i in seq(1:length(reference_lists))){
      reference_lists[i] = paste(i, ". ", reference_lists[i], sep = "")
    }
    return(reference_lists)
  }
}