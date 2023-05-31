############ 
# Packages # 
############
library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(fossil)
library(picante)
library(entropart)
library(ggplot2)
library(stringr)
#library(devtools)
#library(reticulate)
#library(gridExtra)

#############
# Taxonomic #
#############

select.outcome_bin_con <- function(sam_dat, treat.var, sel.treat.var){
  
  ind_bin = which((treat.var != sel.treat.var) & sapply(treat.var, function(x){length(table(sam_dat[,x ]))}) == 2)
  ind_con = which(treat.var != sel.treat.var & sapply(treat.var, function(x){length(table(sam_dat[,x ]))}) != 2 & sapply(treat.var, function(x){length(na.omit(as.numeric(sam_dat[[x]])))}) !=0)
  
  return(treat.var[c(ind_bin, ind_con)])
}

add_NA <- function(taxa.out, tax.tab) {
  taxa.out.ori <- taxa.out
  tax.tab <- as.data.frame(tax.tab)
  for(type in 1:length(taxa.out)) {
    na.num <- 1
    for (rank in 1:6) {
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))) {
        if(substring(colnames(taxa.out[[type]][[rank]])[i], nchar(colnames(taxa.out[[type]][[rank]])[i])-2) == "NA ") {
          colnames(taxa.out[[type]][[rank]])[i] <- paste0(colnames(taxa.out[[type]][[rank]])[i], na.num)
          na.num <- na.num + 1
        }
      }
    }
    for (rank in 1:5) {
      print(rank)
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))){
        print(i)
        try(colnames(taxa.out[[type]][[rank+1]]) <- str_replace(colnames(taxa.out[[type]][[rank+1]]), colnames(taxa.out.ori[[type]][[rank]])[i], colnames(taxa.out[[type]][[rank]])[i]), silent = TRUE) 
      }
    }
    
    for(rank in 1:6) {
      for(i in 1:length(unique(colnames(taxa.out[[type]][[rank]])))) {
        ind <- tax.tab[[rank+1]] == unique(colnames(taxa.out.ori[[type]][[rank]]))[i]
        
        tax.tab[[rank+1]][ind] <- unique(colnames(taxa.out[[type]][[rank]]))[i]
      }
    }
  }
  return(list(taxa.out = taxa.out, tax.tab = tax.tab))
}


tax.trans <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "NANANA") {
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
  
}


taxa.names.rank <- function(taxa.out, include = TRUE){ 
  taxon.names <- list()
  taxon.names <- lapply(taxa.out, function(x) str_split(names(x), ";"))
  
  
  if (include){
    dup.list <- list(NA,NA,NA,NA,NA,NA)
    ranks <- c("K_", "P_", "C_", "O_", "F_", "G_", "S_")
  }else{
    dup.list <- list(NA,NA,NA,NA,NA)
    ranks <- c("K_", "P_", "C_", "O_", "F_", "G_")
  }
  
  
  taxon.names.rank <- list()
  for(rank in 1:(5+include)){
    print(rank)
    taxon <- lapply(taxon.names[[rank]], function(x) str_sub(x,start = 3))
    taxon.names.rank[[rank]] <- sapply(taxon, tail, 1)
    
    if(length(taxon.names.rank[[rank]]) != length(unique(taxon.names.rank[[rank]]))){
      duplicated.taxons <- unique(taxon.names.rank[[rank]][duplicated(taxon.names.rank[[rank]])])
      
      for(i in 1:length(duplicated.taxons)){
        duplicated.taxon <- duplicated.taxons[i]
        ind.dup <- which(taxon.names.rank[[rank]] %in% duplicated.taxon)
        
        for(j in 1:length(ind.dup)){
          duplicated.taxon <- paste(duplicated.taxon,"*",collapse = "")
          taxon.names.rank[[rank]][ind.dup[j]] <- duplicated.taxon 
          dup.list[[rank]][j] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[j]]), sep = ""), collapse = " | "), sep = "")
        }
      }
    }
  }
  names(taxon.names.rank) <- names(taxa.out)[1:(5+include)]
  return(list(names = taxon.names.rank, duplicates = dup.list))
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0 | x == "0.00")
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1 | x == "1.00")
  x[ind.1] <- ">.999"
  return(x)
}

#############
# Mediation #
#############
library(mediation)

cov.func_2 <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind.pri <- colnames(sam.dat) %in% sel.pri.var
  ind.mon.sin.rev <- mon.sin.rev.bin.con$is.mon | mon.sin.rev.bin.con$is.rev | mon.sin.rev.bin.con$is.sin
  return(colnames(sam.dat)[!(ind.pri | ind.mon.sin.rev)])
}

name_split_list <- function(taxa){
  list_split <- list()
  list_mat <- list() 
  
  for (i in 2:length(taxa)){
    name_split <- matrix(unlist(lapply(1:ncol(taxa[[i]]), function(x) strsplit(colnames(taxa[[i]])[[x]], ";"))), ncol = i+1, byrow = TRUE)
    vec <- names(which(table(name_split[,i]) == 1))
    
    pos <- which(name_split[,i] %in% vec)
    list_split[[i]] <- pos 
    list_mat[[i]] <- name_split
  }
  list_split[[1]] <- NA
  return(list_split)
}

med_result_concate <- function(taxa, result_med_total, uniq_list, inc){
  
  new_list <- result_med_total
  
  for (i in 2:(5+inc)){
    print(i)
    taxa.med <- result_med_total[[i]]
    taxon <- taxa[[i]]
    taxon_last <- taxa[[i-1]]
    uniq_list_ind <- uniq_list[[i]] #here
    colname_here <- colnames(taxa[[i]])
    
    for (k in 1:length(taxa.med)){
      print(k)
      taxon.med <- taxa.med[[k]]
      taxa.uniq.name <- colname_here[uniq_list_ind]
      
      #strsplit(taxa.uniq.name[j], split = ";")[[1]][-length(strsplit(taxa.uniq.name[j], split = ";")[[1]])]
      
      start_num <- 1 
      for (j in uniq_list_ind){
        
        uniq.list <- paste0(taxa.uniq.name[start_num], collapse = "|") ==  colnames(taxon)
        uniq.list.last <- paste0(strsplit(taxa.uniq.name[start_num], split = ";")[[1]][-length(strsplit(taxa.uniq.name[start_num], split = ";")[[1]])], collapse = ";") == colnames(taxon_last)
        taxon.med[j, ] <- new_list[[i-1]][[k]][which(uniq.list.last),]
        start_num = start_num + 1 
      }
      new_list[[i]][[k]] <- taxon.med
    }
  }
  return(new_list)
}

med_result_concate_dact <- function(taxa, result_med_total, uniq_list, inc){
  
  list_1 <- result_med_total[[1]]
  list_2 <- result_med_total[[2]]
    
  new_list_1 <- result_med_total[[1]]
  new_list_2 <- result_med_total[[2]]
  
  for (i in 2:(5+inc)){
    
    taxa.med_1 <- list_1[[i]]
    taxa.med_2 <- list_2[[i]]
    
    #taxa.med <- result_med_total[[i]]
    taxon <- taxa[[i]]
    taxon_last <- taxa[[i-1]]
    uniq_list_ind <- uniq_list[[i]] #here
    colname_here <- colnames(taxa[[i]])
    
    taxon.med_1 <- taxa.med_1 
    taxon.med_2 <- taxa.med_2 
    
    taxa.uniq.name <- colname_here[uniq_list_ind]
    
    #strsplit(taxa.uniq.name[j], split = ";")[[1]][-length(strsplit(taxa.uniq.name[j], split = ";")[[1]])]
    
    start_num <- 1 
    for (j in uniq_list_ind){
      
      uniq.list <- paste0(taxa.uniq.name[start_num], collapse = "|") ==  colnames(taxon)
      uniq.list.last <- paste0(strsplit(taxa.uniq.name[start_num], split = ";")[[1]][-length(strsplit(taxa.uniq.name[start_num], split = ";")[[1]])], collapse = ";") == colnames(taxon_last)
      taxon.med_1[j, ] <- new_list_1[[i-1]][which(uniq.list.last),]
      taxon.med_2[j, ] <- new_list_2[[i-1]][which(uniq.list.last),]
      
      start_num = start_num + 1 
    }
    new_list_1[[i]] <- taxon.med_1
    new_list_2[[i]] <- taxon.med_2
  }
  return(list(med = new_list_1, out = new_list_2))
}

######################
# Mediation  Package #
######################

convert_numeric_factor <- function(vec){
  if (sum(is.na(as.numeric(vec))) == sum(is.na(vec)) & length(table(vec)) > 2){
    vec <- as.numeric(vec)
  }else{
    vec <- as.factor(vec)
  }
  return(vec)
}

#response variabel 이 binary인지 continuous 인지 파악하는 함수 
binary_cont <- function(vec){
  if (length(table(vec)) == 2){
    type = "binary"
  }else if(length(table(vec)) > 2){
    type = "continuous"
  }
  return(type)
}


mediation.taxon.ind <- function(sam.dat, taxon_med, taxa_uniq_name,
                                exposure, covariates, 
                                outcome, interaction = TRUE, 
                                regression.method = "logistic", 
                                method = "bootstrap", boot.method, n.sim = 2000, inc, rank){
  
  final_list <- list() 
  
  rep_str <- c("-" = ".", ";" = ".", ":" = ".", " " = "", "\\(" = ".", "\\)" = ".", "\\]" = ".", "\\[" = ".", "\\*" = ".", "^" = ".", "&" = ".", "\\=" = ".")
  taxon.med <- taxon_med
  
  rownames(sam.dat) <- str_replace_all(rownames(sam.dat), rep_str)
  colnames(taxon.med) <- str_replace_all(colnames(taxon.med), rep_str)
  
  dat <<- cbind(taxon.med, sam.dat)  
  
  dat[, exposure] <- convert_numeric_factor(dat[, exposure])
  
  if (length(covariates) != 0){
    for (cov in covariates){
      dat[, cov] <<- convert_numeric_factor(dat[, cov])
    }
  }
  
  if (sum(is.na(dat[[outcome]])) >= 1) {
    ind <- which(is.na(dat[[outcome]]))
    dat <- dat[-ind, ]
  }
  
  mediator.list <- colnames(taxon.med)
  colnames(dat) <- c(mediator.list, colnames(sam.dat))
  
  if (treat.type == "binary"){
    acme.average <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    acme.control <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    acme.treated <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    ade.average <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    total <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
  } else if (treat.type == "continuous"){
    acme.average <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    ade.average <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
    total <- matrix(NA, nrow = length(mediator.list), ncol = 4, dimnames = list(mediator.list, c("Est", "Lower", "Upper", "P.value")))
  }
  
  med.out <- list() 
  list_binary <- list() 
  list_continuous <- list() 
  
  if(length(taxa_uniq_name) !=0){
    
    for (med in mediator.list){
      print(med)
      
      if(med %in% colnames(taxon.med)[taxa_uniq_name]){
        print("inside")
        
        # if (inc){
        #   incProgress(0.1/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus", "Species")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        # 
        # }else{
        #   if (as.numeric(rank) != 6){
        #     incProgress(0.13/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        #   }
        # }
        
        med.out[[med]] <- med.result 
        
        acme.average[med, ] <- rep(NA, 4)
        acme.control[med,] <-  rep(NA, 4)
        acme.treated[med,] <- rep(NA, 4)
        ade.average[med,] <- rep(NA, 4)
        total[med,] <- rep(NA, 4)
        
        
      }else{
        # if (inc){
        #   incProgress(0.1/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus", "Species")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        # 
        # }else{
        #   if (as.numeric(rank) != 6){
        #     incProgress(0.13/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        #   }
        # }
        
        if (length(covariates) != 0){
          f1 <<- as.formula(paste(med, "~", exposure, "+", paste(covariates, collapse = "+"), sep = " "))
          med.fit <- lm(f1, data = dat)
          
          if (interaction == TRUE){
            f2 <<- as.formula(paste(outcome, "~", med, "*", exposure, "+", paste(covariates, collapse = "+"), sep = " "))
          }else{
            f2 <<- as.formula(paste(outcome, "~", med, "+", exposure, "+", paste(covariates, collapse=" + "), sep= " "))
          }
        }else{
          f1 <<- as.formula(paste(med, "~", exposure, sep = " "))
          med.fit <- lm(f1, data = dat)
          
          if (interaction == TRUE){
            f2 <<- as.formula(paste(outcome, "~", med, "*", exposure, sep = " "))
          }else{
            f2 <<- as.formula(paste(outcome, "~", med, "+", exposure, sep= " "))
          }
        }
        
        if (outcome.type == "binary") {
          
          print(f2)
          if (regression.method == "probit") {
            out.fit <- glm(f2, data = dat, family = binomial("probit"))
          }
          else if (regression.method == "logistic") {
            out.fit <- glm(f2, data = dat, family = binomial("logit"))
          } 
          
        }else if (regression.method == "linear") {
          out.fit <- lm(f2, data = dat)
        }
        
        if (method == "quasi-Bayesian") {
          set.seed(0705) 
          med.result <- summary(mediate(med.fit, out.fit, treat = exposure, mediator = med,
                                        robustSE = TRUE, sims = n.sim))
        }else if (method == "bootstrap") {
          set.seed(0705) 
          med.result <- summary(mediate(med.fit, out.fit, treat = exposure, mediator = med,
                                        boot = TRUE, boot.ci.type = boot.method, sims = n.sim))
        }
        
        if(treat.type == "binary"){
          med.out[[med]] <- med.result 
          
          acme.average[med, ] <- c(med.result$d.avg, med.result$d.avg.ci[[1]], med.result$d.avg.ci[[2]], med.result$d.avg.p)
          acme.control[med,] <-  c(med.result$d0, med.result$d0.ci[[1]], med.result$d0.ci[[2]], med.result$d0.p)
          acme.treated[med,] <- c(med.result$d1, med.result$d1.ci[[1]], med.result$d1.ci[[2]], med.result$d1.p)
          ade.average[med,] <- c(med.result$z.avg, med.result$z.avg.ci[[1]], med.result$z.avg.ci[[2]], med.result$z.avg.p)
          total[med,] <- c(med.result$tau.coef, med.result$tau.ci[[1]], med.result$tau.ci[[2]], med.result$tau.p)
          
        } else if(treat.type == "continuous"){
          med.out[[med]] <- med.result 
          acme.average[med,] <- c(med.result$d.avg, med.result$d.avg.ci[[1]], med.result$d.avg.ci[[2]], med.result$d.avg.p)
          ade.average[med,] <- c(med.result$z.avg, med.result$z.avg.ci[[1]], med.result$z.avg.ci[[2]], med.result$z.avg.p)
          total[med,] <- c(med.result$tau.coef, med.result$tau.ci[[1]], med.result$tau.ci[[2]], med.result$tau.p)
        }
        
      }
      
      
    }
  }else{
    
    for (med in mediator.list){
      if (inc){
        incProgress(0.1/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus", "Species")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        
      }else{
        if (as.numeric(rank) != 6){
          incProgress(0.13/length(mediator.list), message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus")[rank], "_", lapply(strsplit(med, split = "_"), tail, n = 1)))
        }
      }
      
      if (length(covariates) != 0){
        f1 <<- as.formula(paste(med, "~", exposure, "+", paste(covariates, collapse = "+"), sep = " "))
        med.fit <- lm(f1, data = dat)
        
        if (interaction == TRUE){
          f2 <<- as.formula(paste(outcome, "~", med, "*", exposure, "+", paste(covariates, collapse = "+"), sep = " "))
        }else{
          f2 <<- as.formula(paste(outcome, "~", med, "+", exposure, "+", paste(covariates, collapse=" + "), sep= " "))
        }
      }else{
        f1 <<- as.formula(paste(med, "~", exposure, sep = " "))
        med.fit <- lm(f1, data = dat)
        
        if (interaction == TRUE){
          f2 <<- as.formula(paste(outcome, "~", med, "*", exposure, sep = " "))
        }else{
          f2 <<- as.formula(paste(outcome, "~", med, "+", exposure, sep= " "))
        }
      }
      
      if (outcome.type == "binary") {
        
        print(f2)
        if (regression.method == "probit") {
          out.fit <- glm(f2, data = dat, family = binomial("probit"))
        }
        else if (regression.method == "logistic") {
          out.fit <- glm(f2, data = dat, family = binomial("logit"))
        } 
        
      }else if (regression.method == "linear") {
        out.fit <- lm(f2, data = dat)
      }
      
      if (method == "quasi-Bayesian") {
        set.seed(0705) 
        med.result <- summary(mediate(med.fit, out.fit, treat = exposure, mediator = med,
                                      robustSE = TRUE, sims = n.sim))
      }else if (method == "bootstrap") {
        set.seed(0705) 
        med.result <- summary(mediate(med.fit, out.fit, treat = exposure, mediator = med,
                                      boot = TRUE, boot.ci.type = boot.method, sims = n.sim))
      }
      
      if(treat.type == "binary"){
        med.out[[med]] <- med.result 
        
        acme.average[med, ] <- c(med.result$d.avg, med.result$d.avg.ci[[1]], med.result$d.avg.ci[[2]], med.result$d.avg.p)
        acme.control[med,] <-  c(med.result$d0, med.result$d0.ci[[1]], med.result$d0.ci[[2]], med.result$d0.p)
        acme.treated[med,] <- c(med.result$d1, med.result$d1.ci[[1]], med.result$d1.ci[[2]], med.result$d1.p)
        ade.average[med,] <- c(med.result$z.avg, med.result$z.avg.ci[[1]], med.result$z.avg.ci[[2]], med.result$z.avg.p)
        total[med,] <- c(med.result$tau.coef, med.result$tau.ci[[1]], med.result$tau.ci[[2]], med.result$tau.p)
        
      } else if(treat.type == "continuous"){
        med.out[[med]] <- med.result 
        acme.average[med,] <- c(med.result$d.avg, med.result$d.avg.ci[[1]], med.result$d.avg.ci[[2]], med.result$d.avg.p)
        ade.average[med,] <- c(med.result$z.avg, med.result$z.avg.ci[[1]], med.result$z.avg.ci[[2]], med.result$z.avg.p)
        total[med,] <- c(med.result$tau.coef, med.result$tau.ci[[1]], med.result$tau.ci[[2]], med.result$tau.p)
      }
    }
  }
  
  if(treat.type == "binary"){
    list.binary <- list(acme_average = data.frame(acme.average), acme_control = data.frame(acme.control), acme_treated = data.frame(acme.treated), ade_average = data.frame(ade.average), total_effect = data.frame(total))
    
    form_regression <- f2   
    
    final_list$med_out <- med.out 
    final_list$list_binary <- list.binary 
    
    for (i in 1:length(final_list$list_binary)){
      rownames(final_list$list_binary[[i]]) <- colnames(taxon_med)
    }}else if(treat.type == "continuous"){
      list.continuous <- list(acme_average = data.frame(acme.average), ade_average = data.frame(ade.average), total_effect = data.frame(total))
      
      final_list$med_out <- med.out 
      final_list$list_continuous <- list.continuous 
      
      for (i in 1:length(final_list$list_continuous)){
        rownames(final_list$list_continuous[[i]]) <- colnames(taxon_med)
      }
      
    }
  
  print(final_list) 
  invisible(final_list)
  
 
}


# sam.dat <- dat_1 
# taxa <- dat_2 
# taxa_unique <- dat_3 
# exposure <- "ecig_status"
# covariates <- NULL
# outcome <- "gingival_inflammation"
# interac <- TRUE 
# regression.method = "logistic"
# method = "bootstrap"
# boot.method.po <- "perc"
# n.sim = 10
# inc = FALSE 

mediation.taxon.total <- function(sam.dat, taxa, taxa_unique, 
                                exposure, covariates, 
                                outcome, interac = TRUE, 
                                regression.method = "logistic", 
                                method = "bootstrap", boot.method.po, n.sim = 2000, inc){
   
  outcome.type <<- binary_cont(sam.dat[,outcome])
  treat.type <<- binary_cont(sam.dat[,exposure])
  
  list_total <- list()
  lets_see <- list()
  
  for (i in 1:(5+inc)){
    
    taxon.med <- taxa[[i]]
    taxa.uniq.pos <- taxa_unique[[i]]
  
    if (is.na(taxa.uniq.pos)){
      if(treat.type == "binary"){
        
        list_total[[i]] <- mediation.taxon.ind(sam.dat, taxon.med, NULL, exposure, covariates, 
                                               outcome, interaction = interac, regression.method, 
                                               method = method, boot.method = boot.method.po, n.sim = n.sim, inc = inc, rank = i)$list_binary
        
        
      }else if(treat.type == "continuous"){
        
        list_total[[i]] <- mediation.taxon.ind(sam.dat, taxon.med, NULL, exposure, covariates, 
                                               outcome, interaction = interac, regression.method, 
                                               method = method, boot.method = boot.method.po,  n.sim = n.sim, inc = inc, rank = i)$list_continuous
      }
    }else{
      if(treat.type == "binary"){
        
        list_total[[i]] <- mediation.taxon.ind(sam.dat, taxon.med, taxa.uniq.pos, exposure, covariates, 
                                               outcome, interaction = interac, regression.method, 
                                               method = method, boot.method = boot.method.po, n.sim = n.sim, inc = inc, rank = i)$list_binary
        
        
      }else if(treat.type == "continuous"){
        
        list_total[[i]] <- mediation.taxon.ind(sam.dat, taxon.med, taxa.uniq.pos, exposure, covariates, 
                                               outcome, interaction = interac, regression.method, 
                                               method = method, boot.method = boot.method.po,  n.sim = n.sim, inc = inc, rank = i)$list_continuous
      }
    }
  }
  names(list_total) <- names(taxa)[1:5+inc]
  invisible(list_total)
}

######################
# P-value Adjustment #
######################

q_convert_tax_med <- function(med_result){
  
              new_list <- list() 
              
              for (i in 1:length(med_result)){
                new_list[[i]] <- list() 
                
                 for (k in 1:length(med_result[[i]])){
                  
                  Q.value <- p.adjust(unlist(med_result[[i]][[k]][[4]]), method = "BH")
                  
                  new_list[[i]][[k]] <- data.frame(matrix(c(as.numeric(unlist(med_result[[i]][[k]][,1])), as.numeric(unlist(med_result[[i]][[k]][,2])), as.numeric(unlist(med_result[[i]][[k]][,3])), as.numeric(unlist(med_result[[i]][[k]][,4])), as.numeric(unlist(Q.value))), nrow=length(unlist(med_result[[i]][[k]][,2]))))
                  
                  rownames(new_list[[i]][[k]]) <- rownames(med_result[[i]][[k]])
                  colnames(new_list[[i]][[k]]) <- c("Est", "Lower", "Upper", "P.val", "Q.val")
                  
                 }
                names(new_list[[i]]) <- names(med_result[[i]])
              }
              names(new_list) <- names(med_result)
              
              return(new_list)
            }
                    
taxa.med.sep <- function(med_result, result_type = "total_effect", inc){
  
  result_1 <<- list() 
  
  for(i in 1:(5+inc)){
    result_1[[i]] <- med_result[[i]][[result_type]]
  }
  names(result_1) <- names(med_result)
  
  return(result_1)
}

#######
# LDM #
#######

med.tax.ldm <- function(sam.dat, taxa, exposure, 
                            covariates, outcome, inc){
  final_list <- list()
  
  assign(exposure, convert_numeric_factor(unlist(sam.dat[, exposure])))  #여기가 문제인듯 
  assign(outcome, convert_numeric_factor(unlist(sam.dat[, outcome])))
  
  if(length(covariates) != 0){
    for (cov in covariates){
      assign(cov, convert_numeric_factor(unlist(sam.dat[,cov])))
      #sam.dat[,cov] <- convert_numeric_factor(data.frame(sam.dat)[,cov])
    }
  }
  
  global_p <- c()
  sig_med <- list() 
  global_p_bivariate <- list() 
  biv_cov1 <- list()
  biv_cov2 <- list() 

  for (i in 1:length(taxa)){
    taxa.ind <<- taxa[[i]]
    
    # if (inc){
    #   incProgress(0.1, message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus", "Species")[i]))
    # }else{
    #   if (i != 6){
    #     incProgress(0.13, message = paste0("Calculating: ", c("Phylum", "Class", "Order", "Family", "Genus")[i]))
    #   }
    # }
    
    if (i == 6 & inc == FALSE){
      incProgress(0.1, message = "Displaying Results in Progress")
    }
    
    if (length(covariates) != 0){
      f1 <<- as.formula(paste("taxa.ind", "|", "(", paste(covariates, collapse = "+"), ")", "~", exposure, "+", outcome, sep = " "))
      data_mani <- data.frame(exposure = sam.dat[, exposure], outcome = sam.dat[, outcome])
      for (cov in covariates){
        data_mani <- cbind(data_mani, cov = sam.dat[, cov])
      }
    }else{
      f1 <<- as.formula(paste("taxa.ind", "~", exposure, "+", outcome, sep = " "))
      data_mani <- data.frame(exposure = sam.dat[, exposure], outcome = sam.dat[,outcome])
    }
    
    fit <- ldm(formula = f1, data = data_mani,seed = 123, n.cores = 4, test.mediation = TRUE)

    global_p <- c(global_p, fit$med.p.global.freq)
    
    sig_med[[i]] <- fit$med.detected.otu.freq
    global_p_bivariate[[i]] <- data.frame(Est = fit$F.global.freq, P.val = fit$p.global.freq) 
    
    matrix_1 <- matrix(NA, nrow = length(taxa.ind), ncol = 5, dimnames = list(names(taxa[[i]]), c("Global_F", "Global_P.val", "F", "P.val", "Q.val")))
    matrix_2 <- matrix(NA, nrow = length(taxa.ind), ncol = 5, dimnames = list(names(taxa[[i]]), c("Global_F", "Global_P.val", "F", "P.val", "Q.val")))
    
    matrix_1[,1] <- fit$F.global.freq[[1]]
    matrix_1[,2] <- fit$p.global.freq[[1]]
    matrix_1[,3] <- fit$F.otu.freq[1,]
    matrix_1[,4] <- fit$p.otu.freq[1,]
    matrix_1[,5] <- p.adjust(fit$p.otu.freq[1,], "BH") 
    
    matrix_2[,1] <- fit$F.global.freq[[2]]
    matrix_2[,2] <- fit$p.global.freq[[2]]
    matrix_2[,3] <- fit$F.otu.freq[2,]
    matrix_2[,4] <- fit$p.otu.freq[2,]
    matrix_2[,5] <- p.adjust(fit$p.otu.freq[2,], "BH") 
    
    biv_cov1[[i]] <- matrix_1
    biv_cov2[[i]] <- matrix_2 
  } 
  names(sig_med) <- names(taxa)
  names(biv_cov1) <- names(taxa)
  names(biv_cov2) <- names(taxa)
  names(global_p_bivariate) <- names(taxa)
  
  final_list$med_role <- global_p  #global p value for each taxa level 
  final_list$sig_med_otu<- sig_med #significant taxa 
  
  final_list$global_bivariate <- global_p_bivariate   # association between treatment and microbiome (cov1) / association between microbiome and outcome (conditional on treatment (cov2))
  final_list$mic_treat <- biv_cov1  #association between treatment and microbiome 
  final_list$mic_outcome <- biv_cov2  #association between microbiome and outcome conditional on treatment 
  
  invisible(final_list)
}

##############
# Sobel Test #
##############

tax.sobel <- function(sam.dat, taxa, exposure, outcome, inc){
  
  final_list <- list() 
  
  for (i in 1:6){
    taxa.ind <- taxa[[i]]
    mat <- matrix(NA, nrow = length(taxa.ind), ncol = 3, dimnames = list(names(taxa.ind), c("Z", "P.val", "Q.val")))
    
    for (j in 1:length(taxa.ind)){
      iv <- convert_numeric_factor (sam.dat[[exposure]])
      dv <- sam.dat[[outcome]] 
      mv <- taxa.ind[,j]
      
      if (i == 6 & inc == FALSE){
        incProgress(0.1, message = "Displaying Results in Progress")
      }

      fit <- mediation.test(mv, iv, dv)
      
      mat[j, 1] <- fit$Sobel[1]
      mat[j, 2] <- fit$Sobel[2]
    }
    
    mat[,3] <- NA
    
    rownames(mat) <- colnames(taxa.ind)
    final_list[[i]] <- mat
  }
  return(final_list)
}

result_med_out <- function(sam.dat, taxa, uniq_list, treat, out, model = "med", interac){
  final_list <- list() 
  
  iv <- sam.dat[[treat]]
  dv <- sam.dat[[out]]
  
  if(model == "med"){
    for (i in 1:6){
      mat <- matrix(NA, nrow = length(taxa[[i]]), ncol = 6, dimnames = list(names(taxa[[i]]), c("Effect", "SE", "P.val", "Q.val", "Lower", "Upper")))
      for (k in 1:length(taxa[[i]])){
        mv <- taxa[[i]][[k]]
        med.fit <- lm(mv ~ iv)
        result.med <- summary(med.fit)
        est <- result.med$coefficients[2, 1]
        se <- result.med$coefficients[2, 2]
        pval <- result.med$coefficients[2, 4]
        conf.int <- confint(med.fit)[2,]
        mat[k, ] <- c(est, se, pval, NA, conf.int) 
      }
      mat[, "Q.val"] <- NA
      final_list[[i]] <- mat 
    }
  }else{
    for (i in 1:6){
      mat <- matrix(NA, nrow = length(taxa[[i]]), ncol = 6, dimnames = list(names(taxa[[i]]), c("Effect", "SE", "P.val", "Q.val", "Lower", "Upper")))
      for (k in 1:length(taxa[[i]])){
        mv <- taxa[[i]][[k]]
        dv <- as.numeric(dv)
        
        if (interac){
          out.fit <- lm(dv ~ iv * mv) 
        }else{
          out.fit <- lm(dv ~ iv + mv) 
        }
        
        result.out <- summary(out.fit)
        est <- result.out$coefficients[3, 1]
        se <- result.out$coefficients[3, 2]
        pval <- result.out$coefficients[3, 4]
        conf.int <- confint(out.fit)[3,]
        mat[k, ] <- c(est, se, pval, NA, conf.int) 
      }
      mat[, "Q.val"] <- NA 
      final_list[[i]] <- mat
    }
  }
  
  names(final_list) <- names(taxa)
  
  return(final_list)
}
# 
# result <- result_sobel_med_ori
# inc <- FALSE 
med_result_concate_sobel <- function(taxa, result, uniq_list, inc){
  new_list <- result 
  for (i in 2:(5+inc)){
    taxa.med <- result[[i]]
    
    taxon <- taxa[[i]]
    taxon_last <- taxa[[i-1]]
    
    uniq_list_ind <- uniq_list[[i]]
    colname_here <- colnames(taxa[[i]])
    
    taxon.med <- taxa.med
    taxa.uniq.name <- colname_here[uniq_list_ind]
    
    start_num <- 1
    for (j in uniq_list_ind){
      uniq.list <- paste0(taxa.uniq.name[start_num], collapse = "|") ==  colnames(taxon)
      uniq.list.last <- paste0(strsplit(taxa.uniq.name[start_num], split = ";")[[1]][-length(strsplit(taxa.uniq.name[start_num], split = ";")[[1]])], collapse = ";") == colnames(taxon_last)
      taxon.med[j,] <- new_list[[i-1]][which(uniq.list.last),]
      
      start_num = start_num + 1 
    }
    new_list[[i]] <- taxon.med
  }
  return(new_list)
}

q_convert_tax_sobel <- function(result, inc){
  for (i in 1:(5+inc)){
    result[[i]][, "Q.val"] <- p.adjust(result[[i]][, "P.val"], "BH")
  }
  return(result)
}

#########
# DACT #
########

# sam.dat <- dat_1 
# taxa <- dat_2 
# taxa_unique <- dat_3 
# exposure <- dat_4 
# outcome <- dat_5 
# covariates <- dat_6
# reg <- "logistic"
# interac = FALSE 

result_med_out_dact <- function (sam.dat, taxa, taxa_unique, exposure, outcome, covariates, reg = "linear", interac) {
  
  final_list_med <- list() 
  final_list_out <- list() 
  
  for (i in 1:6){
    
    mat_med <- matrix(NA, nrow = length(taxa[[i]]), ncol = 6, dimnames = list(names(taxa[[i]]), c("Effect", "SE", "P.val", "Q.val", "Lower", "Upper")))
    mat_out <- matrix(NA, nrow = length(taxa[[i]]), ncol = 6, dimnames = list(names(taxa[[i]]), c("Effect", "SE", "P.val", "Q.val", "Lower", "Upper")))
    
    for (k in 1:length(taxa[[i]])){
      print(k)
      mv <- taxa[[i]][[k]]
      data_mani <- data.frame(exposure = sam.dat[, exposure], outcome = sam.dat[, outcome], med = mv)
      
      if (length(covariates) != 0){
        for (cov in covariates){
          data_mani <- data.frame(cbind(data_mani, cov = sam.dat[, cov]))
        }
        
        f1 <<- as.formula(paste("med ~", exposure, "+", paste(covariates, collapse = "+"), sep = " "))
        
        if (interac == TRUE){
          f2 <<- as.formula(paste(outcome, "~ med", "*", exposure, "+", paste(covariates, collapse=" + "), sep= " "))
        }else{
          f2 <<- as.formula(paste(outcome, "~ med", "+", exposure, "+", paste(covariates, collapse=" + "), sep= " "))
        }
        
      }else{
        f1 <<- as.formula(paste("med ~", exposure, sep = " "))
        if (interac == TRUE){
          f2 <<- as.formula(paste(outcome, "~ med", "*", exposure, sep= " "))
        }else{
          f2 <<- as.formula(paste(outcome, "~ med", "+", exposure, sep= " "))
        }
      }
      
      med.fit <- lm(formula = f1, data = data_mani)
      result.med <- summary(med.fit)
      
      if (reg == "linear"){
        out.fit <- lm(formula = f2, data = data_mani)
      }else{
        out.fit <- glm(formula = f2, data = data_mani, family = binomial(link = "logit"))
      }
      
      result.out <- summary(out.fit)
      
      est_med <- result.med$coefficients[2, 1]
      se_med <- result.med$coefficients[2, 2]
      pval_med <- result.med$coefficients[2, 4]
      conf.int_med <- confint(med.fit)[2,]
      mat_med[k, ] <- c(est_med, se_med, pval_med, NA, conf.int_med) 
      
      
      est_out <- result.out$coefficients[2, 1]
      se_out <- result.out$coefficients[2, 2]
      pval_out <- result.out$coefficients[2, 4]
      conf.int_out <- confint(out.fit)[2,]
      mat_out[k, ] <- c(est_out, se_out, pval_out, NA, conf.int_out) 
    }
    mat_med[, "Q.val"] <- NA 
    mat_out[, "Q.val"] <- NA
    final_list_med[[i]] <- mat_med
    final_list_out[[i]] <- mat_out
  }
  return(list(med = final_list_med, out = final_list_out))
}

q_convert_tax_dact <- function(result, inc){
  final_list_med <- result[[1]]
  final_list_out <- result[[2]]
  
  for (i in 1:(5+inc)){
   final_list_med[[i]][, "Q.val"] <- p.adjust(final_list_med[[i]][, "P.val"], "BH")
   final_list_out[[i]][, "Q.val"] <- p.adjust(final_list_out[[i]][, "P.val"], "BH")
  }
  
  return(list(med = final_list_med, out = final_list_out))
}

tax.dact <- function(dact_result, correc = "JC", inc){
  final_list <- list() 
  
  for (i in 1:6){
    result <- DACT(p_a = dact_result[[1]][[i]][,3], p_b = dact_result[[2]][[i]][,3], correction = correc)
    final_list[[i]] <- data.frame(P.val = result, Q.val = p.adjust(result, "BH"))
  }
  
  return(final_list)
}

tax.dact.new <- function(dact_result, correc = "JC", inc){
  final_list <- list() 
  
  for (i in 1:6){
    
    result <- DACT(p_a = dact_result[[1]][[i]][,4], p_b = dact_result[[2]][[i]][,4], correction = correc)
    final_list[[i]] <- data.frame(P.val = result, Q.val = p.adjust(result, "BH"))
  }
  
  return(final_list)
}

#####################
# Process (Model 4) #
#####################

med.tax.process <- function(sam_dat, taxa, exposure, covariates, outcome, num_sim){
  
  both_list <- list()
  final_list_direct <- list() 
  final_list_indirect <- list() 
  
  for (i in 1:length(taxa)){
    taxa_rank <- taxa[[i]]
    
    mat_direct <- matrix(NA, nrow = length(taxa_rank), ncol = 5, dimnames = list(colnames(taxa_rank), c("Effect", "Lower", "Upper", "P.val", "Q.val")))
    mat_indirect <- matrix(NA, nrow = length(taxa_rank), ncol = 3, dimnames = list(colnames(taxa_rank), c("Effect", "Lower", "Upper")))
    
    for (j in 1:ncol(taxa_rank)){
      #j = 1
      taxa.ind <- taxa_rank[[j]]
      mv <- taxa.ind 
      iv <- convert_numeric_factor(unlist(sam_dat[, exposure]))
      dv <- convert_numeric_factor(unlist(sam_dat[, outcome]))
      dat <- data.frame(cbind(dv, iv, mv))
      colnames(dat) <- c("dv", "iv", "mv")
      
      if (length(covariates) != 0){
        for (cov in covariates){
          dat <- cbind(dat, sam_dat[, cov])
        }
        
        fit <- process(data = dat, y = "dv", x = "iv", m = "mv", cov = covariates, model = 4, boot = num_sim, save = 2)
      }else{
        fit <- process(data = dat, y = "dv", x = "iv", m = "mv", model = 4, boot = num_sim, save = 2)
      }
      
      nrow <- nrow(fit)
      direct <- fit[nrow-1, ]
      indirect <- fit[nrow, ]
      
      mat_direct[j, c(1, 2, 3, 4)] <- direct[c(1, 5, 6, 4)]
      mat_indirect[j, c(1, 2, 3)] <- indirect[c(1, 3, 4)]
    }
    mat_direct[, 5] <- p.adjust(mat_direct[,4], "BH")
    
    final_list_direct[[i]] <- mat_direct
    final_list_indirect[[i]] <- mat_indirect
  }
  
  both_list <- list(final_list_direct, final_list_indirect)
  
  return(both_list)
}

###############
# Forest Plot #
###############

taxa.forest.plot.pages <- function(all.taxa.q.out, species.include, mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.val < 0.05)  
    } else {
      ind.sig <- which(out$P.val < 0.05)
    }
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}


forest_med_tax <- function(med_result, result_type = c("ACME(average)", "ACME(control)", "ACME(treated)", "ADE(average)", "Total Effect")){
  names <- unlist(lapply(strsplit(rownames(med_result), split = ";"), tail, n = 1))
  plot_mat <- as.matrix(rbind(c("Taxon", "Est", "   Q.value"), cbind(names, decimal_adjust(med_result[,"Est"]), p.value.0.1_char(decimal_adjust(med_result[,"Q.value"])))))
  
  plot_mat_list <- list(list(), list(), list())
  plot_mat_list[[1]] <- plot_mat[, 1]
  plot_mat_list[[2]] <- plot_mat[, 2]
  plot_mat_list[[3]] <- plot_mat[, 3]
  
  list_high <- list() 
  list_high[[1]] <- gpar(lwd = 2.5, col = "#000044")
  list_high[[2]] <- gpar(lwd = 2.5, col = "#000044") 
  names(list_high) <- c(2, nrow(med_result) + 2)
  
  Plot <- forestplot(labeltext = plot_mat_list, mean = c(NA, med_result[,"Est"]), lower = c(NA, med_result[,"Lower"]), upper = c(NA, med_result[,"Upper"]), 
                     zero = 1, hrzl_lines = TRUE, new_page = TRUE, boxsize = 0.18, line.margin = unit(0.12, "cm"), grid = 0, lineheight = "lines", colgap = unit(1 , "cm"), graphwidth = unit(8, "cm"),
                     col = fpColors(box = rgb(1, 0, 0, 0.5), line = "black", summary = "red3"), mar = unit(c(0.5,0,0.5,0), "cm"), xlab = "95 % Confidence Interval",
                     txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                    ticks=gpar(fontfamily="", cex=0.75),
                                    xlab=gpar(fontfamily="", cex=0.75)), title = result_type) 
  
  return(Plot)
}

taxa.forest.plot.pages1 <- function(all.taxa.q.out, taxa.names.out, species.include, report.type = "Est", mult.test.cor = TRUE) {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "Effect") {
    report.txt <- "Eff."
  }

  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out[,"Q.val"] < 0.05)  
    } else {
      ind.sig <- which(out[,"P.val"] < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- total - total * (total%/%num.pages) # 수정 mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab <- list()
  all.ci.tab <- list()
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      if (mult.test.cor){
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value", "Q-value"), 1, 6)  
      } else {
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value"), 1, 5)
      }
      ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
    }
    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        if (length(sig.out) == 5){
          info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
          info.data <- rbind(info.data, sig.out)
          info.ci <- rbind(info.ci, cbind(sig.out[report.type], sig.out["Lower"], sig.out["Upper"]))
          info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
        }else{
          info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
          info.data <- rbind(info.data, sig.out)
          info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
          info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
        }
      }
    }
    
    taxa_names_rank <<- cbind(info.rank, info.taxa) 
    
    info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.type)],3), nsmall = 3), p.value.0.1(info.data[, "P.val"]), p.value.0.1(info.data[, "Q.val"])))
    
    #come come
                            
                            
                            #format(p.value.0.1(info.data[, c("P.val", "Q.val")]))))
    
    for(p in 1:num.pages) {
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all)
      
      all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      colnames(info.ci) <- NULL
      all.ci.tab[[p]] <- rbind(as.matrix(ci.tab.all), as.matrix(info.ci[c(initial:(initial+actual.plot.per.page-1)),]))
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
}


taxa.chat.rank.name <- function(all.taxa.q.out, taxa.names.out, species.include, report.type = "Est", mult.test.cor = TRUE) {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "Effect") {
    report.txt <- "Eff."
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out[,"Q.val"] < 0.05)  
    } else {
      ind.sig <- which(out[,"P.val"] < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  rank <- 1
  tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
  
  info.taxa <- numeric()
  info.rank <- numeric()

    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    
    taxa_names_rank <- cbind(info.rank, info.taxa) 
   
  return(taxa_names_rank)
}


taxa.forest.plot.pages1_ind <- function(all.taxa.q.out, taxa.names.out, species.include) {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  report.txt <- "Eff."
  report.type <- "Effect"
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    ind.sig <- which(out[,"Lower"]*out[,"Upper"] > 0)  
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- total - total * (total%/%num.pages) 
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab <- list()
  all.ci.tab <- list()
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    
    if (text.tab.all == 0) {
      text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "Boot Lower", "Boot Upper"), 1, 6)  
      
      ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
    }
    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }

    info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.type, "Lower", "Upper")],3), nsmall = 3)))
    
    for(p in 1:num.pages) {
      
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all)
      
      all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab[[p]] <- rbind(ci.tab.all, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
}


taxa.forest.plot.pages2 <- function(page.taxa.q.out, page) {
  
  text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
  ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
  
  if(is.null(text.tab.all) & is.null(ci.tab.all)){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
    for(i in 1:length(page.taxa.q.out$all.text.tab)){
      str.max[[i]] <- str.max[[i]][,3]
    }
    maxStr <- max(unlist(str.max))
    if(!is.numeric(maxStr)){
      maxStr <- 0
    }
    
    text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
    
    if(nrow(ci.tab.all) <= 45){  
      forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                 hrzl_lines=TRUE, new_page=TRUE, boxsize=0.18, grid=0, colgap = unit(1, "cm"), graphwidth = unit(8, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                 col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                 txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                ticks=gpar(fontfamily="", cex=0.75),
                                xlab=gpar(fontfamily="", cex=0.75)))
    } else {
      forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                 hrzl_lines=TRUE, new_page=TRUE, boxsize=0.18, grid=0, colgap = unit(1, "cm"), graphwidth = unit(8, "cm"), 
                 col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                 txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                ticks=gpar(fontfamily="", cex=0.75),
                                xlab=gpar(fontfamily="", cex=0.75)))
    }
  }
}

taxa.forest.plot.pages1_dact <- function(dact_result, dact_med, dact_out, taxa.names.out, species.include, report.type = "Est", mult.test.cor = TRUE) {
  sig.by.rank <- list()
  
  if(species.include){
    range <- 6
  }else{
    range <- 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "Effect") {
    report.txt <- "Eff."
  }
  
  for(i in 1:5+species.include) {
    out <- dact_result[[i]]
    ind.sig <- which(out[,"Q.val"] < 0.05) 
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- total - total * (total%/%num.pages) # 수정 mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab_1 <- list()   #from here 
  all.ci.tab_1 <- list()
  
  all.text.tab_2 <- list() 
  all.ci.tab_2 <- list()    #to here 
  
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      text.tab.all_1 <- matrix(c("ID", "Rank", "Taxon", "P-value", "Q-value", paste0(report.txt, " (Med)")), 1, 6)   #from here 
      ci.tab.all_1 <- matrix(c(NA, NA, NA), 1, 3)
      
      text.tab.all_2 <- matrix(c(paste0(report.txt, " (Out)")), 1, 1)
      ci.tab.all_2 <- matrix(c(NA, NA, NA), 1, 3)                        #to here 
    }
    
    for (rank in 1:(5+include)) {
      if(!is.null(dact_result[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- dact_result[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.type)],3), nsmall = 3), format(p.value.0.1(info.data[, c("P.val", "Q.val")]))))
    
    for(p in 1:num.pages) {
      
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all_1)
      
      all.text.tab_1[[p]] <- rbind(as.matrix(text.tab.all_1), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab_1[[p]] <- rbind(ci.tab.all_1, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
      all.text.tab_2[[p]] <- rbind(as.matrix(text.tab.all_2), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab_2[[p]] <- rbind(ci.tab.all_2, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab_1 = all.text.tab_1, all.ci.tab_1 = all.ci.tab_1, all.text.tab_2 = all.text.tab_2, all.ci.tab_2 = all.ci.tab_2))
}


taxa.forest.plot.pages1_sobel <- function(sobel_result, sobel_med, sobel_out, taxa.names.out, species.include, report.type = "Est", mult.test.cor = TRUE) {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "Effect") {
    report.txt <- "Eff."
  }
  
  for(i in 1:range) {
    out <- sobel_result[[i]]
    ind.sig <- which(out[,"Q.val"] < 0.05) 
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- total - total * (total%/%num.pages) # 수정 mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab_1 <- list()   #from here 
  all.ci.tab_1 <- list()
  
  all.text.tab_2 <- list() 
  all.ci.tab_2 <- list()    #to here 
  
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.dat_med <- numeric() 
    info.dat_out <- numeric()
    info.rank <- numeric()
    info.ci_sobel <- numeric()
    info.ci_med <- numeric() 
    info.ci_out <- numeric() 
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      text.tab.all_1 <- matrix(c("ID", "Rank", "Taxon", "P-value", "Q-value", paste0(report.txt, " (Med)")), 1, 6)   #from here 
      ci.tab.all_1 <- matrix(c(NA, NA, NA), 1, 3)
      
      text.tab.all_2 <- matrix(c(paste0(report.txt, " (Out)")), 1, 1)
      ci.tab.all_2 <- matrix(c(NA, NA, NA), 1, 3)                       
    }
    
    for (rank in 1:range) {
      if(!is.null(sobel_result[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- sobel_result[[rank]][sig.by.rank[[rank]],]
        med.sig.out <- sobel_med[[rank]][sig.by.rank[[rank]],]
        out.sig.out <- sobel_out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.dat_med <- rbind(info.dat_med, med.sig.out)
        info.dat_out <- rbind(info.dat_out, out.sig.out)
        
        
        if (class(sig.out) == "numeric" | length(med.sig.out) == 3){
          info.ci_sobel <- rbind(info.ci_sobel, cbind(sig.out["P.val"], sig.out["Q.val"]))
          info.ci_med <- rbind(info.ci_med, cbind(med.sig.out["Effect"], med.sig.out["Lower"], med.sig.out["Upper"]))
          info.ci_out <- rbind(info.ci_out, cbind(out.sig.out["Effect"], out.sig.out["Lower"], out.sig.out["Upper"]))
          
        }else{
          info.ci_sobel <- rbind(info.ci_sobel, cbind(sig.out[,"P.val"], sig.out[,"Q.val"]))
          info.ci_med <- rbind(info.ci_med, cbind(med.sig.out[, "Effect"], med.sig.out[, "Lower"], med.sig.out[, "Upper"]))
          info.ci_out <- rbind(info.ci_out, cbind(out.sig.out[, "Effect"], out.sig.out[, "Lower"], out.sig.out[, "Upper"]))
          
        }
        
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    
    hey_taxa_names <<- info.taxa
    info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, 
                            p.value.0.1(info.data[, "P.val"]), p.value.0.1(info.data[, "Q.val"]),
                            format(round(info.dat_med[,"Effect"],3), nsmall = 3),
                            format(round(info.dat_out[,"Effect"],3), nsmall = 3)
                            ))
    
    
    for(p in 1:num.pages) {
      
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all_1)
      
      all.text.tab_1[[p]] <- rbind(as.matrix(text.tab.all_1), info[c(initial:(initial+actual.plot.per.page-1)), 1:6])
      all.ci.tab_1[[p]] <- rbind(ci.tab.all_1, info.ci_med[c(initial:(initial+actual.plot.per.page-1)),])
      all.text.tab_2[[p]] <- rbind(as.matrix(text.tab.all_2), as.matrix(info[c(initial:(initial+actual.plot.per.page-1)), 7]))
      all.ci.tab_2[[p]] <- rbind(ci.tab.all_2, info.ci_out[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab_1 = all.text.tab_1, all.ci.tab_1 = all.ci.tab_1, all.text.tab_2 = all.text.tab_2, all.ci.tab_2 = all.ci.tab_2))
}

#come here
chat_gpt_mediation <- function(api_key, result_names, prim_var, out_var){
  chat_list <- list() 
  
  for (name in result_names){
    Sys.setenv(OPENAI_API_KEY = as.character(api_key))
    past_question <- paste("What is known about", name, "on", prim_var, "and", out_var)
    chat_list[[name]] <- ask_chatgpt(past_question)
  }
  
  names(chat_list) <- result_names 
  return(chat_list)
}
  

  
taxa.forest.plot.pages2_sobel <- function(page.taxa.q.out, page) {

  if(length(page.taxa.q.out$all.text.tab_1) == 0){
    
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
    
  }else{
    text.tab.1 <- page.taxa.q.out$all.text.tab_1[[page]]
    ci.tab.1 <- page.taxa.q.out$all.ci.tab_1[[page]]
    text.tab.2 <- page.taxa.q.out$all.text.tab_2[[page]]
    ci.tab.2 <- page.taxa.q.out$all.ci.tab_2[[page]]
    
    plot.new()
    
    str.max <- lapply(page.taxa.q.out$all.text.tab_1, nchar)
    for(i in 1:length(page.taxa.q.out$all.text.tab_1)){
      str.max[[i]] <- str.max[[i]][,3]
    }
    maxStr <- max(unlist(str.max))
    if(!is.numeric(maxStr)){
      maxStr <- 0
    }
    
    text.tab.1[,3] <- substr(text.tab.1[,3], 1, 55)
    
    
    text_1 <- cbind(text.tab.1, est_out = rep("  ", nrow(text.tab.1)), group = rep("med", nrow(text.tab.1)))
    text_2 <- cbind(rep("   ", length(text.tab.1[,1])), rep("   ", length(text.tab.1[,1])), rep("   ", length(text.tab.1[,1])), rep("   ", length(text.tab.1[,1])), rep("   ", length(text.tab.1[,1])), med_out = rep("  ", nrow(text.tab.1)), text.tab.2,  group = rep("out", nrow(text.tab.1)))
    
    
    new_text <- matrix(NA, nrow = (nrow(text_1) + nrow(text_2)), ncol = ncol(text_1)) 
    eff_med <- c()
    for (i in 1:nrow(text_1)){
      new_text[2*i -1,] <- text_1[i, ]
      new_text[2*i, ] <- text_2[i, ]
    }
    
    new_text[1, 7] <- "Eff. (Out)"
    new_text[2, 7] <- " "
    
    new_text[1,  6] <- "Eff. (Med)"
    new_text[2, 6] <- "  "
    
    ci <- matrix(NA, nrow = (nrow(text_1) + nrow(text_2)), ncol = 3) 
    
    for (i in 1:nrow(text_2)){
      ci[2*i-1, ] <- ci.tab.1[i,]
      ci[2*i,] <- ci.tab.2[i,] 
    }
    
    new_text <- cbind(new_text, ci)
    new_text <- data.frame(new_text)
    colnames(new_text) <- c("id", "rank", "taxon", "pval", "qval", "eff_med", "eff_out", "group", "mean", "lower", "upper")
    
    box_list <- list() 
    for (i in 1:nrow(new_text)){
      if(i %% 2 ==0){
        box_list[[i]] <- gpar(fill = "blue", col = "blue")
      }else{
        box_list[[i]] <- gpar(fill = "red", col = "red")
      }
    }
    
    styles <- fpShapesGp(
      box = box_list)
    
    forestplot(new_text[, c("id", "rank", "taxon", "pval", "qval", "eff_med", "eff_out")], mean = as.numeric(new_text$mean), lower = as.numeric(new_text$lower), upper = as.numeric(new_text$upper), hzrl_lines = TRUE, is.summary  = c(TRUE, rep(FALSE, nrow(new_text)-1)), new_page=TRUE, boxsize=0.18, grid=0, colgap = unit(0.6, "cm"), graphwidth = unit(8.5, "cm"), lineheight = "lines", line.margin = unit(0.16, "cm"), col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.9), gpar(fontfamily="", cex=0.9)), ticks=gpar(fontfamily="", cex=0.9), xlab=gpar(fontfamily="", cex=0.9)), shapes_gp = styles)
  }
}

duplicate.list <- function(duplicate.taxa, taxon.inplot, duplicate.full.list){
  if(length(duplicate.taxa %in% taxon.inplot)>0) {
    duplicate.taxa <- unlist(duplicate.full.list)[duplicate.taxa %in% taxon.inplot]
    par(mar=c(0, 0.5, 0, 0.5))
    text(x=0, y=0.5, paste(duplicate.taxa, collapse = "\n"), cex = 0.75, adj = c(0, NA))
  } else {
    print("?!")
    text(x=0.5, y=0.5, "")
  }
}

taxa.forest.plot.pages <- function(all.taxa.q.out, species.include, mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out[,"Q.val"] < 0.05)  
    } else {
      ind.sig <- which(out[,"P.val"] < 0.05)
    }
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

taxa.forest.plot.pages_pro_ind <- function(all.taxa.q.out, species.include) {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    ind.sig <- which(as.numeric(out[,"Lower"])*as.numeric(out[,"Upper"]) > 0)  
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

############
# Bar Plot #
############

ldm.barplot <- function(ldm.out, out.type = "mic_outcome", rank){
  first_out <<- ldm.out 
  med.out <- ldm.out[[out.type]][[rank]]
  second_out <<- med.out 
  
  if(out.type == "mic_treat"){
    global.out <- paste0("global P: ", p.value.0.1(ldm.out$global_bivariate[[rank]][1, 2]))
  }else{
    global.out <- paste0("global P: ", p.value.0.1(ldm.out$global_bivariate[[rank]][2, 2]))
  }
  
  sig.0.05 <- med.out[,"Q.val"][(0.01 <= med.out[,"Q.val"]) & (med.out[,"Q.val"] <= 0.05)]
  sig.0.01 <- med.out[,"Q.val"][(0.001 <= med.out[,"Q.val"]) & (med.out[,"Q.val"] <= 0.01)]
  sig.0.001 <- med.out[,"Q.val"][med.out[,"Q.val"] <= 0.001]
  
  if(length(c(sig.0.05, sig.0.01, sig.0.001)) == 0){
    plot.new() 
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    
    sig_tax <- unlist(lapply(strsplit(names(sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE)), ";"), function(x) rev(x)[[1]]))
    sig_tax_2 <- names(sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE))
    
    #sig_tax <- rev(strsplit(names(sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE)), ";")[[1]])[1]
    sig_tax_q <- sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE)
    
    if(length(sig_tax) >100){
      sig_tax <- c(sig_tax[1:50], sig_tax[(length(sig_tax) - 49):length(sig_tax)])
      sig_tax_q <- c(sig_tax_q[1:50], sig_tax_q[(length(sig_tax_q) - 49):length(sig_tax_q)])
    }
    
    if (length(sig_tax) <= 50) {
      si <- (1.5-0.1*length(sig_tax)^(0.5))*length(sig_tax)^(0.8)
      #par(mar=c(21 - si, 8, 19 - si, 5))
      par(mar=c(5.1, 13 ,4.1 ,2.1))
      bar <- barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = global.out, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      #abline(v = c(0, 0.01, 0.05),  col = "gold", lty = 1, lwd = 1, cex.axis = 0.75)
      #par(mar=c(5.1, 13 ,4.1 ,2.1))
      barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(0, 0.05), axes = FALSE, add = TRUE, xlab = global.out, font.lab=4)
      text(bar, x = as.numeric(sig_tax_q + 0.001), ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 0.8, col = "dark violet")
      
    }
    
    if (length(sig_tax) > 50 && length(sig_tax) <= 100) {
      par(mfrow = c(1,2))
      ind.1 <- round(length(sig_tax)/2)
      sig_tax_q.2 <- sig_tax_q[(ind.1+1):length(sig_tax_q)]
      sig_tax.2 <- sig_tax[(ind.1+1):length(sig_tax_q)]
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      bar <- barplot(height = sig_tax_q.2, space = 0.5, col.axis = "grey31", cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.2, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = global.out, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      barplot(height = sig_tax_q.2, space = 0.5, col.axis = "grey31", cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.2, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, add = TRUE, xlab = global.out, font.lab=4)
      text(bar, x = sig_tax_q + 0.001, ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 1, col = "dark violet")
      
      sig_tax_q.1 <- sig_tax_q[1:ind.1]
      sig_tax.1 <- sig_tax[1:ind.1]
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      
      bar <- barplot(height = sig_tax_q.1, space = 0.5, col.axis = "grey31",cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.1, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = global.out, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      barplot(height = sig_tax_q.1, space = 0.5, col.axis = "grey31",cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.1, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, add = TRUE, xlab = global.out, font.lab=4)
      text(bar, x = sig_tax_q + 0.001, ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 1, col = "dark violet")
      
    }
  }
}

sobel.barplot <- function(sobel.out, rank, pval = FALSE){
  first_out <- sobel.out 
  med.out <- as.matrix(sobel.out[[rank]])
  second_out <- med.out
  
  if (pval == TRUE){
    sig.0.05 <- med.out[, "P.val"][(0.01 <= med.out[, "P.val"]) & (med.out[, "P.val"] <= 0.05)]
    sig.0.01 <- med.out[, "P.val"][(0.001 <= med.out[, "P.val"]) & (med.out[, "P.val"] <= 0.01)]
    sig.0.001 <- med.out[, "P.val"][med.out[, "P.val"] <= 0.001]
  }else{
    sig.0.05 <- med.out[,"Q.val"][(0.01 <= med.out[,"Q.val"]) & (med.out[,"Q.val"] <= 0.05)]
    sig.0.01 <- med.out[,"Q.val"][(0.001 <= med.out[,"Q.val"]) & (med.out[,"Q.val"] <= 0.01)]
    sig.0.001 <- med.out[,"Q.val"][med.out[,"Q.val"] <= 0.001]
  }
  
  if(length(c(sig.0.05, sig.0.01, sig.0.001)) == 0){
    plot.new() 
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    
    print("here?")
    sig_tax <- unlist(lapply(strsplit(names(sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE)), ";"), function(x) rev(x)[[1]]))
    sig_tax_2 <- names(sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE))
    
    sig_tax_q <- sort(c(sig.0.05, sig.0.01, sig.0.001), decreasing = TRUE)
    
    if (length(sig_tax) == 1){
      si <- (1.5-0.1*length(sig_tax)^(0.5))*length(sig_tax)^(0.8)
      par(mar=c(6.5, 13 ,5.5 ,2.1))
      
      bar <- barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = NULL, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      #abline(v = c(0, 0.01, 0.05),  col = "gold", lty = 1, lwd = 1, cex.axis = 0.75)
      #par(mar=c(5.1, 13 ,4.1 ,2.1))
      barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(0, 0.05), axes = FALSE, add = TRUE, xlab = NULL, font.lab=4)
      text(bar, x = as.numeric(sig_tax_q + 0.001), ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 0.8, col = "dark violet")
      
    }
    
    if(length(sig_tax) >100){
      sig_tax <- c(sig_tax[1:50], sig_tax[(length(sig_tax) - 49):length(sig_tax)])
      sig_tax_q <- c(sig_tax_q[1:50], sig_tax_q[(length(sig_tax_q) - 49):length(sig_tax_q)])
    }
    
    if (1 < length(sig_tax) & length(sig_tax) <= 50) {
      si <- (1.5-0.1*length(sig_tax)^(0.5))*length(sig_tax)^(0.8)
      par(mar=c(5.1, 13 ,4.1 ,2.1))
      bar <- barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = NULL, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      #abline(v = c(0, 0.01, 0.05),  col = "gold", lty = 1, lwd = 1, cex.axis = 0.75)
      #par(mar=c(5.1, 13 ,4.1 ,2.1))
      barplot(height = sig_tax_q, width = 2, col.axis = "grey31", cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = "pink", border = NA, horiz = TRUE, names = sig_tax, las = 2, xlim = c(0, 0.05), axes = FALSE, add = TRUE, xlab = NULL, font.lab=4)
      text(bar, x = as.numeric(sig_tax_q + 0.001), ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 0.8, col = "dark violet")
      
    }
    
    if (length(sig_tax) > 50 && length(sig_tax) <= 100) {
      par(mfrow = c(1,2))
      ind.1 <- round(length(sig_tax)/2)
      sig_tax_q.2 <- sig_tax_q[(ind.1+1):length(sig_tax_q)]
      sig_tax.2 <- sig_tax[(ind.1+1):length(sig_tax_q)]
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      bar <- barplot(height = sig_tax_q.2, space = 0.5, col.axis = "grey31", cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.2, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = global.out, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      barplot(height = sig_tax_q.2, space = 0.5, col.axis = "grey31", cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.2, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, add = TRUE, xlab = global.out, font.lab=4)
      text(bar, x = sig_tax_q + 0.001, ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 1, col = "dark violet")
      
      sig_tax_q.1 <- sig_tax_q[1:ind.1]
      sig_tax.1 <- sig_tax[1:ind.1]
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      
      bar <- barplot(height = sig_tax_q.1, space = 0.5, col.axis = "grey31",cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.1, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, xlab = global.out, font.lab=4)
      par(lwd=1.5)
      axis(1, col = "grey31", lty = 1, lwd = 1, cex.axis = 0.75, at = c(0, 0.01, 0.05), labels= c("0", "0.01", "0.05"), col.axis = "grey31")
      
      par(mar=c(5.1, max(4.1,max(length(sig_tax))/1.8) ,4.1 ,2.1))
      barplot(height = sig_tax_q.1, space = 0.5, col.axis = "grey31",cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = "pink", border = NA, horiz = TRUE, names = sig_tax.1, las = 2, xlim = c(-0.001, 0.0501), axes = FALSE, add = TRUE, xlab = global.out, font.lab=4)
      text(bar, x = sig_tax_q + 0.001, ifelse(sig_tax_2 %in% names(sig.0.05), "*", ifelse(sig_tax_2 %in% names(sig.0.01), "**", ifelse(sig_tax_2 %in% names(sig.0.001), "***", " "))) ,cex = 1, col = "dark violet")
      
    }
  }
}

#############
# Bean Plot #
#############

data_bean_plot_bin <- function(sam.dat, taxa, exposure, outcome, rank, q_val_result){
  taxon <- taxa[[rank]]
  q_val_result_ind <- q_val_result[[rank]]
  
  exp.cat <- names(table(sam.dat[,exposure]))
  out.cat <- names(table(sam.dat[,outcome]))
  
  ind <- list() 
  
  num <- 0 
  for (i in 1:length(exp.cat)){
    for (k in 1:length(out.cat)){
      num <- num + 1 
      ind[[num]] <- which((sam.dat[,exposure] == exp.cat[i]) & (sam.dat[,outcome] == out.cat[k]))
    }
  }
  
  sig_otu_list <- list() 
  num_2 <- 0 
  for (j in 1:length(taxon)){
    dat <- c() 
    if (j %in% which(q_val_result_ind[,"Q.val"] < 0.05)){
      num_2 <- num_2 + 1
      for (p in 1:length(ind)){
        dat <- qpcR:::cbind.na(dat, taxon[ind[[p]], j])
      }
      sig_otu_list[[num_2]] <- dat[,2:ncol(dat)]
    }
  }
  names(sig_otu_list) <- names(which(q_val_result_ind[,"Q.val"] < 0.05))
  return(sig_otu_list)
}

bean_plot_bin <- function(sam.dat, taxa, exposure, outcome, rank, q_val_result){
  
  sig.num <- length(which(q_val_result[[rank]][,"Q.val"] < 0.05)) 
  nrow <- ceiling(sig.num/4)
  
  dat.bin <- data_bean_plot_bin(sam.dat, taxa, exposure, outcome, rank, q_val_result)
  q_val_result_ind <- q_val_result[[rank]]
  
  if(nrow > 0){
    par(mfrow = c(nrow, 4),  mai = c(0.6, 0.6, 0.6, 0.6))#, mai = c(0.3, 0.3, 0.3, 0.3)
    
    xlab.v<- paste0("*q:", p.value.0.1_char(decimal_adjust(as.numeric(q_val_result_ind[,"Q.val"][q_val_result_ind[,"Q.val"] < 0.05]))), sep = "")
    
    ylab.v <- c() 
    for (r in 1:length(which(q_val_result_ind[,"Q.val"] < 0.05))){
      ylab.v[r] <- rev(strsplit(names(dat.bin)[[r]], ";")[[1]])[1]
    }
    
    for (s in 1:length(which(q_val_result_ind[,"Q.val"] < 0.05))){
      try(beanplot(data.frame(dat.bin[[s]]), side = 'both', border = 'NA', what = c(1, 1, 1, 0), col = list(c("#ffffb2", "black", "black", "green"), c("#fd8d3c", "white", "black", "blue")), names = names(table(sam.dat[,exposure])), ylab=ylab.v[s], xlab = xlab.v[s]), silent = TRUE)
    }
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom', fill=c('#ffffb2','#fd8d3c'), legend = names(table(sam.dat[,outcome])), horiz = TRUE)
    #, lwd = 5, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n'
    
  }else{
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot.new()
    text(x = 0.5, y = 0.5, paste("No significant taxa are found in ", ranks[rank], sep = ""), cex = 1.2, col = "black")
  }
  
}

####################
#### Dendrogram ####
####################


taxa.dend.dat <- function(bin.var, taxa, q.out) {
  
  names <- names(table(bin.var))
  
  ind.ref <- which(bin.var == name[1])
  ind.com <- which(bin.var == name[2])
  
  taxa.added <- list()
  for(i in 1:6) {
    if(!is.null(taxa[[i]])){
      n.tax <- length(taxa[[i]])
      med.diff <- numeric()
      for(j in 1:n.tax) {  
        taxon <- taxa[[i]][,j]
        med.ref <- median(taxon[ind.ref], na.rm = TRUE)
        med.com <- median(taxon[ind.com], na.rm = TRUE)
        med.diff[j] <- med.com - med.ref
      }
      q.out[[i]] <- cbind(Est_dend = med.diff, q.out[[i]])
    }
  }
  return(q.out)
}




taxa.sig.dend <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  names(out) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  
  ci.tab.all <- matrix(NA, 1, 1)
  
  for (i in 1:6) {
    
    ind.sig <- which(out[[i]][, "Q.val"] < 0.05)
    if (length(ind.sig) > 1) {
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(sig.out)
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.val"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.val"], 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }else if (length(ind.sig) == 1) {
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(out[[i]])[ind.sig]
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[1], 3), nsmall = 3), 
                        format(round(sig.out["P.val"], 3), nsmall = 3), 
                        format(round(sig.out["Q.val"], 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out[1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }

  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}


taxa.sig.dend.genus <- function(out, tax.tab, layout.type = "twopi", species.include = FALSE) {

  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  
  ci.tab.all <- matrix(NA, 1, 1)
  
  for (i in 1:5) {
    ind.sig <- which(out[[i]][, "Q.val"] < 0.05)
    if (length(ind.sig) > 1) {
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(sig.out)
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.val"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.val"], 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }else if (length(ind.sig) == 1) {
      
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(out[[i]])[ind.sig]
      
      text.tab <- c(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[[1]], 3), nsmall = 3), 
                    format(round(sig.out[["P.val"]], 3), nsmall = 3), 
                    format(round(sig.out[["Q.val"]], 3), nsmall = 3))
      
      
      ci.tab <- cbind(sig.out[[1]])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  
  for (i in 1:5) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
    
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])

        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
  phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
  Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
  
  Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
  pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
  neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
  total.group <- c(Class, Order, Family, Genus)
  total.group <<- total.group[!total.group %in% "NA"]
  na.taxa <<- total.group[grepl("NA", total.group)]
  non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  
  
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}



taxa.sig.dend.neutral <- function(out, tax.tab, layout.type = "twopi", species.include = FALSE) {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  
  ci.tab.all <- matrix(NA, 1, 1)
  
  for (i in 1:6) {
    ind.sig <- which(out[[i]][, "Q.val"] < 0.05)
    if (length(ind.sig) > 1) {
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(sig.out)
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.val"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.val"], 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }else if (length(ind.sig) == 1) {
      
      sig.out <- out[[i]][as.vector(ind.sig),]
      taxa <- rownames(out[[i]])[ind.sig]
      
      text.tab <- c(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[[1]], 3), nsmall = 3), 
                    format(round(sig.out[["P.val"]], 3), nsmall = 3), 
                    format(round(sig.out[["Q.val"]], 3), nsmall = 3))
        
      
      ci.tab <- cbind(sig.out[[1]])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <<- total.group[!total.group %in% "NA"]
    na.taxa <<- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}

tax.tab.cleanNA <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  # print(head(tax.tab.c))
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- "NA "
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  
  ranks <- c("p_", "c_", "o_", "f_", "g_", "s_")
  for (i in 1:6) {
    na.num <- 1
    ind <- which(tax.tab.h[,i+1] != na.code)
    ind_omit <- which(tax.tab.h[,i+1] != na.code & tax.tab.h[,i] == na.code)
    tax.tab.h[ind ,i+1] <- paste(tax.tab.h[ind,i],paste(ranks[i],tax.tab.h[ind,i+1], sep = ""), sep = ";")
    
    
    if(length(ind_omit)!=0) tax.tab.h[ind_omit,c((i+1):7)] = na.code
  }
  
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

tax.tab.cleanNA <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  # print(head(tax.tab.c))
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- "NA "
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  
  ranks <- c("p_", "c_", "o_", "f_", "g_", "s_")
  for (i in 1:6) {
    na.num <- 1
    ind <- which(tax.tab.h[,i+1] != na.code)
    ind_omit <- which(tax.tab.h[,i+1] != na.code & tax.tab.h[,i] == na.code)
    tax.tab.h[ind ,i+1] <- paste(tax.tab.h[ind,i],paste(ranks[i],tax.tab.h[ind,i+1], sep = ""), sep = ";")
    
    if(length(ind_omit)!=0) tax.tab.h[ind_omit,c((i+1):7)] = na.code
  }
  
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

taxa.sig.dend.DACT <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value"), 1, 4)
  
  ci.tab.all <- matrix(NA, 1, 1)
  for (i in 1:6) {
    ind.sig <- which(out[[i]] < 0.05)
    if (length(ind.sig) > 1) {
      sig.out <- as.matrix(out[[i]])[as.vector(ind.sig),]
      taxa <- names(sig.out)
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out, 3), nsmall = 3), 
                        format(round(sig.out, 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out)
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }else if (length(ind.sig) == 1) {
      sig.out <- as.matrix(out[[i]])[as.vector(ind.sig),]
      taxa <- names(out[[i]])[ind.sig]
      
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[1], 3), nsmall = 3), 
                        format(round(sig.out, 3), nsmall = 3))
      
      ci.tab <- cbind(sig.out[1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = mediumpurple1, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = mediumpurple1, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}

tax.tab.cleanNA <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  # print(head(tax.tab.c))
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- "NA "
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  
  ranks <- c("p_", "c_", "o_", "f_", "g_", "s_")
  for (i in 1:6) {
    na.num <- 1
    ind <- which(tax.tab.h[,i+1] != na.code)
    ind_omit <- which(tax.tab.h[,i+1] != na.code & tax.tab.h[,i] == na.code)
    tax.tab.h[ind ,i+1] <- paste(tax.tab.h[ind,i],paste(ranks[i],tax.tab.h[ind,i+1], sep = ""), sep = ";")
    
    if(length(ind_omit)!=0) tax.tab.h[ind_omit,c((i+1):7)] = na.code
  }
  
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

biom.cleanSNA <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, surv.time, censor, follow, subgroup.var = NULL, subgroup = NULL) {
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  if (!is.null(subgroup.var)){
    ind <- sam.dat [[subgroup.var]] == subgroup
    sam.dat  <- sam.dat [ind, ]
  }
  
  otu.tab <- otu.tab[,colnames(otu.tab) %in% rownames(sam.dat)]
  
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
  ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  if (tax.tab.c) {
    tax.tab <- tax.tab.cleanNA(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}



biom.cleanSNA.no.both <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, surv.time, censor, follow, subgroup.var = NULL, subgroup = NULL) {
  
  otu.tab <- otu_table(biom)
  sam.dat <- sample_data(biom)
  
  if (!is.null(subgroup.var)){
    ind <- sam.dat [[subgroup.var]] == subgroup
    sam.dat  <- sam.dat [ind, ]
  }
  
  otu.tab <- otu.tab[,colnames(otu.tab) %in% rownames(sam.dat)]
  
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
  
  ind.com.otu <- rownames(otu.tab)
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  
  biom <- merge_phyloseq(otu.tab, sam.dat)
  
  biom <- otu.tab.clean.no.both(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

  
biom.cleanSNA.no.tree <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, surv.time, censor, follow, subgroup.var = NULL, subgroup = NULL) {
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  sam.dat <- sample_data(biom)
  
  if (!is.null(subgroup.var)){
    ind <- sam.dat [[subgroup.var]] == subgroup
    sam.dat  <- sam.dat [ind, ]
  }
  
  otu.tab <- otu.tab[,colnames(otu.tab) %in% rownames(sam.dat)]
  
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
  
  if (tax.tab.c) {
    tax.tab <- tax.tab.cleanNA(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, sam.dat)
  
  biom <- otu.tab.clean.no.tree(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

biom.cleanSNA.no.tax.tab <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05, surv.time, censor, follow, subgroup.var = NULL, subgroup = NULL) {
  
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  if (!is.null(subgroup.var)){
    ind <- sam.dat [[subgroup.var]] == subgroup
    sam.dat  <- sam.dat [ind, ]
  }
  
  otu.tab <- otu.tab[,colnames(otu.tab) %in% rownames(sam.dat)]
  
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

otu.tab.clean <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
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
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}


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



tax.trans.na <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "") {
  
  # tax.tab <- matrix(nrow = nrow(tax.tab.ori), ncol = ncol(tax.tab.ori))
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    
    na.num <- 1
    for(i in 1:length(tax)) {
      na.taxa <- substring(tax[i], nchar(tax[i])-2)
      if(na.taxa == "_NA") {
        new.tax <- paste0(tax[i], na.num)
        ind.na.tab <- which(tax.tab[,j+1] == tax[i])
        tax.tab[ind.na.tab, j+1] <<- new.tax
        rare.tax.tab[ind.na.tab, j+1] <<- new.tax
        tax[i] <- new.tax
        na.num <- na.num + 1
      }
      
    }
    
    
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
}
