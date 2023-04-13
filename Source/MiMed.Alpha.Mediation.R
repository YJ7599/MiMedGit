library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(quantreg)
library(fossil)
library(picante)
library(entropart)
library(mediation)
library(forestplot)
library(dplyr)
library(bda)
library(DACT)
library(dashboardthemes)

########################################
#  Quality Control and Transformation  #
########################################

rem.tax.d <- c("", "metagenome", "gut metagenome", "mouse gut metagenome")
rem.tax.str.d <- c("uncultured", "incertae", "Incertae", "unidentified", "unclassified", "unknown")

tax.tab.clean <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- na.code
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
  ind <- which(tax.tab.h[,2] != na.code)
  ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  ind <- which(tax.tab.h[,3] != na.code)
  ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  ind <- which(tax.tab.h[,4] != na.code)
  ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  ind <- which(tax.tab.h[,5] != na.code)
  ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  ind <- which(tax.tab.h[,6] != na.code)
  ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  ind <- which(tax.tab.h[,7] != na.code)
  ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
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

biom.clean <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
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
    tax.tab <- tax.tab.clean(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

num.tax.rank <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.cleaned <- tax.tab.clean(tax.tab, rem.tax, rem.tax.str, na.code = na.code)
  num.taxa <- c()
  for (i in 1:6) {
    taxa <- unique(tax.tab.cleaned[,i+1])
    uni.taxa <- sum(taxa == na.code)
    num.taxa[i] <- nrow(taxa) - uni.taxa
  }
  return(num.taxa)
}

lib.size.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  lib.size.sum <- c(mean(lib.size), quantile(lib.size))
  names(lib.size.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(lib.size = lib.size, lib.size.sum = lib.size.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

mean.prop.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  mean.prop.sum <- c(mean(mean.prop), quantile(mean.prop))
  names(mean.prop.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(mean.prop = mean.prop, mean.prop.sum = mean.prop.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

rarefy.func <- function(biom, cut.off, multi.rarefy = FALSE) {
  
  if (!multi.rarefy | multi.rarefy == 1) {
    biom <- rarefy_even_depth(biom, cut.off, rngseed = 487)
  } else {
    otu.tab <- otu_table(biom)
    tax.tab <- tax_table(biom)
    tree <- phy_tree(biom)
    sam.dat <- sample_data(biom)
    
    otu.tab.list <- list()
    for (i in 1:multi.rarefy) {
      otu.tab.list[[i]] <- otu_table(rarefy_even_depth(biom, cut.off, rngseed = i), taxa_are_rows = TRUE)
    }
    
    sum.otu.tab <- otu.tab.list[[1]]
    for (i in 2:multi.rarefy) {
      sum.otu.tab <- sum.otu.tab + otu.tab.list[[i]]
    }
    otu.tab <- otu_table(round(sum.otu.tab/multi.rarefy), taxa_are_rows = TRUE)
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat) 
  }
  
  return(biom)
}

#######################################
#     Alpha Diversity Calculation     #
#######################################

alpha.pe.pqe.func <- function(x, tree, norm = TRUE) {
  ind <- which(x != 0)
  s.tree <- prune_taxa(names(x[ind]), tree)
  pe <- AllenH(x[ind], 1, s.tree, Normalize = norm, CheckArguments = FALSE)
  pqe <- AllenH(x[ind], 2, s.tree, Normalize = norm, CheckArguments = FALSE)
  return(c(pe, pqe))
}

alpha.v1.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice, alpha.pd))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  
  return(alpha.div)
}

alpha.v1.func.no.tree <- function(biom) {
  
  otu.tab <- otu_table(biom)
  biom <- merge_phyloseq(round(otu.tab))
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE")
  
  return(alpha.div)
}

alpha.v2.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  alpha.pe.pqe <- t(apply(t.prop.otu.tab, 1, function(x) alpha.pe.pqe.func(x, tree)))
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE")], alpha.ice, alpha.pd, alpha.pe.pqe))
  colnames(alpha.div) <- c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE", "ICE", "PD", "PE", "PQE")
  
  return(alpha.div)
}

#######################################
#          Data Manipulation          #
#######################################

## Which Treatment variable do you want to select?

# Is mon(), sin(), rev, binary, continuous?
is.mon.sin.rev.bin.con <- function(sam.dat) {
  
  n.var <- ncol(sam.dat)
  n.sam <- nrow(sam.dat)
  is.mon <- logical()
  is.rev <- logical()
  is.bin <- logical()
  is.con <- logical()
  is.sin <- logical()
  
  for (i in 1:n.var) {
    sam.var <- as.matrix(sam.dat[,i])
    if (length(table(sam.var)) == 1) {
      is.mon[i] <- TRUE
    }
    if (length(table(sam.var)) != 1) {
      is.mon[i] <- FALSE
    }
    if (length(table(sam.var)) == 2 & any(table(sam.var)==1)) {
      is.sin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2 | !any(table(sam.var)==1)) {
      is.sin[i] <- FALSE
    }
    if (length(table(sam.var)) == n.sam & sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.rev[i] <- TRUE
    }
    if (length(table(sam.var)) != n.sam | sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.rev[i] <- FALSE
    }
    if (length(table(sam.var)) == 2) {
      is.bin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2) {
      is.bin[i] <- FALSE
    }
    if (length(table(sam.var)) != 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- TRUE
    }
    if (length(table(sam.var)) == 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- FALSE
    }
    if (sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.con[i] <- FALSE
    }
    
  }
  return(list(is.mon = is.mon, is.sin = is.sin, is.rev = is.rev, is.bin = is.bin, is.con = is.con))
}


## Extract only Binary or Continuous data

extract.bin.con.func <- function(sam.dat, mon.sin.rev.bin.con) {
  colnames(sam.dat)[(mon.sin.rev.bin.con$is.bin | mon.sin.rev.bin.con$is.con) & !mon.sin.rev.bin.con$is.mon & !mon.sin.rev.bin.con$is.sin]
}


## Is the data Binary or Continuous?

is.bin.con.pri <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind <- which(colnames(sam.dat) == sel.pri.var)
  if(length(ind) != 0){
    if (mon.sin.rev.bin.con$is.bin[ind]) {
      out <- "Binary"
    } else {
      out <- "Continuous"
    }
  }else {
    out = "Neither"
  }
  return(out)
}

## Which outcome do you want to select?

select.outcome.func <- function(treat.var, sel.treat.var) {
  ind = which(treat.var == sel.treat.var)
  return(treat.var[-ind])
}

## Which covariate(s) do you want to select?

select.covariates.func <- function(sam.dat, mon.sin.rev.bin.con, sel.treat.var, sel.out.var) {
  ind.selected <- (colnames(sam.dat) == sel.treat.var | colnames(sam.dat) == sel.out.var)
  ind.mon.sin.rev <- mon.sin.rev.bin.con$is.mon | mon.sin.rev.bin.con$is.rev | mon.sin.rev.bin.con$is.sin
  return(colnames(sam.dat)[!(ind.selected | ind.mon.sin.rev)])
}


alpha.bin.cat.func <- function(sam.dat, sel.bin.var) {
  # bin.var <- unlist(sam.dat[,sel.bin.var])
  # bin.var.no.na <- bin.var[!is.na(bin.var)]
  # bin.cat <- unique(bin.var.no.na)
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  
  return(bin.cat)
}

is.binary <- function(sam.dat, sel.bin.var) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  bin.var.no.na <- bin.var[!is.na(bin.var)]
  bin.cat <- unique(bin.var.no.na)
  if (length(bin.cat) != 2) {
    return(FALSE)
  }
  return(TRUE)
}

alpha.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

alpha.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

alpha.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, alpha.div) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  return(list(bin.var = bin.var, alpha.div = alpha.div))
}

alpha.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x, na.rm = TRUE), quantile(x, na.rm = TRUE))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

alpha.bin.sum.func <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  ref.sum <- matrix(NA, n.alpha, 7)
  com.sum <- matrix(NA, n.alpha, 7)
  print(1)
  for (i in 1:n.alpha) {
    ind.alpha <- alpha.div[,i]
    sum.out <- tapply(ind.alpha, bin.var, alpha.ind.sum.func)
    ref.sum[i,] <- sum.out[[1]]
    com.sum[i,] <- sum.out[[2]]
    print(2)
  }
  rownames(ref.sum) <- colnames(alpha.div)
  colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  rownames(com.sum) <- colnames(alpha.div)
  colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
  names(out) <- levels(bin.var)
  return(out)
}

#######################################
#          Regression Model           #
#######################################

Regression = function(data.na, mediator.list,
                      treatment, outcome, is.covariates, covariates,
                      interaction = c("TRUE", "FALSE"),
                      treatment.type = c("Binary", "Continuous"),
                      outcome.type = c("Binary", "Continuous"),
                      regression.model = c("Probit Regression", "Logistic Regression", "Linear Regression")) {
  
  set.seed(521)
  
  Med.fit <- Out.fit <- Regression.result <- list()
  Mediator.result <- Outcome.result <- data.frame(matrix(nrow=length(mediator.list), ncol=5))
  
  rownames(Mediator.result) <- rownames(Outcome.result) <- mediator.list
  colnames(Mediator.result) <- colnames(Outcome.result) <- c("Estimate", "SE", "P-value", "LLCI", "ULCI")
  
  if (treatment.type == "Binary") {
    data.na[[treatment]] <- as.factor(data.na[[treatment]])
  } else {
    data.na[[treatment]] <- as.numeric(data.na[[treatment]])
  }
  
  for (med in mediator.list) {
    
    if (is.covariates == "None") {
      f1 <<- as.formula(paste(med, "~", treatment, sep=" "))
      Med.fit <- lm(f1, data = data.na)
      
      if (interaction == "TRUE") {
        f2 <<- as.formula(paste(outcome, "~", med, "*", treatment, sep=" "))
      }
      else if (interaction == "FALSE"){
        f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, sep=" "))
      }
    } 
    else {
      f1 <<- as.formula(paste(med, "~", treatment, "+", covariates, sep=" "))
      Med.fit <- lm(f1, data = data.na)
      
      if (interaction == "TRUE") {
        f2 <<- as.formula(paste(outcome, "~", med, "*", treatment, "+", covariates, sep=" "))
      }
      else if (interaction == "FALSE"){
        f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, "+", covariates, sep=" "))
      }
    }
    
    Mediator.result[med, ] <- c(summary(Med.fit)$coefficients[2, -3],
                                confint(Med.fit, rownames(summary(Med.fit)$coefficients)[2], level=0.95))
    
    if (outcome.type == "Binary") {
      
      if (regression.model == "Probit Regression") {
        Out.fit <- glm(f2, data = data.na, family = binomial("probit"))
      }
      else if (regression.model == "Logistic Regression") {
        Out.fit <- glm(f2, data = data.na, family = binomial("logit"))
      }
      
    }
    else if (outcome.type == "Continuous") {
      Out.fit <- lm(f2, data = data.na)
    }
    
    Outcome.result[med, ] <- c(summary(Out.fit)$coefficients[2, -3], confint(Out.fit, med, level=0.95))
    
  }
  
  Regression.result[[1]] <- Mediator.result
  Regression.result[[2]] <- Outcome.result
  
  invisible(Regression.result)
  
}

#######################################
#      Forest Plot for Regression     #
#######################################

Forestplot_Regression = function(regression.result,
                                 result.type = c("Mediator", "Outcome")) {
  
  if (result.type == "Mediator") {
    
    Mean = regression.result[[1]]$Estimate
    SE = regression.result[[1]]$SE
    P_value = regression.result[[1]]$`P-value`
    Lower = regression.result[[1]]$LLCI
    Upper = regression.result[[1]]$ULCI
    
  } else {
    
    Mean = regression.result[[2]]$Estimate
    SE = regression.result[[2]]$SE
    P_value = regression.result[[2]]$`P-value`
    Lower = regression.result[[2]]$LLCI
    Upper = regression.result[[2]]$ULCI
    
  }
  
  Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", round(Mean, 3)))
  SE <- ifelse(abs(round(SE, 3)) == 0, "<.001", sprintf("%.3f", round(SE, 3)))
  P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", round(P_value, 3)))
  
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = rownames(regression.result[[1]]),
                             estimate = as.character(Estimate),
                             se = as.character(SE),
                             p_value = as.character(P_value))
  
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          estimate = "Estimate",
                          se = "SE",
                          p_value = "P-value",
                          summary = TRUE)
  
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  output_df %>% 
    forestplot(labeltext = c(alpha_div, estimate, se, p_value),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1.1)),
               hrzl_lines = list("2" = gpar(lty = 2), 
                                 "11" = gpar(lwd = 1, columns = 1:4, col = "#000044")),
               title = if (result.type == "Mediator") {"Mediator Model"}
               else if (result.type == "Outcome") {"Outcome Model"},
               col = fpColors(box = if (result.type == "Mediator") {"#FFBE3C"}
                              else if (result.type == "Outcome") {"#8f9940"},
                              line = if (result.type == "Mediator") {"golden rod"}
                              else if (result.type == "Outcome") {"darkgreen"}),
               boxsize = 0.15,
               alpha = 0.75)
}

#######################################
#             Sobel Test              #
#######################################

Sobel_Test = function(data.na, mediator.list, treatment, outcome) {
  
  set.seed(521)
  
  Sobel.result <- data.frame(matrix(nrow=length(mediator.list), ncol=2))
  
  rownames(Sobel.result) <- mediator.list
  colnames(Sobel.result) <- c("Z-value", "P-value")
  
  if (class(data.na[[outcome]]) == "factor") {
    data.na[[outcome]] <- as.character(data.na[[outcome]])
  }
  
  for (med in mediator.list) {
    Sobel.result[med, ] <- mediation.test(data.na[, med], data.na[, treatment], data.na[, outcome])[, 1]
  }
  
  Sobel.result[["Z-value"]] <- ifelse(abs(round(Sobel.result[["Z-value"]], 3)) == 0, 
                                      "<.001", sprintf("%.3f", Sobel.result[["Z-value"]]))
  Sobel.result[["P-value"]] <- ifelse(abs(round(Sobel.result[["P-value"]], 3)) == 0, 
                                      "<.001", sprintf("%.3f", Sobel.result[["P-value"]]))
  
  Sobel.result <- cbind("Alpha Diversity" = mediator.list, Sobel.result)
  rownames(Sobel.result) <- 1:length(mediator.list)
  
  invisible(Sobel.result)
  
}

#######################################
#      Forest Plot for Sobel Test     #
#######################################

Forestplot_Sobel = function(regression.result, sobel.result, is.tree) {
  
  Med_mean = regression.result[[1]]$Estimate
  Med_lower = regression.result[[1]]$LLCI
  Med_upper = regression.result[[1]]$ULCI
  Out_mean = regression.result[[2]]$Estimate
  Out_lower = regression.result[[2]]$LLCI
  Out_upper = regression.result[[2]]$ULCI
  
  Mean = Lower = Upper = c()
  for (i in 1:length(Med_mean)) {
    Mean[2*i-1] = Med_mean[i]
    Lower[2*i-1] = Med_lower[i]
    Upper[2*i-1] = Med_upper[i]
  }
  
  for (i in 1:length(Out_mean)) {
    Mean[2*i] = Out_mean[i]
    Lower[2*i] = Out_lower[i]
    Upper[2*i] = Out_upper[i]
  }
  
  Med_estimate <- rep(ifelse(abs(round(Med_mean, 3)) == 0, "<.001", sprintf("%.3f", round(Med_mean, 3))), each = 2)
  Out_estimate <- rep(ifelse(abs(round(Out_mean, 3)) == 0, "<.001", sprintf("%.3f", round(Out_mean, 3))), each = 2)
  P_value <- rep(sobel.result$`P-value`, each=2)
  Alpha_div <- rep(rownames(regression.result[[1]]), each=2)
  
  for (i in 1:length(P_value)) {
    if (i %% 2 == 0) {
      Med_estimate[i] = NA
      P_value[i] = NA
      Alpha_div[i] = NA
    } else {
      Out_estimate[i] = NA
    }
  }
  
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = Alpha_div,
                             p_value = P_value,
                             med_estimate = Med_estimate,
                             out_estimate = Out_estimate)
  
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          p_value = "P-value",
                          med_estimate = "Est. (Med)",
                          out_estimate = "Est. (Out)",
                          summary = TRUE)
  
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  
  box_list <- list()
  for (i in 1:(nrow(output_df))){
    if(i %% 2 == 1){
      box_list[[i]] <- gpar(fill = "#8f9940", col = "#8f9940")
    }else{
      box_list[[i]] <- gpar(fill = "#FFBE3C", col = "#FFBE3C")
    }
  }
  
  styles <- fpShapesGp(box = box_list)
  
  output_df %>% 
    forestplot(labeltext = c(alpha_div, p_value, med_estimate, out_estimate),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1),
                                label = gpar(fontfamily="", cex=1.15)),
               hrzl_lines = if (is.tree == "withTree") {list("2" = gpar(lty = 2), 
                                                             "20" = gpar(lwd = 1, columns = 1:4, col = "#000044"))} 
                            else if (is.tree == "withoutTree") {list("2" = gpar(lty = 2), 
                                                                     "18" = gpar(lwd = 1, columns = 1:4, col = "#000044"))},
               col = fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"),
               boxsize = 0.15,
               alpha = 0.75,
               shapes_gp = styles)
}

#######################################
# Preacher and Hayes Bootstrap Method #
#######################################

Hayes_Mediation = function(data.na, mediator.list,
                           treatment, outcome, outcome.type,
                           is.covariates, hayes.covariates, n.sim) {
  
  set.seed(521)
  
  Indirect.result <- data.frame(matrix(nrow=length(mediator.list), ncol=4))
  Direct.result <- Total.result <- data.frame(matrix(nrow=length(mediator.list), ncol=6))
  
  rownames(Indirect.result) <- rownames(Direct.result) <- rownames(Total.result) <- mediator.list
  colnames(Indirect.result) <- c("Effect", "BootSE", "BootLLCI", "BootULCI")
  colnames(Direct.result) <- colnames(Total.result) <- c("Effect", "SE", "z", "p", "LLCI", "ULCI")
  
  data.na[[treatment]] <- as.numeric(data.na[[treatment]])
  data.na[[outcome]] <- as.numeric(data.na[[outcome]])
  
  if (outcome.type == "Binary") {
    
    for (med in mediator.list) {
      
      incProgress(1/10, message = paste0("Calculating: ", med))
      
      if (is.covariates == "None") {
        result.saved <- process(data=data.na, y=outcome, x=treatment, m=med, boot=n.sim, model=4, total=1, save=2)
      }
      else {
        
        for (cov in hayes.covariates) {
          data.na[[cov]] <- as.numeric(data.na[[cov]])
        }
        
        result.saved <- process(data=data.na, y=outcome, x=treatment, m=med, boot=n.sim, cov=hayes.covariates,
                                model=4, total=1, save=2)
      }
      
      Indirect.result[med, ] <- result.saved[nrow(result.saved), 1:4]
      Direct.result[med, ]<- result.saved[nrow(result.saved)-1, 1:6]
      
    }
    
    Hayes.result <- list(Indirect = Indirect.result, Direct = Direct.result)
    
  } else if (outcome.type == "Continuous") {
    
    for (med in mediator.list) {
      
      incProgress(1/10, message = paste0("Calculating: ", med))
      
      if (is.covariates == "None") {
        result.saved <- process(data=data.na, y=outcome, x=treatment, m=med, boot=n.sim, model=4, total=1, save=2)
      }
      else {
        
        for (cov in hayes.covariates) {
          data.na[[cov]] <- as.numeric(data.na[[cov]])
        }
        
        result.saved <- process(data=data.na, y=outcome, x=treatment, m=med, boot=n.sim, cov=hayes.covariates,
                                model=4, total=1, save=2)
      }
      
      Indirect.result[med, ] <- result.saved[nrow(result.saved), 1:4]
      Direct.result[med, ]<- result.saved[nrow(result.saved)-1, 1:6]
      
    }
    
    Hayes.result <- list(Indirect = Indirect.result, Direct = Direct.result)
    
  }
  
  invisible(Hayes.result)
  
}

#######################################
#    Forest Plot for Hayes' Method    #
#######################################

Forestplot_Hayes <- function(hayes.result, result.type = c("Indirect", "Direct"), is.tree) {
  
  Mean = c(); SE = c(); Lower = c(); Upper = c(); P_value = c()
  
  if (result.type == "Indirect") {
    
    Mean = hayes.result[[1]]$Effect
    SE = hayes.result[[1]]$BootSE
    Lower = hayes.result[[1]]$BootLLCI
    Upper = hayes.result[[1]]$BootULCI
    
    Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", Mean))
    
    output_base_data <- tibble(mean = Mean,
                               lower = Lower,
                               upper = Upper,
                               alpha_div = rownames(hayes.result[[1]]),
                               estimate = as.character(Estimate))
    
    output_header <- tibble(alpha_div = "Alpha Diversity",
                            estimate = "Estimate",
                            summary = TRUE)
    
    output_df <- bind_rows(output_header,
                           output_base_data)
    
    output_df %>% 
      forestplot(labeltext = c(alpha_div, estimate),
                 is.summary = summary,
                 xlab = "95% Confidence Interval",
                 txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                  ticks = gpar(cex=0.75),
                                  xlab = gpar(cex = 1),
                                  label = gpar(fontfamily="", cex=1.15)),
                 hrzl_lines = if (is.tree == "withTree") {list("2" = gpar(lty = 2), 
                                                               "11" = gpar(lwd = 1, columns = 1:3, col = "#000044"))} 
                              else if (is.tree == "withoutTree") {list("2" = gpar(lty = 2), 
                                                               "10" = gpar(lwd = 1, columns = 1:3, col = "#000044"))},
                 title = "Indirect (Mediation) Effect",
                 col = fpColors(box = "#D55E00", line="black"),
                 boxsize = 0.15,
                 alpha = 0.75)
    
  } else if (result.type == "Direct") {
    
    Mean = hayes.result[[2]]$Effect
    SE = hayes.result[[2]]$SE
    P_value = hayes.result[[2]]$p
    Lower = hayes.result[[2]]$LLCI
    Upper = hayes.result[[2]]$ULCI
    
    Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", Mean))
    P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", P_value))
    
    output_base_data <- tibble(mean = Mean,
                               lower = Lower,
                               upper = Upper,
                               alpha_div = rownames(hayes.result[[2]]),
                               estimate = as.character(Estimate),
                               p_value = as.character(P_value))
    
    output_header <- tibble(alpha_div = "Alpha Diversity",
                            estimate = "Estimate",
                            p_value = "P-value",
                            summary = TRUE)
    
    output_df <- bind_rows(output_header,
                           output_base_data)
    
    output_df %>% 
      forestplot(labeltext = c(alpha_div, estimate, p_value),
                 is.summary = summary,
                 xlab = "95% Confidence Interval",
                 txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                  ticks = gpar(cex=0.75),
                                  xlab = gpar(cex = 1),
                                  label = gpar(fontfamily="", cex=1.15)),
                 hrzl_lines = if (is.tree == "withTree") {list("2" = gpar(lty = 2), 
                                                               "11" = gpar(lwd = 1, columns = 1:3, col = "#000044"))} 
                              else if (is.tree == "withoutTree") {list("2" = gpar(lty = 2), 
                                                                       "10" = gpar(lwd = 1, columns = 1:3, col = "#000044"))},
                 title = "Direct Effect",
                 col = fpColors(box = "#084081", line="black"),
                 boxsize = 0.15,
                 alpha = 0.75)
    
  }

}

########################################
# Divide-Aggregate Composite-Null Test #
########################################

DACT_Mediation = function(regression.result, mediator.list) {
  
  DACT.result <- as.data.frame(matrix(nrow=length(mediator.list), ncol=2))
  colnames(DACT.result) <- c("Alpha Diversity", "P-value")
  
  p_a <- regression.result[[1]]$`P-value`
  p_b <- regression.result[[2]]$`P-value`
  
  DACT.p <- DACT(p_a, p_b, correction="JC")
  
  DACT.result[["Alpha Diversity"]] <- mediator.list
  DACT.result[["P-value"]] <- ifelse(abs(round(DACT.p, 3)) == 0, "<.001", sprintf("%.3f", DACT.p))
  
  invisible(DACT.result)
  
}

#######################################
#         Forest Plot for DACT        #
#######################################

Forestplot_DACT = function(regression.result, DACT.result, is.tree) {
  
  Med_mean = regression.result[[1]]$Estimate
  Med_lower = regression.result[[1]]$LLCI
  Med_upper = regression.result[[1]]$ULCI
  Out_mean = regression.result[[2]]$Estimate
  Out_lower = regression.result[[2]]$LLCI
  Out_upper = regression.result[[2]]$ULCI
  
  Mean = Lower = Upper = c()
  for (i in 1:length(Med_mean)) {
    Mean[2*i-1] = Med_mean[i]
    Lower[2*i-1] = Med_lower[i]
    Upper[2*i-1] = Med_upper[i]
  }
  
  for (i in 1:length(Out_mean)) {
    Mean[2*i] = Out_mean[i]
    Lower[2*i] = Out_lower[i]
    Upper[2*i] = Out_upper[i]
  }
  
  Med_estimate <- rep(ifelse(abs(round(Med_mean, 3)) == 0, "<.001", sprintf("%.3f", round(Med_mean, 3))), each = 2)
  Out_estimate <- rep(ifelse(abs(round(Out_mean, 3)) == 0, "<.001", sprintf("%.3f", round(Out_mean, 3))), each = 2)
  
  for (i in 1:nrow(DACT.result)) {
    if (is.na(DACT.result$`P-value`[i])) {
      DACT.result$`P-value`[i] = "NA"
    }
  }
  
  P_value <- rep(DACT.result$`P-value`, each=2)
  Alpha_div <- rep(rownames(regression.result[[1]]), each=2)
  
  for (i in 1:length(P_value)) {
    if (i %% 2 == 0) {
      Med_estimate[i] = NA
      P_value[i] = NA
      Alpha_div[i] = NA
    } else {
      Out_estimate[i] = NA
    }
  }
  
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = Alpha_div,
                             p_value = P_value,
                             med_estimate = Med_estimate,
                             out_estimate = Out_estimate)
  
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          p_value = "P-value",
                          med_estimate = "Est. (Med)",
                          out_estimate = "Est. (Out)",
                          summary = TRUE)
  
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  
  box_list <- list()
  for (i in 1:(nrow(output_df))){
    if(i %% 2 == 1){
      box_list[[i]] <- gpar(fill = "#8f9940", col = "#8f9940")
    }else{
      box_list[[i]] <- gpar(fill = "#FFBE3C", col = "#FFBE3C")
    }
  }
  
  styles <- fpShapesGp(box = box_list)
  
  output_df %>% 
    forestplot(labeltext = c(alpha_div, p_value, med_estimate, out_estimate),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1),
                                label = gpar(fontfamily="", cex=1.15)),
               hrzl_lines = if (is.tree == "withTree") {list("2" = gpar(lty = 2), 
                                                             "20" = gpar(lwd = 1, columns = 1:4, col = "#000044"))} 
                            else if (is.tree == "withoutTree") {list("2" = gpar(lty = 2), 
                                                                     "18" = gpar(lwd = 1, columns = 1:4, col = "#000044"))},
               col = fpColors(box=rgb(1,0,0), line="black"),
               boxsize = 0.15,
               alpha = 0.75,
               shapes_gp = styles)
  
}

########################################
#     Potential Outcomes Framework     #
########################################

Mediation = function(data.na, mediator.list,
                     treatment, outcome, is.covariates, covariates,
                     interaction = c("TRUE", "FALSE"),
                     treatment.type = c("Binary", "Continuous"),
                     outcome.type = c("Binary", "Continuous"),
                     regression.model = c("Probit Regression", "Logistic Regression", "Linear Regression"),
                     method = c("Quasi-Bayesian", "Bootstrap"), 
                     boot.ci.type = c("BCa CI", "Percentile CI (Default)"), n.sim) {
  
  set.seed(521)
  
  ## Data Manipulation
  if (treatment.type == "Binary") {
    data.na[[treatment]] <- as.factor(data.na[[treatment]])
  } else {
    data.na[[treatment]] <- as.numeric(data.na[[treatment]])
  }
  
  if (outcome.type == "Binary") {
    data.na[[outcome]] <- as.factor(data.na[[outcome]])
  } else {
    data.na[[outcome]] <- as.numeric(data.na[[outcome]])
  }
  
  print(class(data.na[[outcome]]))
  
  ## Potential Outcomes Framework Mediation Analysis
  Med.fit <- list(); Out.fit <- list(); Med.out <- list()
  
  for (med in mediator.list) {
    
    incProgress(1/10, message = paste0("Calculating: ", med))
    
    if (is.covariates == "None") {
      f1 <<- as.formula(paste(med, "~", treatment, sep=" "))
      Med.fit <- lm(f1, data = data.na)
      
      if (interaction == "TRUE") {
        f2 <<- as.formula(paste(outcome, "~", med, "*", treatment, sep=" "))
        # f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, "+", paste(med, treatment, sep = ":"), sep=" "))
      }
      else if (interaction == "FALSE"){
        f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, sep=" "))
      }
    } 
    else {
      f1 <<- as.formula(paste(med, "~", treatment, "+", covariates, sep=" "))
      Med.fit <- lm(f1, data = data.na)
      
      if (interaction == "TRUE") {
        f2 <<- as.formula(paste(outcome, "~", med, "*", treatment, "+", covariates, sep=" "))
        # f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, "+", paste(med, treatment, sep = ":"), "+", covariates, sep=" "))
      }
      else if (interaction == "FALSE"){
        f2 <<- as.formula(paste(outcome, "~", med, "+", treatment, "+", covariates, sep=" "))
      }
    }
    
    outcome.type = match.arg(outcome.type)
    
    if (outcome.type == "Binary") {
      
      if (regression.model == "Probit Regression") {
        Out.fit <- glm(f2, data = data.na, family = binomial("probit"))
      }
      else if (regression.model == "Logistic Regression") {
        Out.fit <- glm(f2, data = data.na, family = "binomial")
      } 
      
    }
    else if (outcome.type == "Continuous") {
      Out.fit <- lm(f2, data = data.na)
    }
    
    method = match.arg(method)
    
    if (method == "Quasi-Bayesian") {
      Med.out[[med]] <- mediate(Med.fit, Out.fit, treat = treatment, mediator = med,
                                robustSE = TRUE, sims = n.sim)
    }
    else if (method == "Bootstrap") {
      if (boot.ci.type == "Percentile CI (Default)") {
        Med.out[[med]] <- mediate(Med.fit, Out.fit, treat = treatment, mediator = med,
                                  boot = TRUE, boot.ci.type = "perc", sims = n.sim)
      }
      else if (boot.ci.type == "BCa CI") {
        Med.out[[med]] <- mediate(Med.fit, Out.fit, treat = treatment, mediator = med,
                                  boot = TRUE, boot.ci.type = "bca", sims = n.sim)
      }
    }
  }
  
  invisible(Med.out)
  
}


###############################################
# Forestplot for Potential Outcomes Framework #
###############################################

Forestplot = function(med.out, is.tree,
                      outcome.type = c("Binary", "Continuous"),
                      treatment.type =  c("Binary", "Continuous"),
                      interaction = c("TRUE", "FALSE"),
                      effect = c("ACME (Average)", "ADE", "Total Effect",
                                 "ACME (Control)", "ACME (Treated)")) {
  
  Mean = c(); Lower = c(); Upper = c(); P_value = c()
  
  if (effect == "ACME (Average)") {
    for (i in 1:length(med.out)) {
      Mean[i] = summary(med.out[[i]])$d.avg
      Lower[i] = summary(med.out[[i]])$d.avg.ci[1]
      Upper[i] = summary(med.out[[i]])$d.avg.ci[2]
      P_value[i] = summary(med.out[[i]])$d.avg.p
    }
  }
  else if (effect == "ADE") {
    for (i in 1:length(med.out)) {
      Mean[i] = summary(med.out[[i]])$z.avg
      Lower[i] = summary(med.out[[i]])$z.avg.ci[1]
      Upper[i] = summary(med.out[[i]])$z.avg.ci[2]
      P_value[i] = summary(med.out[[i]])$z.avg.p
    }
  }
  else if (effect == "Total Effect") {
    for (i in 1:length(med.out)) {
      Mean[i] = summary(med.out[[i]])$tau.coef
      Lower[i] = summary(med.out[[i]])$tau.ci[1]
      Upper[i] = summary(med.out[[i]])$tau.ci[2]
      P_value[i] = summary(med.out[[i]])$tau.p
    }
  } 
  
  outcome.type = match.arg(outcome.type)
  
  if (treatment.type == "Binary") {
    
    if (effect == "ACME (Control)") {
      for (i in 1:length(med.out)) {
        Mean[i] = summary(med.out[[i]])$d0
        Lower[i] = summary(med.out[[i]])$d0.ci[1]
        Upper[i] = summary(med.out[[i]])$d0.ci[2]
        P_value[i] = summary(med.out[[i]])$d0.p
      }
    }
    else if (effect == "ACME (Treated)") {
      for (i in 1:length(med.out)) {
        Mean[i] = summary(med.out[[i]])$d1
        Lower[i] = summary(med.out[[i]])$d1.ci[1]
        Upper[i] = summary(med.out[[i]])$d1.ci[2]
        P_value[i] = summary(med.out[[i]])$d1.p
      }
    }
    
  }
  
  else if (treatment.type == "Continuous") {
    
    if (effect == "ACME (Control)" | effect == "ACME (Treated)") {
      stop(paste(ifelse(effect == "ACME (Control)", "ACME (Control)", "ACME (Treated)"), 
                 "is not provided for treatment.type = \"Continuous\""))
    }
  }
  
  Estimate <- ifelse(abs(round(Mean, 3)) == 0, "<.001", sprintf("%.3f", Mean))
  P_value <- ifelse(abs(round(P_value, 3)) == 0, "<.001", sprintf("%.3f", P_value))
  
  output_base_data <- tibble(mean = Mean,
                             lower = Lower,
                             upper = Upper,
                             alpha_div = names(med.out),
                             estimate = as.character(Estimate),
                             p_value = as.character(P_value))
  
  output_header <- tibble(alpha_div = "Alpha Diversity",
                          estimate = "Estimate",
                          p_value = "P-value",
                          summary = TRUE)
  
  output_df <- bind_rows(output_header,
                         output_base_data)
  
  output_df %>% 
    forestplot(labeltext = c(alpha_div, estimate, p_value),
               is.summary = summary,
               xlab = "95% Confidence Interval",
               txt_gp = fpTxtGp(title = gpar(cex=1.35),
                                ticks = gpar(cex=0.75),
                                xlab = gpar(cex = 1),
                                label = gpar(fontfamily="", cex=1.15)),
               hrzl_lines = if (is.tree == "withTree") {list("2" = gpar(lty = 2), 
                                                             "11" = gpar(lwd = 1, columns = 1:3, col = "#000044"))} 
                            else if (is.tree == "withoutTree") {list("2" = gpar(lty = 2), 
                                                                     "10" = gpar(lwd = 1, columns = 1:3, col = "#000044"))},
               title = effect,
               col = fpColors(box = if (effect == "ACME (Average)") {"#D55E00"}
                              else if (effect == "ACME (Control)") {"#f0b135"}
                              else if (effect == "ACME (Treated)") {"orange2"}
                              else if (effect == "ADE") {"#084081"}
                              else if (effect == "Total Effect") {"#046B2D"},
                              line="black"),
               boxsize = 0.15,
               alpha = 0.75)
  
}

########################
# Sensitivity Analysis #
########################

Sensitivity_analysis <- function(med.out, n.sims = 100) {
  
  mediator.list <- names(med.out)
  
  sens.out <- list(); sens.plot <- list()
  
  for (med in mediator.list) {
    sens.out[[med]] <- medsens(med.out[[med]], rho.by = 0.1, effect.type = "indirect", sims = 100)
    
    par(mfrow = c(2,2))
    
    plot(sens.out[[med]], sens.par = "rho", main = med, ylim = c(-0.2, 0.2))
    plot(sens.out[[med]], sens.par = "R2", r.type = "total", sign.prod = "positive")
    sens.plot[[med]] <- recordPlot()
  }
  
  dev.off()
  
  invisible(sens.plot)
  
}


########################
#        Theme         #
########################

### creating custom theme object
customTheme <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(0,0,0)"
  ,primaryFontColor = "rgb(0,0,0)"
  ,infoFontColor = "rgb(0,0,0)"
  ,successFontColor = "rgb(0,0,0)"
  ,warningFontColor = "rgb(0,0,0)"
  ,dangerFontColor = "rgb(0,0,0)"
  ,bodyBackColor = "rgb(255,255,255)"
  
  ### header
  ,logoBackColor = "rgb(35, 49, 64)"
  
  ,headerButtonBackColor = "rgb(238,238,238)"
  ,headerButtonIconColor = "rgb(75,75,75)"
  ,headerButtonBackColorHover = "rgb(210,210,210)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(238,238,238)"
  ,headerBoxShadowColor = "#aaaaaa"
  ,headerBoxShadowSize = "0px 0px 0px"
  
  ### sidebar
  ,sidebarBackColor = "rgb(52,62,72)"
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = ""
  ,sidebarShadowColor = "0px 0px 0px"
  
  ,sidebarUserTextColor = "rgb(35, 49, 64)"
  
  ,sidebarSearchBackColor = "rgb(255, 255, 255)"
  ,sidebarSearchIconColor = "rgb(35, 49, 64)"
  ,sidebarSearchBorderColor = "rgb(35, 49, 64)"
  
  ,sidebarTabTextColor = "rgb(205,205,205)"
  ,sidebarTabTextSize = 14
  ,sidebarTabBorderStyle = "none"
  ,sidebarTabBorderColor = "none"
  ,sidebarTabBorderWidth = 0
  
  ,sidebarTabBackColorSelected = "rgb(70,80,90)"
  ,sidebarTabTextColorSelected = "rgb(255,255,255)"
  ,sidebarTabRadiusSelected = "0px"
  
  ,sidebarTabBackColorHover = "rgb(55,65,75)"
  ,sidebarTabTextColorHover = "rgb(255,255,255)"
  ,sidebarTabBorderStyleHover = "none"
  ,sidebarTabBorderColorHover = "none"
  ,sidebarTabBorderWidthHover = 0
  ,sidebarTabRadiusHover = "0px"
  
  ### boxes
  ,boxBackColor = "rgb(245,245,245)"
  ,boxBorderRadius = 3
  ,boxShadowSize = "0px 0px 0px"
  ,boxShadowColor = "rgba(0,0,0,0)"
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(210,214,220)"
  ,boxPrimaryColor = "rgb(35, 49, 64)"
  ,boxInfoColor = "rgb(210,214,220)"
  # Progress Bar Color
  ,boxSuccessColor = "rgb(112,173,71)"
  ,boxWarningColor = "rgb(244,156,104)"
  ,boxDangerColor = "rgb(255,88,55)"
  
  ,tabBoxTabColor = "rgb(255,255,255)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(0,0,0)"
  ,tabBoxTabTextColorSelected = "rgb(35, 49, 64)"
  ,tabBoxBackColor = "rgb(255,255,255)"
  ,tabBoxHighlightColor = "rgb(44, 59, 65)"
  ,tabBoxBorderRadius = 0
  
  ### inputs
  ,buttonBackColor = "rgb(245,245,245)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(35, 49, 64)"
  ,buttonBorderRadius = 3
  
  ,buttonBackColorHover = "rgb(227,227,227)"
  ,buttonTextColorHover = "rgb(100,100,100)"
  ,buttonBorderColorHover = "rgb(200,200,200)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(200,200,200)"
  ,textboxBorderRadius = 0
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(200,200,200)"
  
  ### tables
  ,tableBackColor = "rgb(255, 255, 255)"
  ,tableBorderColor = "rgb(245, 245, 245)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)
