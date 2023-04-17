######################################
# Quality control and transformation #
######################################

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
library(lme4)
library(lmerTest)
library(dirmult) 
library(robustbase)
library(robCompositions) 
library(BiasedUrn)
library(CompQuadForm)
library(GUniFrac) 
library(ecodist) 
library(MiRKAT)
library(GLMMMiRKAT)
library(proxy)
library(matrixStats)

#######################################
#          Data Manipulation          #
#######################################

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  return(bin.cat)
}

beta.bin.cat.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

beta.bin.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  
  return(sam.dat)
}

beta.con.recode.func <- function(sam.dat, sel.con.var, rename.con.var) {
  sam.dat[[sel.con.var]] <- as.numeric(sam.dat[[sel.con.var]])
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  
  return(sam.dat)
}

Renaming = function(data.na, treatment, outcome, treatment.type, outcome.type,
                    bin_ref1 = NA, bin_com1 = NA, 
                    bin_ref2 = NA, bin_com2 = NA,
                    con1 = NA, con2 = NA) {
  
  if (treatment.type == "Binary") {
    beta.bin.ori.cat_treatvar <- beta.bin.cat.ori.func(data.na, treatment)
    
    beta.Data_treat <- beta.bin.recode.func(data.na, treatment, beta.bin.ori.cat_treatvar,
                                            bin_ref1, bin_com1)
    
    if (outcome.type == "Binary") {
      beta.bin.ori.cat_outvar <- beta.bin.cat.ori.func(data.na, outcome)
      
      beta.Data <- beta.bin.recode.func(beta.Data_treat, outcome,
                                        beta.bin.ori.cat_outvar,
                                        bin_ref2, bin_com2)
      
    } else if (outcome.type == "Continuous") {
      
      beta.Data <- beta.con.recode.func(beta.Data_treat, outcome, con2)
      
    }
    
  } else if (treatment.type == "Continuous") {
    beta.Data_treat <- beta.con.recode.func(data.na, treatment, con1)
    
    if (outcome.type == "Binary") {
      beta.bin.ori.cat_outvar <- beta.bin.cat.ori.func(data.na, outcome)
      
      beta.Data <- beta.bin.recode.func(beta.Data_treat, outcome,
                                        beta.bin.ori.cat_outvar,
                                        bin_ref2, bin_com2)
      
    } else if (outcome.type == "Continuous") {
      
      beta.Data <- beta.con.recode.func(beta.Data_treat, outcome, con2)
      
    }
    
  }
  
  invisible(beta.Data)
  
}

#######################################
#           Beta Diversity            #
#######################################

Ds.Ks.func <- function(rare.biom, biom.after.qc, is.tree) {
  
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  
  if (is.tree == "withTree") {
    
    no.rare.tree <- phy_tree(biom.after.qc)
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
    u.unif <- unifs[, , "d_UW"]
    g.unif <- unifs[, , "d_0.5"]
    w.unif <- unifs[, , "d_1"]
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    u.unif.k <- D2K(u.unif)
    g.unif.k <- D2K(g.unif)
    w.unif.k <- D2K(w.unif)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
    rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
    rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
    
    return(
      list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
           Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k))
    )
  } 
  else if (is.tree == "withoutTree") {
    
    jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
    bc <- as.matrix(bcdist(t(rare.otu.tab)))
    
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    
    rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
    rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
    
    return(
      list(Ds = list(Jaccard = jac, Bray.Curtis = bc),
           Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k))
    )
  }
  
}


beta.bin.out.func <- function(sam.dat, Ds.Ks, sel.bin.var, sel.ref, sel.com) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.ref)
  ind.com <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  
  return(list(bin.var = bin.var, Ds = Ds, Ks = Ks))
}

beta.con.out.func <- function(sam.dat, Ds.Ks, sel.con.var, rename.con.var, ind.missing) {
  con.var <- as.data.frame(as.matrix(sam.dat[[rename.con.var]]))
  colnames(con.var) <- rename.con.var
  
  if ("NA" %in% ind.missing) {
    Ds <- Ds.Ks$Ds
    Ks <- Ds.Ks$Ks
  } else {
    Ds <- Ds.Ks$Ds
    Ks <- Ds.Ks$Ks
    for (i in 1:length(Ds)) {
      Ds[[i]] <- Ds[[i]][-ind.missing, -ind.missing]
    }
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][-ind.missing, -ind.missing]
    }
  }
  
  return(list(con.var = con.var, Ds = Ds, Ks = Ks))
}

#######################################
#               MedTest               #
#######################################

MedTest = function(data.na, beta.mediator.na, medtest.covariates, treatment, outcome, n.perm) {
  
  set.seed(521)
  cov = c()
  
  if (length(medtest.covariates) == 0) {
    
    out <- MedOmniTest(x=as.numeric(data.na[[treatment]]), y=as.numeric(data.na[[outcome]]), m.list=beta.mediator.na, 
                       z=NULL, nperm = n.perm)
    
  } else if (length(medtest.covariates) == 1) {
    
    out <- MedOmniTest(x=as.numeric(data.na[[treatment]]), y=as.numeric(data.na[[outcome]]), m.list=beta.mediator.na, 
                       z=as.numeric(data.na[[medtest.covariates]]), nperm = n.perm)
    print(out)
    
  } else if (length(medtest.covariates) > 1) {
    
    for (i in 1:length(medtest.covariates)) {
      cov = cbind(cov, as.numeric(data.na[[medtest.covariates[[i]]]]))
    }
    print(head(cov))
    print(ncol(cov))
    print(class(cov))
    
    out <- MedOmniTest(x=as.numeric(data.na[[treatment]]), y=as.numeric(data.na[[outcome]]), m.list=beta.mediator.na, 
                       z=cov, nperm = n.perm)
  }
  
  return(out)
  
}

#######################################
#            Visualization            #
#######################################

## Without tree --------------------------------------------------------------------------------------

MedTest.bin.bin.plot = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    if (med.test.out$margPs[i] < 0.05) {
      sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    if (med.test.out$margPs[i] >= 0.05) {
      sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    
    treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
    out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
    
    plot(treat.mod, ellipse = TRUE, hull = FALSE, 
         main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
         col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
    plot(out.mod, ellipse = TRUE, hull = FALSE, 
         main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
         col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("left", title = "Treatment", legend = levels(beta.treatvar.out$bin.var),
         fil = c("golden rod", "darkgreen", cex=2.0, box.lty=0), bty = "n", cex=1.2)
  legend("right", title = "Outcome", legend = levels(beta.outvar.out$bin.var),
         fil = c("orangered3", "royalblue4", cex=2.0, box.lty=0), bty = "n", cex=1.2)
  legend("bottom right", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.2)
  
}

MedTest.con.bin.plot = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    if (med.test.out$margPs[i] < 0.05) {
      sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    if (med.test.out$margPs[i] >= 0.05) {
      sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    
    treat.mod <- cmdscale(beta.treatvar.out$Ds[[i]])
    out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
    
    temp.col <- colorRampPalette(c("yellow2", "darkgreen"))
    color.palette <- temp.col(nrow(treat.mod))
    
    plot(treat.mod, ellipse = FALSE, hull = FALSE, 
         main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
         col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
    plot(out.mod, ellipse = TRUE, hull = FALSE, 
         main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
         col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image <- as.raster(matrix(color.palette, ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=0.45, y = seq(0.37,0.52,l=2), labels = seq(min(beta.treatvar.out$con.var), max(beta.treatvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image, 0.1, 0.34, 0.25, 0.55)
  text(x=0.3, y=0.66, labels = "Treatment", cex=1.6)
  
  legend("right", legend = levels(beta.outvar.out$bin.var), title = "Outcome",
         fil = c("orangered3", "royalblue4", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)
  
}

MedTest.bin.con.plot = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    if (med.test.out$margPs[i] < 0.05) {
      sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    if (med.test.out$margPs[i] >= 0.05) {
      sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    
    treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
    out.mod <- cmdscale(beta.outvar.out$Ds[[i]])
    
    temp.col <- colorRampPalette(c("red2", "darkblue"))
    color.palette <- temp.col(nrow(out.mod))
    
    plot(treat.mod, ellipse = TRUE, hull = FALSE, 
         main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
         col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
    plot(out.mod, ellipse = TRUE, hull = FALSE, 
         main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
         col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image <- as.raster(matrix(color.palette, ncol=1))
  plot(c(0,4),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=3.65, y = seq(0.37,0.52,l=2), labels = seq(min(beta.outvar.out$con.var), max(beta.outvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image, 2.65, 0.34, 2.95, 0.55)
  text(x=3.4, y=0.66, labels = "Outcome", cex=1.6)
  
  legend("left", legend = levels(beta.treatvar.out$bin.var), title = "Treatment",
         fil = c("golden rod", "darkgreen", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)
  
}

MedTest.con.con.plot = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    if (med.test.out$margPs[i] < 0.05) {
      sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    if (med.test.out$margPs[i] >= 0.05) {
      sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
    }
    
    treat.mod <- cmdscale(beta.treatvar.out$Ds[[i]])
    out.mod <- cmdscale(beta.outvar.out$Ds[[i]])
    
    temp.col1 <- colorRampPalette(c("yellow2", "darkgreen"))
    color.palette1 <- temp.col1(nrow(treat.mod))
    temp.col2 <- colorRampPalette(c("red2", "darkblue"))
    color.palette2 <- temp.col2(nrow(out.mod))
    
    plot(treat.mod, ellipse = TRUE, hull = FALSE, 
         main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
         col = color.palette1, mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
    plot(out.mod, ellipse = TRUE, hull = FALSE, 
         main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
         col = color.palette2, mgp=c(2.5,1,0), pch=c(16,17),
         cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image1 <- as.raster(matrix(color.palette1, ncol=1))
  plot(c(0,4),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=0.8, y = seq(0.37,0.52,l=2), labels = seq(min(beta.treatvar.out$con.var), max(beta.treatvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image1, 0.1, 0.34, 0.4, 0.55)
  text(x=0.6, y=0.66, labels = "Treatment", cex=1.6)
  
  legend_image2 <- as.raster(matrix(color.palette2, ncol=1))
  text(x=3.65, y = seq(0.37,0.52,l=2), labels = seq(min(beta.outvar.out$con.var), max(beta.outvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image2, 2.65, 0.34, 2.95, 0.55)
  text(x=3.4, y=0.66, labels = "Outcome", cex=1.6)
  
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)
  
}



## With tree ----------------------------------------------------------------------------------------

MedTest.bin.bin.plot1 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
  
    # Jaccard, U.UniFrac, W.UniFrac
    if (i %% 2 == 1) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
      out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.5, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4) 
      
    }
  }
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.35)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[3], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[5], side = 3, line = -50.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[3], side = 3, line = -72.5, outer = TRUE, cex = 1.05)
  
}

MedTest.bin.bin.plot2 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Bray.Curtis, G.UniFrac
    if (i %% 2 == 0) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
      out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
    }
  }
  
  print(sub.title)
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[4], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("left", title = "Treatment", legend = levels(beta.treatvar.out$bin.var),
         fil = c("golden rod", "darkgreen", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("right", title = "Outcome", legend = levels(beta.outvar.out$bin.var),
         fil = c("orangered3", "royalblue4", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)
  
}

MedTest.con.bin.plot1 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Jaccard, U.UniFrac, W.UniFrac
    if (i %% 2 == 1) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- cmdscale(beta.treatvar.out$Ds[[i]])
      out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
      
      temp.col <- colorRampPalette(c("yellow2", "darkgreen"))
      color.palette <- temp.col(nrow(treat.mod))
      
      plot(treat.mod, ellipse = FALSE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
    }
  }
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.35)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[3], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[5], side = 3, line = -50.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[3], side = 3, line = -72.5, outer = TRUE, cex = 1.05)
  
}

MedTest.con.bin.plot2 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Bray.Curtis, G.UniFrac
    if (i %% 2 == 0) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- cmdscale(beta.treatvar.out$Ds[[i]])
      out.mod <- betadisper(as.dist(beta.outvar.out$Ds[[i]]), beta.outvar.out$bin.var)
      
      temp.col <- colorRampPalette(c("yellow2", "darkgreen"))
      color.palette <- temp.col(nrow(treat.mod))
      
      plot(treat.mod, ellipse = FALSE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("orangered3", "royalblue4"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
    }
  }
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[4], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image <- as.raster(matrix(color.palette, ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=0.45, y = seq(0.37,0.52,l=2), labels = seq(min(beta.treatvar.out$con.var), max(beta.treatvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image, 0.1, 0.34, 0.25, 0.55)
  text(x=0.3, y=0.66, labels = "Treatment", cex=1.6)

  legend("right", legend = levels(beta.outvar.out$bin.var), title = "Outcome",
         fil = c("orangered3", "royalblue4", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)

}

MedTest.bin.con.plot1 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Jaccard, U.UniFrac, W.UniFrac
    if (i %% 2 == 1) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
      out.mod <- cmdscale(beta.outvar.out$Ds[[i]])
      
      temp.col <- colorRampPalette(c("red2", "darkblue"))
      color.palette <- temp.col(nrow(out.mod))
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
    }
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.35)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[3], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[5], side = 3, line = -50.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[3], side = 3, line = -72.5, outer = TRUE, cex = 1.05)
  
}

MedTest.bin.con.plot2 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Bray.Curtis, G.UniFrac
    if (i %% 2 == 0) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- betadisper(as.dist(beta.treatvar.out$Ds[[i]]), beta.treatvar.out$bin.var)
      out.mod <- cmdscale(beta.outvar.out$Ds[[i]])
      
      temp.col <- colorRampPalette(c("red2", "darkblue"))
      color.palette <- temp.col(nrow(out.mod))
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = c("golden rod", "darkgreen"), mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
    }
  }
  
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[4], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image <- as.raster(matrix(color.palette, ncol=1))
  plot(c(0,4),c(0,1),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=3.65, y = seq(0.37,0.52,l=2), labels = seq(min(beta.outvar.out$con.var), max(beta.outvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image, 2.65, 0.34, 2.95, 0.55)
  text(x=3.4, y=0.66, labels = "Outcome", cex=1.6)
  
  legend("left", legend = levels(beta.treatvar.out$bin.var), title = "Treatment",
         fil = c("golden rod", "darkgreen", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)

}

MedTest.con.con.plot1 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Jaccard, U.UniFrac, W.UniFrac
    if (i %% 2 == 1) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <- cmdscale(beta.treatvar.out$Ds[[i]])
      out.mod <- cmdscale(beta.outvar.out$Ds[[i]])
      
      temp.col1 <- colorRampPalette(c("yellow2", "darkgreen"))
      color.palette1 <- temp.col1(nrow(treat.mod))
      temp.col2 <- colorRampPalette(c("red2", "darkblue"))
      color.palette2 <- temp.col2(nrow(out.mod))
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette1, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette2, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
       
    }
  }
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[1], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.35)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[3], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[5], side = 3, line = -50.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[3], side = 3, line = -72.5, outer = TRUE, cex = 1.05)
  
}

MedTest.con.con.plot2 = function(med.test.out, beta.treatvar.out, beta.outvar.out) {
  
  par(mfrow = c(3, 2))
  par(oma = c(5, 0, 5, 0))
  par(mar = c(7, 4, 4, 1))
  
  sub.title = c()
  
  for (i in 1:length(beta.treatvar.out$Ds)) {
    
    # Bray.Curtis, G.UniFrac
    if (i %% 2 == 0) {
      
      if (med.test.out$margPs[i] < 0.05) {
        sub.title <- c(sub.title, paste("*p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      if (med.test.out$margPs[i] >= 0.05) {
        sub.title <- c(sub.title, paste("p: ", P_value(med.test.out$margPs[i]), sep=""))
      }
      
      treat.mod <<- cmdscale(beta.treatvar.out$Ds[[i]])
      out.mod <<- cmdscale(beta.outvar.out$Ds[[i]])
      
      temp.col1 <<- colorRampPalette(c("yellow2", "darkgreen"))
      color.palette1 <<- temp.col1(nrow(treat.mod))
      temp.col2 <<- colorRampPalette(c("red2", "darkblue"))
      color.palette2 <<- temp.col2(nrow(out.mod))
      
      plot(treat.mod, ellipse = TRUE, hull = FALSE, 
           main = "Treatment", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette1, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.5, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
      plot(out.mod, ellipse = TRUE, hull = FALSE, 
           main = "Outcome", xlab="PC 1", ylab="PC 2", sub="", 
           col = color.palette2, mgp=c(2.5,1,0), pch=c(16,17),
           cex=1.4, cex.main=1.4, cex.sub=1.4, label.cex=1.4, cex.lab=1.4)
      
    }
  }
  print(sub.title)
  
  mtext(names(beta.treatvar.out$Ds)[2], side = 3, line = -0.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[1], side = 3, line = -22, outer = TRUE, cex = 1.05)
  mtext(names(beta.treatvar.out$Ds)[4], side = 3, line = -25.5, outer = TRUE, font = 2, cex = 1.25)
  mtext(sub.title[2], side = 3, line = -47, outer = TRUE, cex = 1.05)
  
  legend_image1 <- as.raster(matrix(color.palette1, ncol=1))
  plot(c(0,3),c(0,0.5),type = 'n', axes = F, xlab = '', ylab = '', main = '')
  text(x=0.8, y = seq(0.37,0.52,l=2), labels = seq(min(beta.treatvar.out$con.var), max(beta.treatvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image1, 0.1, 0.34, 0.4, 0.55)
  text(x=0.6, y=0.66, labels = "Treatment", cex=1.6)
  
  legend_image2 <<- as.raster(matrix(color.palette2, ncol=1))
  text(x=3.65, y = seq(0.37,0.52,l=2), labels = seq(min(beta.outvar.out$con.var), max(beta.outvar.out$con.var), length.out=2), cex=1.6)
  rasterImage(legend_image2, 2.65, 0.34, 2.95, 0.55)
  text(x=3.4, y=0.66, labels = "Outcome", cex=1.6)
  
  legend("bottom", paste("Omnibus MedTest: ", P_value(med.test.out$permP), sep=""), bty = "n", cex=1.8)
  
}

###################
# Other functions #
###################

P_value = function(x) {
  
  x <- ifelse(abs(round(x, 3)) == 0, "<.001", 
              ifelse(abs(round(x, 3)) == 1, ">.999", sprintf("%.3f", x)))
  
  return(x)
}
