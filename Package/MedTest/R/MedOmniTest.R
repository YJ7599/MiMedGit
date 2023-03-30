MedOmniTest <-
		function(x, y, m.list, z=NULL, nperm=9999) {
	
	n <- length(y)
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	perm.mat <-  matrix(NA, nrow=nperm, ncol=length(m.list))
	
	obs.vec <- numeric(length(m.list))
	
	if (is.null(z)) {
		x.adj <- resid(lm(x ~ 1))
		z1 <- rep(1, n)
		Px <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
	} else {
		x.adj <- resid(lm(x ~ z))
		z1 <- cbind(1, z)
		Px <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
	}
	
	if (is.null(z)) {
		y.adj <- resid(lm(y ~ x)) 
		z1 <- cbind(1, x)
		Py <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
		
	} else {
		y.adj <- resid(lm(y ~ x + z))
		z1 <- cbind(1, x, z)
		Py <- diag(n) - z1 %*% solve(t(z1) %*% z1) %*% t(z1) 
	}
	for (j in 1:length(m.list)) {
		assign(".Random.seed", save.seed, .GlobalEnv)
		
		dmat <- as.matrix(m.list[[j]])
		n <- nrow(dmat)
		dmat <- -dmat^2/2
		
		G <- mean(dmat) + dmat - rowMeans(dmat) - matrix(rep(1, n), ncol=1) %*% colMeans(dmat)
		
		obj <- eigen(G, symmetric=TRUE)
		lambda <- obj$values
		u <- obj$vectors
		u <- u[, lambda > 1E-8, drop = FALSE]
		lambda <- lambda[lambda > 1E-8]
		if (is.null(z)) {
			u.adj <- apply(u, 2, function (xx) resid(lm(xx ~ 1)))
			u.x.adj <- apply(u, 2, function (xx) resid(lm(xx ~ x)))
		} else {
			u.adj <- apply(u, 2, function (xx) resid(lm(xx ~ z)))
			u.x.adj <- apply(u, 2, function (xx) resid(lm(xx ~ z + x)))
		}
		
		obs.vec[j] <- sum(lambda * abs(t(u.adj) %*% x.adj) * abs(t(u.x.adj) %*% y.adj))
		perm.mat[, j] <- sapply(1:nperm, function(i) {
					x.adj.p <- Px %*% x.adj[sample(n)]
					y.adj.p <- Py %*% y.adj[sample(n)]
					f1.f2.p <- sum(lambda * abs(t(u.adj) %*% x.adj) * abs(t(u.x.adj) %*% y.adj.p))
					f2.f1.p <- sum(lambda * abs(t(u.adj) %*% x.adj.p) * abs(t(u.x.adj) %*% y.adj))
					f1.p.f2.p <- sum(lambda * abs(t(u.adj) %*% x.adj.p) * 
									abs(t(u.x.adj) %*% y.adj.p))
					max(f1.f2.p, f2.f1.p, f1.p.f2.p)
				})
	}
# marginal p values
	margPs <- colMeans(rbind(sweep(perm.mat, 2, obs.vec, '>'), rep(1, length(m.list)))) 
	
	
	perm.mat.p <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
	
	
	list(margPs = margPs, 
#margP = min(min(margPs * length(m.list)), 1),
			permP = mean(c(rowMins(perm.mat.p) <= min(margPs), 1)))
	
}
