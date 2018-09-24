# Inputs:
## spike.counts: list of 5 objects: Acounts, Bcounts, ABcounts, bin.mids and bin.width. Acounts is a matrix of spike counts, each column is a single A trial, each row is a bin. Bcounts and ABcounts are defined similarly. bin.mids gives mid points of the bins used for spike counting. bin.width is a scalar giving the width of the bins. Time is measured in ms.
## lengthScale: a vector giving the length-scale values to consider for the GP prior on the eta = log(alpha / (1 - alpha)) curves. Measured in ms. Larger length-scale produces smoother/flatter curves, whereas smaller length-scales produce wavy curves.
## lsPrior: a vector of positive number of the same length as lengthScale. Gives the prior precision (= sum(lsPrior)) and prior expectation (= lsPrior / sum(lsPrior)) for ls.probs which is the truplet's unknown vector of weights on the length-scale to be inferred from data.
## (prec, w1, w2): the overall value of eta is modeled as a mixture of Gaussians, whose behavior is governed by (perc, w1, w2). Higher prec gives larger number of relatively important mixture components, smaller prec gives one or two dominant components.

dapp <- function(spike.counts, lengthScale = NULL, lsPrior = NULL, hyper = list(prec = c(1,1), sig0 = 1.87), burnIn = 1e3, nsamp = 1e3, thin = 4, plot = FALSE, verbose = TRUE, remove.zeros = FALSE){
    
    if(is.null(hyper$prec)) hyper$prec <- c(1,1)
    if(is.null(hyper$sig0)) hyper$sig0 <- 1.87
    
    x1 <- spike.counts$Acounts
    x2 <- spike.counts$Bcounts
    x3 <- spike.counts$ABcounts
    bin.mids <- spike.counts$bin.mids
    bin.width <- spike.counts$bin.width
    
    resp.horiz <- length(bin.mids) * bin.width
    if(is.null(lengthScale)){
        lengthScale <- sort(0.16 * resp.horiz / c(4, 3, 2, 1, 0.5, 0.1))
        lsPrior <- 1:length(lengthScale) ##rep(1, length(lengthScale))
        lsPrior <- 2 * lsPrior / sum(lsPrior)
    }
    if(is.null(lsPrior)) lsPrior <- rep(1, length(lengthScale))/length(lengthScale)
    

    if(remove.zeros){
      x1 <- x1[,colSums(x1) > 0, drop = FALSE]
      x2 <- x2[,colSums(x2) > 0, drop = FALSE]
      x3 <- x3[,colSums(x3) > 0, drop = FALSE]
    }
        
    sig0 <- hyper$sig0
    sigSq0 <- sig0^2
    prec.a <- hyper$prec[1]
    prec.b <- hyper$prec[2]
    prec <- rtruncgamma(1, shape = prec.a, rate = prec.b, low = 1e-8)
    w1 <- 1; w2 <- 0e-4 + prec
    
    nbins <- length(bin.mids)
    if(nrow(x3) != nbins) stop("dimension mismatch between spike counts and bins")
    
    n1 <- ncol(x1)
    n2 <- ncol(x2)
    n3 <- ncol(x3)
    
    #get.hyper1 <- smoogam(x1);
    #get.hyper2 <- smoogam(x2);
    #m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, nbins); s.1 <- sqrt(am.1)/bm.1
    #m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, nbins); s.2 <- sqrt(am.2)/bm.2
    
    x1.smu <- sapply(1:n1, function(i) supsmu(bin.mids, x1[,i])$y)
    x2.smu <- sapply(1:n2, function(i) supsmu(bin.mids, x2[,i])$y)
    x3.smu <- sapply(1:n3, function(i) supsmu(bin.mids, x3[,i])$y)
    
    x.max <- max(max(x1.smu), max(x2.smu), max(x3.smu))
    
    m.1 <- apply(x1.smu, 1, mean); s.1 <- apply(x1.smu, 1, sd)/sqrt(n1)
    am.1 <- n1 * m.1;  bm.1 <- rep(n1, nbins); s.1 <- sqrt(am.1)/bm.1
    m.2 <- apply(x2.smu, 1, mean); s.2 <- apply(x2.smu, 1, sd)/sqrt(n2);
    am.2 <- n2 * m.2; bm.2 <- rep(n2, nbins); s.2 <- sqrt(am.2)/bm.2
    
    rA <- am.1; sA <- bm.1
    rB <- am.2; sB <- bm.2
    
    lambda.A <- rgamma(nbins, shape = rA, rate = sA)
    lambda.B <- rgamma(nbins, shape = rB, rate = sB)
    
    L <- length(lengthScale) ## number of distinct length scale values

    # covariance matrix of pinned GP eta
    # calculate cov matrix, its Cholesky factor and inverse for all length scales
    K.mat = K.inv = K.chol = u.pin = Kinv.u = u.Kinv.u = list()
    dist.sq <- as.matrix(dist(bin.mids))^2
    
    for(l in 1:L){
        K.full <- exp(-0.5 * dist.sq/(lengthScale[l]^2))
        u.pin[[l]] <- rep(sig0, nbins)
        K.pin <- K.full ## no pinning downs here
        
        K.mat[[l]] <- sigSq0 * (K.pin + diag(1e-10, nbins))
        K.chol[[l]] <- chol(K.mat[[l]])
        K.inv[[l]] <- chol2inv(K.chol[[l]])
        Kinv.u[[l]] <- K.inv[[l]] %*% u.pin[[l]]
        u.Kinv.u[[l]] <- sum(c(backsolve(K.chol[[l]], u.pin[[l]], transpose = TRUE))^2)
    }

    ##ls.index <- sample(1:L, n3, replace = TRUE, prob = lsPrior)
    ls.index <- rep(L, n3)
    ##ls.clusters <- lapply(1:L, function(l) which(ls.index == l))
    ##ls.clust.sizes <- sapply(ls.clusters, length)
    ##ls.probs <- rDirichlet(1, lsPrior + ls.clust.sizes)
    
    x3.ave <- colMeans(x3)
    x3.alpha <- pmin(0.99, pmax(0.01, (x3.ave - mean(m.2)) / (0.01 + mean(m.1 - m.2))))
    if(plot) hist(x3.alpha, breaks = seq(0,1,.1), freq = FALSE, col = tcol(3, .2), border = "white", bty = "n", xlab = "alpha.bar", ylab = "density", main = "")
    eta <- outer(rep(1, nbins), qlogis(x3.alpha))
    
    j.eta.Kinv.eta <- rep(NA, n3)
    j.eta.Kinv.u <- rep(NA, n3)
    j.u.Kinv.u <- rep(NA, n3)
    for(j in 1:n3){
        j.ls <- ls.index[j]
        j.eta.Kinv.eta[j] <- c(crossprod(eta[,j], K.inv[[j.ls]] %*% eta[,j]))
        j.eta.Kinv.u[j] <- sum(eta[,j] * Kinv.u[[j.ls]])
        j.u.Kinv.u[j] <- u.Kinv.u[[j.ls]]
    }

    gamma <- 1:n3
    eta.clust <- split(1:n3, gamma)
    nclust <- length(eta.clust)
    eta.clust.size <- c(sapply(eta.clust, length), rep(0, n3 + 1 - nclust))
    occupied.clusts <- which(eta.clust.size > 0)
    clust.sum.eta.Kinv.eta <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.eta[eta.clust[[g]]]))
    clust.sum.eta.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.u[eta.clust[[g]]]))
    clust.sum.u.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.u.Kinv.u[eta.clust[[g]]]))
    
    b.psi <- 1e-2
    sdFn <- function(psi, U) return(sqrt(psi*(1-psi)/(psi + U*(1-psi))))
    meanFn <- function(psi, U, UZ) return((1-psi)*UZ/(psi + U*(1-psi)))
    lgFn <- function(psi, U, Z) return(dnorm(Z, 0, sqrt(psi/U + 1 - psi), log = TRUE) + dbeta(psi, w1 + 1, w2, log = TRUE) + b.psi * psi)
    myTable <- function(x, x.targ) return(sapply(x.targ, function(val) sum(x == val)))

    eta.psi <- rep(NA, n3+1); eta.phi <- rep(NA, n3+1); eta.pi <- matrix(NA, L, n3+1)
    eta.psi[occupied.clusts] <- rtruncbeta(nclust, w1, w2, low = 1e-8, up = 1 - 1e-8)
    eta.phi[occupied.clusts] <- rnorm(nclust, meanFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u, clust.sum.eta.Kinv.u), sdFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u))
    eta.pi[,occupied.clusts] <- sapply(occupied.clusts, function(this.clust) rDirichlet(1, lsPrior + myTable(ls.index[eta.clust[[this.clust]]], 1:L)))

    # empty objects to store select parameter draws
    LSPROB_POST <- matrix(nrow = nsamp, ncol = L)
    dimnames(LSPROB_POST)[[2]] <- lengthScale
    LAMBDA.A_POST <- matrix(nrow = nsamp, ncol = nbins)
    LAMBDA.B_POST <- matrix(nrow = nsamp, ncol = nbins)
    ALPHA <- array(dim = c(nsamp, nbins, n3))
    A_POST <- array(dim = c(nsamp, nbins, n3))
    ALPHA_PRED <- matrix(nrow = nsamp, ncol = nbins)
    PSL_PRED <- matrix(nrow = nsamp, ncol = 3)
    PREC <- rep(NA, nsamp)
    
    
    
    ######## MCMC Loop Begins #########
    
    alpha <- plogis(eta) ## dim(alpha) = c(nbins, n3)
    X <- x3 ## to be consistent with paper. dim(X) = dim(alpha) = c(nbins, n3)
    nTot <- n3*nbins
    
    nextra <- 20
    
    samp.count <- 0
    nacpt <- 0
    niter <- 0
    
    ticker <- (thin*nsamp + burnIn) / 10
    ntick <- 0

    for(iter in -burnIn:(nsamp*thin)){
        
        niter <- niter + 1
        
        # Block 1: (A, A_star, B, B_star | --)
        rate.A <- alpha * lambda.A / (alpha * lambda.A + (1 - alpha) * lambda.B)
        A <- matrix(rbinom(nTot, X, rate.A), nrow = nbins)
        B <- X - A
        A_star <- A + matrix(rpois(nTot, lambda.A * (1 - alpha)), nrow = nbins)
        B_star <- B + matrix(rpois(nTot, lambda.B * alpha), nrow = nbins)
        
        # Block 2: (lambda.A, lambda.B | --)
        lambda.A <- rgamma(nbins, rA + rowSums(A_star), sA + n3)
        lambda.B <- rgamma(nbins, rB + rowSums(B_star), sB + n3)
        
        # Block 3: (Omega | --)
        N <- A_star + B_star
        Omega <- matrix(sapply(1:nTot, function(jj) rpg.devroye(num = 1, n = N[jj], z = eta[jj])), nrow = nbins)
        Omega[Omega == 0] <- 1e-12


        # Block 4: (eta | --)
        kappa <- (A + B_star - B) - 0.5 * (A_star + B_star)
        
        for(j in 1:n3){
            Omega.j <- diag(Omega[,j], nbins)
            Omega.j.inv <- diag(1/Omega[,j], nbins)
            kappa.j <- kappa[,j]
            z.j <- kappa.j / Omega[,j]
            j.clust <- gamma[j]
            
            ## update length-scale
            chol_C_plus_Omega.inv <- lapply(1:L, function(l) return(chol(eta.psi[j.clust] * K.mat[[l]] + Omega.j.inv)))
            logP.ls <- sapply(1:L, function(l) return(-sum(log(diag(chol_C_plus_Omega.inv[[l]]))) - 0.5*sum(c(backsolve(chol_C_plus_Omega.inv[[l]], z.j - eta.phi[j.clust]*u.pin[[l]], transpose = TRUE))^2)))
            logP.ls <- logP.ls + log(eta.pi[,j.clust])
            ls.index[j] <- j.ls <- sample(L, 1, prob = exp(logP.ls - max(logP.ls)))
            
            ## update eta[,j]
            C.tilde.inv <- K.inv[[j.ls]] + eta.psi[j.clust] * Omega.j
            R.tilde <- chol(C.tilde.inv)
            m.tilde <- backsolve(R.tilde, backsolve(R.tilde, eta.psi[j.clust]*kappa.j + eta.phi[j.clust] * Kinv.u[[j.ls]], transpose = TRUE))
            eta[,j] <-  m.tilde + sqrt(eta.psi[j.clust]) * c(backsolve(R.tilde, rnorm(nbins)))
            
        }
        alpha <- plogis(eta)
        
        # Block 5: (ls.probs | --)
        occupied.clusts <- which(eta.clust.size > 0)
        nclust <- length(occupied.clusts)
        eta.pi[,occupied.clusts] <- sapply(occupied.clusts, function(this.clust) rDirichlet(1, lsPrior + myTable(ls.index[eta.clust[[this.clust]]], 1:L)))
        ##ls.clusters <- lapply(1:L, function(l) which(ls.index == l))
        ##ls.clust.sizes <- sapply(ls.clusters, length)
        ##ls.probs <- rDirichlet(1, lsPrior + ls.clust.sizes)
        
        
        # Block 6: (gamma, eta.clust, eta.phi, eta.psi, eta.pi | --)
        #cat("prec =", signif(prec, 7), "\n")
        for(j in 1:n3){
            j.ls <- ls.index[j]

            psi.extra <- rtruncbeta(nextra, w1, w2, low = 1e-8, up = 1 - 1e-8)
            if(any((1 - psi.extra) < 1e-8)) cat("j =", j, "w1 =", w1, "w2 =", w2, "1 - psi.extra =", 1 - psi.extra, "\n")
            phi.extra <- rnorm(nextra, 0, sqrt(1 - psi.extra))
            pi.extra <- t(rDirichlet(nextra, lsPrior))
            
            ## remove j from its current cluster gamma[j]
            eta.clust.size[gamma[j]] <- eta.clust.size[gamma[j]] - 1
            eta.clust[[gamma[j]]] <- drop.item(eta.clust[[gamma[j]]], j)
            if(eta.clust.size[gamma[j]] == 0){
                phi.extra[1] <- eta.phi[gamma[j]]
                psi.extra[1] <- eta.psi[gamma[j]]
                pi.extra[,1] <- eta.pi[,gamma[j]]
            }
            
            occupied.clusts.j <- which(eta.clust.size > 0)
            nclust.j <- length(occupied.clusts.j)
            first.free.clust <- min(which(eta.clust.size == 0))
            phi.all <- c(eta.phi[occupied.clusts.j], phi.extra)
            psi.all <- c(eta.psi[occupied.clusts.j], psi.extra)
            pi.all <- c(eta.pi[j.ls, occupied.clusts.j], pi.extra[j.ls,])
            size.all <- c(eta.clust.size[occupied.clusts.j], rep(prec / nextra, nextra))
            
            j.eta.Kinv.eta[j] <- c(crossprod(eta[,j], K.inv[[j.ls]] %*% eta[,j]))
            j.eta.Kinv.u[j] <- sum(eta[,j] * Kinv.u[[j.ls]])
            j.u.Kinv.u[j] <- u.Kinv.u[[j.ls]]
    
            log.prob <- log(size.all) - 0.5*nbins*log(psi.all) - 0.5 * (j.eta.Kinv.eta[j] - 2*phi.all*j.eta.Kinv.u[j] + phi.all^2*j.u.Kinv.u[j])/ psi.all + log(pi.all)
            pick <- sample(length(log.prob), 1, prob = exp(log.prob - max(log.prob)))
            
            if(pick > nclust.j){ ## new cluster
                gamma[j] <- first.free.clust
                eta.clust.size[gamma[j]] <- 1
                eta.clust[[gamma[j]]] <- j
                eta.psi[gamma[j]] <- psi.extra[pick - nclust.j]
                eta.phi[gamma[j]] <- phi.extra[pick - nclust.j]
                eta.pi[,gamma[j]] <- pi.extra[,pick - nclust.j]
            } else { ## add to existing cluster
                gamma[j] <- occupied.clusts.j[pick]
                eta.clust.size[gamma[j]] <- eta.clust.size[gamma[j]] + 1
                eta.clust[[gamma[j]]] <- c(eta.clust[[gamma[j]]], j)
            }
        }
        
        ## Recolate clusters
        occupied.clusts <- which(eta.clust.size > 0)
        nclust <- length(occupied.clusts)
        
        ## update pi
        eta.pi[,occupied.clusts] <- sapply(occupied.clusts, function(this.clust) rDirichlet(1, lsPrior + myTable(ls.index[eta.clust[[this.clust]]], 1:L)))

        ## update psi for clusters with size > 1
        clust.sum.eta.Kinv.eta <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.eta[eta.clust[[g]]]))
        clust.sum.eta.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.u[eta.clust[[g]]]))
        clust.sum.u.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.u.Kinv.u[eta.clust[[g]]]))
        
        high.occu <- which(eta.clust.size[occupied.clusts] > 1)
        #cat("high clust sizes = ", eta.clust.size[occupied.clusts[high.occu]], "\n")
        n.high <- length(high.occu)
        if(n.high > 0){
            clust.U <- clust.sum.u.Kinv.u[high.occu]
            clust.Z <- clust.sum.eta.Kinv.u[high.occu] / clust.U
            clust.V <- pmax(0, clust.sum.eta.Kinv.eta[high.occu] - clust.U * clust.Z^2)
            shape.new <- (eta.clust.size[occupied.clusts[high.occu]]*nbins - 1)/2
            rate.new <- clust.V/2 + b.psi
            ##cat("Mayday: n.high", n.high, "shape.new = ", signif(shape.new), "rate.new = ", signif(rate.new, 4), "psi.new = ")
            psi.new <- 1/(psi.new.inv <- rtruncgamma(n.high, low = 1 + 1e-6, up = 1e6, shape = shape.new, rate = rate.new))
            ##cat(signif(psi.new, 4), "\n")
            if(any(is.na(psi.new))){
                cat("psi[c(", paste(which(is.na(psi.new)), ","),")] are NA\n")
                cat("shape: "); print(shape.new)
                cat("rate: "); print(rate.new)
            }
            if(any(1 - psi.new < 1e-6)){
                cat("psi[c(", paste(which(1 - psi.new < 1e-6), ","),")] are > 1 - 1e-6\n")
                cat("shape: "); print(shape.new)
                cat("rate: "); print(rate.new)
            }
            l1 <- lgFn(psi.new, clust.U, clust.Z)
            l2 <- lgFn(eta.psi[occupied.clusts[high.occu]], clust.U, clust.Z)
            acpt.new <- exp(l1 - l2)
            if(any(is.na(acpt.new))){
                cat("acpt.new =", signif(acpt.new, 4), "\n")
                cat("1 - eta.psi[high.occu] =", signif(1 - eta.psi[occupied.clusts[high.occu]], 4), "\n")
                cat("1 - psi.new =", signif(1 - psi.new, 4), "\n")
                cat("(prec, w1, w2) = ", signif(c(prec, w1, w2), 4), "\n")
            }
            to.acpt <- (runif(n.high) < acpt.new)
            if(any(to.acpt)){
                eta.psi[occupied.clusts[high.occu]][to.acpt] <- psi.new[to.acpt]
                nacpt <- nacpt + 1
            }
            #cat("1 - eta.psi[occupied] =", signif(1 - eta.psi[occupied.clusts], 4), "\n")
        }
        ## update phi for all occupied clusters
        eta.phi[occupied.clusts] <- rnorm(nclust, meanFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u, clust.sum.eta.Kinv.u), sdFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u))
        
        ## Block 7: (prec | --)
        prec.an <- prec.a + 2 * nclust
        prec.bn <- prec.b - sum(log1p(-eta.psi[occupied.clusts]))
        prec.w <- rtruncbeta(1, prec, n3, low = 1e-8, up = 1 - 1e-8)
        ##cat("prec.an = ", signif(prec.an, 4), "prec.bn = ", signif(prec.bn, 4), "\n")
        prec <- rtruncgamma(1, shape = prec.an, rate = prec.bn - log(prec.w), low = 1e-8)
        w2 <- 0e-4 + prec
        #if(prec < 1e-8){
        #    cat("prec.an = ", signif(prec.an, 4), "prec.bn = ", signif(prec.bn, 4), "prec.w = ", signif(prec.w, 4), "eta.psi =", signif(eta.psi[occupied.clusts], 4), "\n")
        #    cat("clust size = ", eta.clust.size[occupied.clusts], "\n")
        #}
        
        if((iter %% ticker) == 0){
            ntick <- ntick + 1
            if(verbose) cat("iter =", iter, "prec =", round(prec, 4), " : ", paste("(", eta.clust.size[occupied.clusts], "|", round(eta.phi[occupied.clusts], 2), "|", round(eta.psi[occupied.clusts], 2), ")", sep = "", collapse = " + "), "\n")
            if(plot){
                alpha1.grid <- seq(0.01,0.99,.01)
                dens.alpha1 <- sapply(alpha1.grid, function(z) return(sum(c(eta.clust.size[occupied.clusts], prec) * dnorm(qlogis(z), c(eta.phi[occupied.clusts], 0), sqrt(c(eta.psi[occupied.clusts], 1))) / (z*(1-z))))) / (n3 + prec)
                lines(alpha1.grid, dens.alpha1, ty = "l", col = tcol(4, min(1, ntick/15)))
            }
        }


        ## Storage after burn-in
        if(iter > 0 & (iter %% thin) == 0){
            samp.count <- samp.count + 1
            
            ## Key model parameters
            ##LSPROB_POST[samp.count, ] <- ls.probs
            LAMBDA.A_POST[samp.count, ] <- lambda.A
            LAMBDA.B_POST[samp.count, ] <- lambda.B
            ALPHA[samp.count, , ] <- alpha
            A_POST[samp.count, , ] <- A
            PREC[samp.count] <- prec
            
            ## Posterior predictive
            occupied.clusts <- which(eta.clust.size > 0)
            nclust <- length(occupied.clusts)
            new.clust <- sample(nclust + 1, 1, prob = c(eta.clust.size[occupied.clusts], prec))
            if(new.clust > nclust){
                new.psi <- rtruncbeta(1, w1, w2, low = 1e-8, up = 1 - 1e-8)
                new.phi <- rnorm(1, 0, sqrt(1 - new.psi))
                new.pi <- rDirichlet(1, lsPrior)
            } else {
                new.psi <- eta.psi[occupied.clusts[new.clust]]
                new.phi <- eta.phi[occupied.clusts[new.clust]]
                new.pi <- eta.pi[,occupied.clusts[new.clust]]
            }
            LSPROB_POST[samp.count,] <- new.pi
            new.l <- sample(L, 1, prob = new.pi)
            PSL_PRED[samp.count,] <- c(new.phi, new.psi, lengthScale[new.l])
            new.eta <- new.phi * u.pin[[new.l]] + sqrt(new.psi) * c(crossprod(K.chol[[new.l]], rnorm(nbins)))
            ALPHA_PRED[samp.count, ] <- plogis(new.eta)
            
        }
    }

    ################### End of MCMC loop ########################
    
    OUT <- list(lsprob = LSPROB_POST, lambda.A = LAMBDA.A_POST, lambda.B = LAMBDA.B_POST, alpha = ALPHA, A = A_POST, prec = PREC, alpha.pred = ALPHA_PRED, psl.pred = PSL_PRED, details = c(niter = niter, nsamp = nsamp, burnIn = burnIn, thin = thin, acpt = nacpt/niter), hyper = hyper, lengthScale = lengthScale, lsPrior = lsPrior, bin.mids = bin.mids, bin.width = bin.width, mcmc = c(burnIn = burnIn, thin = thin, nsamp = nsamp))
    class(OUT) <- "dapp"
    return(OUT)
                
}

summary.dapp <- function(object, cut.width = 0.1, tilt.prior = FALSE, mesh.tilt = 0.1, nprior = object$mcmc["nsamp"], ...){
    nsamp <- object$mcmc["nsamp"]
    alpha.post.pred <- object$alpha.pred
    minmax.post.pred <- apply(alpha.post.pred, 1, function(a) diff(range(a)))
    
    sim.prior <- dapp.simulate(length(object$bin.mids) * object$bin.width, object$bin.width, lengthScale = object$lengthScale, lsPrior = object$lsPrio, hyper = object$hyper, nsamp = nprior)
    alpha.prior.pred <- sim.prior$alpha.pred
    minmax.prior.pred <- apply(alpha.prior.pred, 1, function(a) diff(range(a)))

    if(tilt.prior){
        dmm.prior.inv <- 1/(1e-10 + hist(minmax.prior.pred, seq(0,1,mesh.tilt), plot = FALSE)$counts)
        tilt.wt.prior <- dmm.prior.inv[cut(minmax.prior.pred, seq(0,1,mesh.tilt), include = TRUE, labels = FALSE)]
        tilt.wt.prior <- nprior * tilt.wt.prior / sum(tilt.wt.prior) ## preserve total weight at nsamp
        tilt.wt.post <- dmm.prior.inv[cut(minmax.post.pred, seq(0,1,mesh.tilt), include = TRUE, labels = FALSE)]
        tilt.wt.post <- nsamp * tilt.wt.post / sum(tilt.wt.post) ## preserve total weight at nsamp
    } else {
        tilt.wt.prior <- rep(1, nprior)
        tilt.wt.post <- rep(1, nsamp)
    }
    
    breaks <- seq(0,1,cut.width)
    cut.names <- levels(cut(seq(0,1,cut.width/2), breaks, include = TRUE))
    ## Minmax
    h.minmax <- cbind('minmax-prior' = wtd.hist(minmax.prior.pred, breaks, plot = FALSE, weight = tilt.wt.prior)$dens * cut.width, 'minmax-post' = wtd.hist(minmax.post.pred, breaks, plot = FALSE, weight = tilt.wt.post)$dens * cut.width)
    dimnames(h.minmax)[[1]] <- cut.names
    
    ## Average alpha.pred
    h.ave <- cbind('ave-prior' = wtd.hist(rowMeans(alpha.prior.pred), breaks, plot = FALSE, weight = tilt.wt.prior)$dens * cut.width, 'ave-post' = wtd.hist(rowMeans(alpha.post.pred), breaks, plot = FALSE, weight = tilt.wt.post)$dens * cut.width)
    dimnames(h.ave)[[1]] <- cut.names

    sm <- cbind(minmax = round(h.minmax,2), ave = round(h.ave,2))
    attr(sm, "tilted") <- tilt.prior
    return(sm)
    
}

plot.dapp <- function(x, add.prior = TRUE, synth.data = NULL, tilt.prior = FALSE, mesh.tilt = 0.1, nprior = x$mcmc["nsamp"], ncurves = 10, ...){
    
    object <- x
    bin.mids <- object$bin.mids
    bin.width <- object$bin.width
    nbins <- length(bin.mids)
    nsamp <- object$mcmc["nsamp"]
    
    if(add.prior|tilt.prior){
        sim.prior <- dapp.simulate(length(bin.mids)*bin.width, bin.width, lengthScale = object$lengthScale, lsPrior = object$lsPrior, hyper = object$hyper, nsamp = nprior)
        lsprob.prior <- sim.prior$lsprob
        alpha.prior.pred <- sim.prior$alpha.pred
        minmax.prior.pred <- apply(alpha.prior.pred, 1, function(a) diff(range(a)))
    }

    lsprob.post <- object$lsprob
    alpha.post <- object$alpha
    alpha.post.pred <- object$alpha.pred
    minmax.post.pred <- apply(alpha.post.pred, 1, function(a) diff(range(a)))
    
    if(tilt.prior){
        dmm.prior.inv <- 1/(1e-10 + hist(minmax.prior.pred, seq(0,1,mesh.tilt), plot = FALSE)$counts)

        tilt.wt.prior <- dmm.prior.inv[cut(minmax.prior.pred, seq(0,1,mesh.tilt), include = TRUE, labels = FALSE)]
        tilt.wt.prior <- nprior * tilt.wt.prior / sum(tilt.wt.prior) ## preserve total weight at nprior
        tilt.wt.post <- dmm.prior.inv[cut(minmax.post.pred, seq(0,1,mesh.tilt), include = TRUE, labels = FALSE)]
        tilt.wt.post <- nsamp * tilt.wt.post / sum(tilt.wt.post) ## preserve total weight at nsamp
        
        resamp.prior <- sample(nprior, nprior, prob = tilt.wt.prior, replace = TRUE)
        resamp.post <- sample(nsamp, nsamp, prob = tilt.wt.post, replace = TRUE)
        
        lsprob.prior <- lsprob.prior[resamp.prior,]
        alpha.prior.pred <- alpha.prior.pred[resamp.prior,]
        minmax.prior.pred <- minmax.prior.pred[resamp.prior]
        
        lsprob.post <- lsprob.post[resamp.post,]
        alpha.post <- alpha.post[resamp.post, , ]
        alpha.post.pred <- alpha.post.pred[resamp.post,]
        minmax.post.pred <- minmax.post.pred[resamp.post]
    }
    
    par(mfrow = c(2,3))
    
    ## Estimated alpha
    alpha.post.mean <- apply(alpha.post, c(2,3), mean)
    plot(bin.mids, 0*bin.mids, ty = "n", ylim = c(0,1), xlab = "Time (ms)", ylab = "Value", bty = "n")
    for(j in 1:ncol(alpha.post.mean)) lines(bin.mids, alpha.post.mean[,j], col = tcol(j, 0.5), lwd = 2)
    title(main = bquote(paste("Estimated ", alpha(t))), font.main = 1)
    
    
    ## Posterior predictive draws for alpha
    ss <- sample(nsamp, ncurves)
    alpha.post.pred.ord <- alpha.post.pred[order(minmax.post.pred, decreasing = TRUE),]

    ##colors <- terrain.colors(length(minmax.post.pred), alpha = 0.8)
    #colors <- tcol(rep(1, nsamp), 0.6)
    #colors <- rainbow(length(minmax.post.pred), alpha = 0.6, end = 2/3)
    colors <- rainbow(7, alpha = 0.6, end = 2/3)[ceiling(7 * (1 - minmax.post.pred))]
    #colors <- tcol(1:length(minmax.post.pred), 0.5)
    plot(bin.mids, 0*bin.mids, ty = "n", ylim = c(0,1), xlab = "Time (ms)", ylab = "Value", bty = "n")
    #for(i in ss) lines(bin.mids, alpha.post.pred.ord[i,], col = colors[i])
    i.ix <- 1
    for(i in ss){
        lines(bin.mids, alpha.post.pred[i,], col = colors[i], lwd = 1, ty = "o", pch = paste(i.ix), cex = 0.5)
        i.ix <- i.ix + 1
    }
    title(main = bquote(paste(alpha(t), ": posterior predictive draws")), font.main = 1)
    
    ## lambda.A and lambda.B estimates
    
    count2rate.factor <- 1000 / bin.width
    m.1 <- count2rate.factor * colMeans(object$lambda.A); s.1 <- count2rate.factor * apply(object$lambda.A, 2, sd)
    m.2 <- count2rate.factor * colMeans(object$lambda.B); s.2 <- count2rate.factor * apply(object$lambda.B, 2, sd)
    
    plot(bin.mids, 0*bin.mids, ylim = c(0, max(max(m.1+2*s.1), max(m.2+2*s.2))), ann = FALSE, ty = "n", bty = "n")
    polygon(bin.mids[c(1:nbins, nbins:1)], c(m.1 - 2*s.1, (m.1 + 2*s.1)[nbins:1]), col = tcol("orange", .5), border = tcol("orange", .5))
    polygon(bin.mids[c(1:nbins, nbins:1)], c(m.2 - 2*s.2, (m.2 + 2*s.2)[nbins:1]), col = tcol("cyan", .5), border = tcol("cyan", .5))
    abline(v = 0, lty = 2)
    title(xlab = "Time (ms)", ylab = "Firing rate (Hz)")
    title(main = "Single stimulus response", font.main = 1)
    
    ## LS-Prob (to be removed by Azeem)
    #boxplot(lsprob.post, outline = FALSE, col = tcol(2, .3), border = tcol(2, .5), xlab = "length-scale (ms)", ylab = "probability", ylim = c(0,1))
    #if(add.prior) boxplot(lsprob.prior, outline = FALSE, col = tcol(4, .1), border = tcol(4, .5), add = TRUE)
    #title(main = "Length-scale distn")

    ## Minmax
    h1 <- hist(minmax.post.pred, breaks = seq(0,1,.1), plot = FALSE)
    h2 <- hist(minmax.prior.pred, breaks = seq(0,1,.1), plot = FALSE)
    plot(h1, freq = FALSE, col = tcol(2, .3), border = tcol(2, .3), bty = "n", ann = FALSE, ylim = c(0, max(max(h1$density), max(h2$density))))
    if(add.prior) plot(h2, freq = FALSE, col = tcol(4, .1), border = tcol(4, .1), add = TRUE)
    title(xlab = bquote(max[t] * alpha(t) - min[t] * alpha(t)), ylab = "Density", main = bquote(paste("Predicted range of ", alpha(t))), font.main = 1)
    
    ## Average alpha.pred
    h1 <- hist(rowMeans(alpha.post.pred), breaks = seq(0,1,.1), plot = FALSE)
    h2 <- hist(rowMeans(alpha.prior.pred), breaks = seq(0,1,.1), plot = FALSE)
    plot(h1, freq = FALSE, col = tcol(2, .3), border = tcol(2, .3), bty = "n", ann = FALSE, ylim = c(0, max(max(h1$density), max(h2$density))))
    if(add.prior) plot(h2, freq = FALSE, col = tcol(4, .1), border = tcol(4, .1), add = TRUE)
    title(xlab = bquote(bar(alpha) == integral(alpha(t)*dt,0,T)/T), ylab = "Density", main = bquote(paste("Predicted long term average of ", alpha(t))), font.main = 1)
    
    ## Swings
    #swing.cuts <- c(0.2, 0.4, 0.6)
    #nswing.prior.pred <- nswing.post.pred <- list()
    #for(i in 1:length(swing.cuts)) nswing.post.pred[[i]] <- apply(alpha.post.pred, 1, swing.counter, width = swing.cuts[i])
    #names(nswing.post.pred) <- swing.cuts
    #boxplot(nswing.post.pred, xlab = "swing magnitude", ylab = "counts", col = tcol(2, .3), border = tcol(2, .5), pch = 19)
    
    #if(add.prior){
    #    for(i in 1:length(swing.cuts)) nswing.prior.pred[[i]] <- apply(alpha.prior.pred, 1, swing.counter, width = swing.cuts[i])
    #    names(nswing.prior.pred) <- swing.cuts
    #    boxplot(nswing.prior.pred, col = tcol(4, .1), border = tcol(4, .5), pch = 19, add = TRUE)
    #}
    #title(main = "Swings")

    ## Length-scale
    lsprobs <- rbind(colMeans(lsprob.post), colMeans(lsprob.prior))
    resp.horiz <- bin.width * nbins
    E.swing <- 0.16 * resp.horiz / object$lengthScale
    swing.ix <- order(E.swing)
    dimnames(lsprobs)[[2]] <- round(E.swing, 1)
    barplot(lsprobs[,swing.ix], beside = TRUE,
        col = c(tcol(2,.3), tcol(4,.1)), border = c(tcol(2,.3), tcol(4,.1)), ylim = c(0,1),
        xlab = bquote(paste(0.16*T, "/\u2113 = E(#up-crossing)")), ylab = "Probability")
    title(main = bquote(paste("Predicted waviness of ", alpha(t))), font.main = 1)
    
    if(!is.null(synth.data)){
        ss <- sort(sample(ncol(alpha.post.mean), 6))
        for(j in ss){
            plot(bin.mids, 0*bin.mids, ty = "n", ylim = c(0,1), xlab = "t (ms)", ylab = bquote(alpha(t)), bty = "n", main = paste("AB trial #", j), font.main = 1)
            lines(bin.mids, alpha.post.mean[,j], col = tcol(1, 0.5), lwd = 4)
            lines(synth.data$time.pts, synth.data$alphas[,j], col = 1, lty = 2)
            #legend("topright", c("true", "estimated"), lty = c(2,1), lwd = c(1,4), col = c(j, tcol(j,.5)), bty = "n")
        }
    }
}

dapp.simulate <- function(horizon = 1000, bin.width = 25, lengthScale, lsPrior = rep(1/length(lengthScale),length(lengthScale)), hyper = list(prec = c(1,1), sig0 = 1.87), nsamp = 1e3){
    
    horizon <- bin.width * (horizon %/% bin.width)
    bin.mids <- seq(bin.width/2, horizon, bin.width)
    sig0 <- hyper$sig0
    sigSq0 <- sig0^2
    prec.a <- hyper$prec[1]
    prec.b <- hyper$prec[2]
    
    nbins <- length(bin.mids)
    if(missing(lengthScale)) lengthScale <- sort(0.16 * horizon / c(4, 3, 2, 1, 0.5, 0.1))
    L <- length(lengthScale) ## number of distinct length scale values
    
    # covariance matrix of pinned GP eta
    # calculate cov matrix, its Cholesky factor and inverse for all length scales
    K.mat = K.inv = K.chol = u.pin = Kinv.u = u.Kinv.u = list()
    dist.sq <- as.matrix(dist(bin.mids))^2
    
    for(l in 1:L){
        K.full <- exp(-0.5 * dist.sq/(lengthScale[l]^2))
        #u.pin[[l]] <- K.full[,1] / K.full[1,1]
        #K.pin <- K.full - tcrossprod(K.full[,1]) / K.full[1,1]
        u.pin[[l]] <- rep(sig0, nbins)
        K.pin <- K.full ## no pinning downs here
        
        K.mat[[l]] <- sigSq0 * (K.pin + diag(1e-10, nbins))
        K.chol[[l]] <- chol(K.mat[[l]])
        K.inv[[l]] <- chol2inv(K.chol[[l]])
        Kinv.u[[l]] <- K.inv[[l]] %*% u.pin[[l]]
        u.Kinv.u[[l]] <- sum(c(backsolve(K.chol[[l]], u.pin[[l]], transpose = TRUE))^2)
    }


    ALPHA_PRED <- matrix(nrow = nsamp, ncol = nbins)
    LSPROB <- matrix(nrow = nsamp, ncol = L)
    dimnames(LSPROB)[[2]] <- lengthScale
    PREC <- rep(NA, nsamp)

    for(samp.count in 1:nsamp){
        PREC[samp.count] <- prec <- rtruncgamma(1, shape = prec.a, rate = prec.b, low = 1e-8)
        w1 <- 1; w2 <- 0e-4 + prec
        LSPROB[samp.count,] <- ls.probs <- rDirichlet(1, lsPrior)
        psi <- rtruncbeta(1, w1, w2, low = 1e-8, up = 1 - 1e-8)
        phi <- rnorm(1, 0, sqrt(1 - psi))
        l <- sample(L, 1, prob = ls.probs)
        eta <- phi * u.pin[[l]] + sqrt(psi) * c(crossprod(K.chol[[l]], rnorm(nbins)))
        ALPHA_PRED[samp.count, ] <- plogis(eta)
    }

    return(list(lsprob = LSPROB, alpha.pred = ALPHA_PRED, prec = PREC))
}

rPoiProc <- function(lambda, time.bins = 0:1000){
    ## time measured in milliseconds. lambda is a vector of values of
    ## a rate function (in Hz) recorded at the mid-points of time.bins
    
    bin.widths <- diff(time.bins)
    nbins <- length(bin.widths)
    
    nspikes.bins <- rpois(nbins, bin.widths * lambda / 1e3)
    bins.with.spikes <- rep(1:nbins, nspikes.bins)
    return(runif(sum(nspikes.bins), time.bins[bins.with.spikes], time.bins[1+bins.with.spikes]))
}

flatFn <- function(intervals = list(c(0,1)), wts = 1, time.pts = 1:1e3 - .5){
    k <- sample(length(wts), 1, prob = wts)
    return(rep(runif(1, intervals[[k]][1], intervals[[k]][2]), length(time.pts)))
}

sineFn <- function(span = c(0, 1), period.range = c(400, 1000), time.pts = 1:1e3 - 0.5){
    period <- runif(1, period.range[1], period.range[2])
    jit <- runif(1, 0, period)
    return(span[1] + diff(span) * (1 + sin(2 * pi * (jit + time.pts) / period)) / 2)
}

## synthetic triplet data (A, B, AB) generation
###### INPUT #####
## ntrials: 3-vector giving nubers of trials for each condition
## time.bins: a fine grid of time points (in ms) to measure spike times
## lambda.A: flat average intensity of A trials
## lambda.B: flat average intensity of B trials
## pr.flat: probability with which an AB trial is flat
## intervals, wts: flat AB trial intensities are generated according to a mixture of uniform distributions on the segments specified by the list 'intervals' according to the weights given by wts. For example, one could specify intervals = list(c(0, 0.25)), wts = 1 to get flat AB trials that are close to B. One could also do intervals = list(c(0, 0.25), c(0.75, 1)), wts = c(2/3, 1/3) will get flat intensities for AB trials uniformly from the two intervals in a 2:1 ratio
## span, period.range: Non-flat AB trials are generated from according to sine curves restricted to the range specified by 'span' and with periods uniformly generated from period.range. These curves are given random phase shifts uniformly generated from (0, period).
####### OUTPUT #####
## returns a list containing
## spiketimes: a list containing 3 sublists, one for each of A, B and AB conditions. Each sublist containts a further sublist of vectors of spike times
## alphas: alpha curves for the AB trials
## lambdas: the corresponding lambda.AB intensity functions
## time.pts: the time points (mid points of the bins in time.bins) at which alphas and lambdas are evaluated

synthesis.dapp <- function(ntrials = c(10, 10, 10), time.bins = 0:1000, lambda.A = 400, lambda.B = 100, pr.flat = 0.5, intervals = list(c(0,1)), wts = 1, span = c(0,1), period.range = c(400, 1000)){
    spiketimes <- list()
    spiketimes$A <- replicate(ntrials[1], rPoiProc(lambda.A, time.bins))
    spiketimes$B <- replicate(ntrials[2], rPoiProc(lambda.B, time.bins))
    
    nflat <- rbinom(1, size = ntrials[3], prob = pr.flat)
    alphas.flat <- numeric(0); alphas.sine <- numeric(0)
    time.pts <- time.bins[-1] - diff(time.bins)/2
    if(nflat > 0) alphas.flat <- replicate(nflat, flatFn(intervals, wts, time.pts))
    if(nflat < ntrials[3]) alphas.sine <- replicate(ntrials[3] - nflat, sineFn(span, period.range, time.pts))
    alphas <- cbind(alphas.flat, alphas.sine)[,sample(ntrials[3])]
    lambdas <- apply(alphas, 2, function(alpha) return(lambda.A * alpha + lambda.B * (1 - alpha)))
    
    spiketimes$AB <- lapply(1:ntrials[3], function(j) rPoiProc(lambdas[,j], time.bins = time.bins))
    return(list(spiketimes = spiketimes, alphas = alphas, lambdas = lambdas, time.pts = time.pts))
}


## auxiliary functions

## already defined in poisson_analysis.R
## bin.counter <- function(x, b) return(diff(sapply(b, function(a) sum(x <= a))))

## already defined in poisson_analysis.R
## logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

rDirichlet <- function(n, alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- rowSums(x)
    return(x/sm)
}

rinvgamma <-function(n,shape, scale = 1) return(1/rgamma(n = n, shape = shape, rate = scale))

rtruncgamma <- function(n, low = 0, up = Inf, ...){
    unifs <- runif(n)
    if(length(low) < n) low <- rep(low, ceiling(n / length(low)))[1:n]
    if(length(up) < n) up <- rep(up, ceiling(n / length(up)))[1:n]
    if(any(low > up)) stop("lower truncation point is larger than the upper truncation point")
    plo <- pgamma(low, ...)
    pup <- pgamma(up, ...)
    vals <- qgamma(plo + unifs * (pup - plo), ...)
    if(any(is.na(vals) | (vals == Inf) | (vals == 0))){
        pos1 <- pgamma(up, ...) < 0.5
        pos2 <- pgamma(low, ...) > 0.5
        ## those with pos1 = TRUE
        if(any(pos1)){
            lplo <- pgamma(low, ..., log.p = TRUE, lower.tail = TRUE)
            lpup <- pgamma(up, ..., log.p = TRUE, lower.tail = TRUE)
            lp <- apply(cbind(lplo, log(unifs) + lpup + log1p(-exp(lplo - lpup))), 1, logsum)
            vals[pos1] <- qgamma(lp, ..., log.p = TRUE, lower.tail = TRUE)[pos1]
        }
        ## now those with pos2 = TRUE
        if(any(pos2)){
            lplo <- pgamma(low, ..., log.p = TRUE, lower.tail = FALSE)
            lpup <- pgamma(up, ..., log.p = TRUE, lower.tail = FALSE)
            lp <- lplo + log1p(unifs * expm1(lpup - lplo))
            lp <- apply(cbind(lplo, log(unifs) + lplo + log1p(-exp(lpup - lplo))), 1, logsum)
            vals[pos2] <- qgamma(lp, ..., log.p = TRUE, lower.tail = FALSE)[pos2]
        }
    }
    return(vals)
}

rtruncbeta <- function(n, low = 0, up = Inf, ...){
    unifs <- runif(n)
    vals <- rep(NA, n)
    pos <- (pbeta(low, ...) < 0.5)
    ## those with pos = TRUE
    lplo <- pbeta(low, ..., log.p = TRUE, lower.tail = TRUE)
    lpup <- pbeta(up, ..., log.p = TRUE, lower.tail = TRUE)
    lp <- lplo + log1p(unifs * expm1(lpup - lplo))
    vals[pos] <- qbeta(lp, ..., log.p = TRUE, lower.tail = TRUE)[pos]
    ## now those with pos = FALSE
    lplo <- pbeta(low, ..., log.p = TRUE, lower.tail = FALSE)
    lpup <- pbeta(up, ..., log.p = TRUE, lower.tail = FALSE)
    lp <- lplo + log1p(unifs * expm1(lpup - lplo))
    vals[!pos] <- qbeta(lp, ..., log.p = TRUE, lower.tail = FALSE)[!pos]
    return(vals)
}

#smoogam <- function(x){
#    T <- nrow(x)
#    n <- ncol(x)
#    if(T > 1){
#        x.dd <- data.frame(cts = c(x), time = rep(1:T, n))
#        x.gam <- gam(cts ~ s(time, bs = "ad"), data = x.dd, family = poisson(link = "log"))
#        x.pred <- predict(x.gam, data.frame(cts = NA, time = 1:T), se.fit = TRUE)
#        mu <- x.pred$fit; sig <- x.pred$se.fit
#    } else {
#        mu <- log(mean(c(x))); sig <- log(sd(c(x)))
#    }
#    firingrate.mean <- exp(mu)
#    firingrate.vari <- expm1(sig^2/2) * firingrate.mean^2
#    return(list(a = 1 / expm1(sig^2/2), b = 1/(expm1(sig^2/2) * firingrate.mean),
#    mean = firingrate.mean, vari = firingrate.vari))
#}

drop.item <- function(x, j) return(x[-match(j, x, nomatch = length(x) + 1)])

tcol <- Vectorize(function(col, alpha = 1) {x <- col2rgb(col)/255; return(rgb(x[1],x[2],x[3],alpha = alpha))}, "col")


