# Inputs:
## spike.counts: list of 5 objects: Acounts, Bcounts, ABcounts, bin.mids and bin.width. Acounts is a matrix of spike counts, each column is a single A trial, each row is a bin. Bcounts and ABcounts are defined similarly. bin.mids gives mid points of the bins used for spike counting. bin.width is a scalar giving the width of the bins. Time is measured in ms.
## lengthScale: a vector giving the length-scale values to consider for the GP prior on the eta = log(alpha / (1 - alpha)) curves. Measured in ms. Larger length-scale produces smoother/flatter curves, whereas smaller length-scales produce wavy curves.
## lsPrior: a vector of positive number of the same length as lengthScale. Gives the prior precision (= sum(lsPrior)) and prior expectation (= lsPrior / sum(lsPrior)) for ls.probs which is the truplet's unknown vector of weights on the length-scale to be inferred from data.
## (prec, w1, w2): the overall value of eta is modeled as a mixture of Gaussians, whose behavior is governed by (perc, w1, w2). Higher prec gives larger number of relatively important mixture components, smaller prec gives one or two dominant components.

dynamic.model.fit <- function(spike.counts, lengthScale = c(75, 125, 200, 300, 500), lsPrior = rep(1/length(lengthScale),length(lengthScale)), hyper = list(prec = c(1,1), sig0 = 1.87), burnIn = 1e3, nsamp = 1e3, thin = 4, plot = FALSE, verbose = TRUE){
    
    require(BayesLogit)
    if(is.null(hyper$prec)) hyper$prec <- c(1,1)
    if(is.null(hyper$sig0)) hyper$sig0 <- 1.87
    
    x1 <- spike.counts$Acounts
    x2 <- spike.counts$Bcounts
    x3 <- spike.counts$ABcounts
    bin.mids <- spike.counts$bin.mids
    bin.width <- spike.counts$bin.width
    
    sig0 <- hyper$sig0
    sigSq0 <- sig0^2
    prec.a <- hyper$prec[1]
    prec.b <- hyper$prec[2]
    prec <- rgamma(1, prec.a, prec.b)
    w1 <- prec; w2 <- 1
    
    nbins <- length(bin.mids)
    if(nrow(x3) != nbins) stop("dimension mismatch between spike counts and bins")
    
    n1 <- ncol(x1)
    n2 <- ncol(x2)
    n3 <- ncol(x3)
    
    get.hyper1 <- smoogam(x1);
    get.hyper2 <- smoogam(x2);
    m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, nbins); s.1 <- sqrt(am.1)/bm.1
    m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, nbins); s.2 <- sqrt(am.2)/bm.2
    
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

    ls.index <- rep(L, n3)
    ls.clusters <- lapply(1:L, function(l) which(ls.index == l))
    ls.clust.sizes <- sapply(ls.clusters, length)
    ls.probs <- rDirichlet(1, lsPrior + ls.clust.sizes)
    
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
    lgFn <- function(psi, U, Z) return(dnorm(Z, 0, psi/U + 1 - psi, log = TRUE) + dbeta(psi, w1+1, w2, log = TRUE) + b.psi * psi)

    eta.psi <- rep(NA, n3+1); eta.phi <- rep(NA, n3+1)
    eta.psi[occupied.clusts] <- rbeta(nclust, w1, w2)
    eta.phi[occupied.clusts] <- rnorm(nclust, meanFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u, clust.sum.eta.Kinv.u), sdFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u))

    # empty objects to store select parameter draws
    LSPROB_POST <- matrix(nrow = nsamp, ncol = L)
    dimnames(LSPROB_POST)[[2]] <- lengthScale
    LAMBDA.A_POST <- matrix(nrow = nsamp, ncol = nbins)
    LAMBDA.B_POST <- matrix(nrow = nsamp, ncol = nbins)
    ALPHA <- array(dim = c(nsamp, nbins, n3))
    A_POST <- array(dim = c(nsamp, nbins, n3))
    ALPHA_PRED <- matrix(nrow = nsamp, ncol = nbins)
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
            logP.ls <- logP.ls + log(ls.probs)
            ls.index[j] <- j.ls <- sample(L, 1, prob = exp(logP.ls - max(logP.ls)))
            
            ## update eta[,j]
            C.tilde.inv <- K.inv[[j.ls]] + eta.psi[j.clust] * Omega.j
            R.tilde <- chol(C.tilde.inv)
            m.tilde <- backsolve(R.tilde, backsolve(R.tilde, eta.psi[j.clust]*kappa.j + eta.phi[j.clust] * Kinv.u[[j.ls]], transpose = TRUE))
            eta[,j] <-  m.tilde + sqrt(eta.psi[j.clust]) * c(backsolve(R.tilde, rnorm(nbins)))
            
        }
        alpha <- plogis(eta)
        
        # Block 5: (ls.probs | --)
        ls.clusters <- lapply(1:L, function(l) which(ls.index == l))
        ls.clust.sizes <- sapply(ls.clusters, length)
        ls.probs <- rDirichlet(1, lsPrior + ls.clust.sizes)
        
        
        # Block 6: (gamma, eta.clust, eta.phi, eta.psi | --)
        for(j in 1:n3){
            psi.extra <- rbeta(nextra, w1, w2)
            phi.extra <- rnorm(nextra, 0, sqrt(1 - psi.extra))
            
            ## remove j from its current cluster gamma[j]
            eta.clust.size[gamma[j]] <- eta.clust.size[gamma[j]] - 1
            eta.clust[[gamma[j]]] <- drop.item(eta.clust[[gamma[j]]], j)
            
            occupied.clusts <- which(eta.clust.size > 0)
            nclust <- length(occupied.clusts)
            first.free.clust <- min(which(eta.clust.size == 0))
            phi.all <- c(eta.phi[occupied.clusts], phi.extra)
            psi.all <- c(eta.psi[occupied.clusts], psi.extra)
            size.all <- c(eta.clust.size[occupied.clusts], rep(prec / nextra, nextra))
            
            j.ls <- ls.index[j]
            j.eta.Kinv.eta[j] <- c(crossprod(eta[,j], K.inv[[j.ls]] %*% eta[,j]))
            j.eta.Kinv.u[j] <- sum(eta[,j] * Kinv.u[[j.ls]])
            j.u.Kinv.u[j] <- u.Kinv.u[[j.ls]]
    
            log.prob <- log(size.all) - 0.5*nbins*log(psi.all) - 0.5 * (j.eta.Kinv.eta[j] - 2*phi.all*j.eta.Kinv.u[j] + phi.all^2*j.u.Kinv.u[j])/ psi.all
            pick <- sample(length(log.prob), 1, prob = exp(log.prob - max(log.prob)))
            
            if(pick > nclust){ ## new cluster
                gamma[j] <- first.free.clust
                eta.clust.size[gamma[j]] <- 1
                eta.clust[[gamma[j]]] <- j
                eta.psi[gamma[j]] <- psi.extra[pick - nclust]
                eta.phi[gamma[j]] <- phi.extra[pick - nclust]
                
                #psi.extra[pick - nclust] <- rbeta(1, w1, w2)
                #phi.extra[pick - nclust] <- rnorm(1, 0, sqrt(1 - psi.extra[pick - nclust]))
            } else { ## add to existing cluster
                gamma[j] <- occupied.clusts[pick]
                eta.clust.size[gamma[j]] <- eta.clust.size[gamma[j]] + 1
                eta.clust[[gamma[j]]] <- c(eta.clust[[gamma[j]]], j)
            }
        }
        
        ## update psi for clusters with size > 1
        occupied.clusts <- which(eta.clust.size > 0)
        nclust <- length(occupied.clusts)
        clust.sum.eta.Kinv.eta <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.eta[eta.clust[[g]]]))
        clust.sum.eta.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.eta.Kinv.u[eta.clust[[g]]]))
        clust.sum.u.Kinv.u <- sapply(occupied.clusts, function(g) sum(j.u.Kinv.u[eta.clust[[g]]]))
        
        high.occu <- which(nbins * eta.clust.size[occupied.clusts] > 1)
        n.high <- length(high.occu)
        if(n.high > 0){
            clust.U <- clust.sum.u.Kinv.u[high.occu]
            clust.Z <- clust.sum.eta.Kinv.u[high.occu] / clust.U
            clust.V <- pmax(0, clust.sum.eta.Kinv.eta[high.occu] - clust.U * clust.Z^2)
            shape.new <- (eta.clust.size[occupied.clusts[high.occu]]*nbins - 1)/2
            rate.new <- clust.V/2 + b.psi
            psi.new <- 1/(psi.new.inv <- rtruncgamma(n.high, low = 1 + 1e-6, up = 1e10, shape = shape.new, rate = rate.new))
            if(any(is.na(psi.new))){
                print(shape.new)
                print(rate.new)
            }
            l1 <- lgFn(psi.new, clust.U, clust.Z)
            l2 <- lgFn(eta.psi[occupied.clusts[high.occu]], clust.U, clust.Z)
            acpt.new <- exp(l1 - l2)
            to.acpt <- (runif(n.high) < acpt.new)
            if(any(to.acpt)){
                eta.psi[occupied.clusts[high.occu]][to.acpt] <- psi.new[to.acpt]
                nacpt <- nacpt + 1
            }
            
        }
        ## update phi for all occupied clusters
        eta.phi[occupied.clusts] <- rnorm(nclust, meanFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u, clust.sum.eta.Kinv.u), sdFn(eta.psi[occupied.clusts], clust.sum.u.Kinv.u))
        
        ## Block 7: (prec | --)
        prec.p <- -sum(log(eta.psi[occupied.clusts]))
        prec.aux <- rbeta(1, prec + 1, n3)
        prec.pi.odds <- (prec.a + 2 * nclust - 1)/ (n3 * (prec.b + prec.p - log(prec.aux)))
        prec.2nd <- (runif(1) < 1/(1 + prec.pi.odds))
        prec <- rgamma(1, prec.a + 2*nclust - prec.2nd, prec.b + prec.p - log(prec.aux))
        w1 <- prec
        
        if((iter %% ticker) == 0){
            ntick <- ntick + 1
            if(verbose) cat("iter =", iter, " : ", paste("(", eta.clust.size[occupied.clusts], "|", round(eta.phi[occupied.clusts], 2), "|", round(eta.psi[occupied.clusts], 2), ")", sep = "", collapse = " + "), "\n")
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
            LSPROB_POST[samp.count, ] <- ls.probs
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
                new.psi <- rbeta(1, w1, w2)
                new.phi <- rnorm(1, 0, sqrt(1 - new.psi))
            } else {
                new.psi <- eta.psi[occupied.clusts[new.clust]]
                new.phi <- eta.phi[occupied.clusts[new.clust]]
            }
            
            new.l <- sample(L, 1, prob = ls.probs)
            new.eta <- new.phi * u.pin[[new.l]] + sqrt(new.psi) * c(crossprod(K.chol[[new.l]], rnorm(nbins)))
            ALPHA_PRED[samp.count, ] <- plogis(new.eta)
            
        }
    }

    ################### End of MCMC loop ########################
    
    OUT <- list(lsprob = LSPROB_POST, lambda.A = LAMBDA.A_POST, lambda.B = LAMBDA.B_POST, alpha = ALPHA, A = A_POST, prec = PREC, alpha.pred = ALPHA_PRED, details = c(niter = niter, nsamp = nsamp, burnIn = burnIn, thin = thin, acpt = nacpt/niter), hyper = hyper, lengthScale = lengthScale, lsPrior = lsPrior, bin.mids = bin.mids, bin.width = bin.width, mcmc = c(burnIn = burnIn, thin = thin, nsamp = nsamp))
    class(OUT) <- "neurodyn"
    return(OUT)
                
}

summary.neurodyn <- function(fit, cut.width = 0.1, tilt.prior = FALSE, mesh.tilt = 0.1, nprior = fit$mcmc["nsamp"]){
    nsamp <- fit$mcmc["nsamp"]
    alpha.post.pred <- fit$alpha.pred
    minmax.post.pred <- apply(alpha.post.pred, 1, function(a) diff(range(a)))
    
    sim.prior <- dynamic.model.simu(fit$bin.mids, fit$bin.width, lengthScale = fit$lengthScale, lsPrior = fit$lsPrio, hyper = fit$hyper, nsamp = nprior)
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
    
    require(weights)
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

plot.neurodyn <- function(fit, add.prior = TRUE, true.alphas = NULL, tilt.prior = FALSE, mesh.tilt = 0.1, nprior = fit$mcmc["nsamp"]){
    
    bin.mids <- fit$bin.mids
    bin.width <- fit$bin.width
    nbins <- length(bin.mids)
    nsamp <- fit$mcmc["nsamp"]
    
    if(add.prior|tilt.prior){
        sim.prior <- dynamic.model.simu(bin.mids, bin.width, lengthScale = fit$lengthScale, lsPrior = fit$lsPrior, hyper = fit$hyper, nsamp = nprior)
        lsprob.prior <- sim.prior$lsprob
        alpha.prior.pred <- sim.prior$alpha.pred
        minmax.prior.pred <- apply(alpha.prior.pred, 1, function(a) diff(range(a)))
    }

    lsprob.post <- fit$lsprob
    alpha.post <- fit$alpha
    alpha.post.pred <- fit$alpha.pred
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
    
    par(mfrow = c(3,2))
    
    ## Estimated alpha
    alpha.post.mean <- apply(alpha.post, c(2,3), mean)
    plot(bin.mids, 0*bin.mids, ty = "n", ylim = c(0,1), xlab = "time (ms)", ylab = "alpha", bty = "n")
    for(j in 1:ncol(alpha.post.mean)) lines(bin.mids, alpha.post.mean[,j], col = tcol(j, 0.5), lwd = 4)
    if(length(true.alphas) > 0){
        for(j in 1:ncol(alpha.post.mean)) lines(synth.data$time.pts, synth.data$alphas[,j], col = j, lty = 2)
        title(main = "Alpha (dashed = true, solid = est)")
    } else {
        title(main = "Alpha estimates")
    }
    
    ## Posterior predictive draws for alpha
    ss <- sample(nsamp, 100)
    alpha.post.pred.ord <- alpha.post.pred[order(minmax.post.pred),]

    colors <- terrain.colors(length(minmax.post.pred), alpha = 0.8)
    plot(bin.mids, 0*bin.mids, ty = "n", ylim = c(0,1), xlab = "time (ms)", ylab = "alpha", bty = "n")
    for(i in ss) lines(bin.mids, alpha.post.pred.ord[i,], col = colors[i])
    title(main = "Posterior predictive draws of Alpha")
    
    ## lambda.A and lambda.B estimates
    
    count2rate.factor <- 1000 / bin.width
    m.1 <- count2rate.factor * colMeans(fit$lambda.A); s.1 <- count2rate.factor * apply(fit$lambda.A, 2, sd)
    m.2 <- count2rate.factor * colMeans(fit$lambda.B); s.2 <- count2rate.factor * apply(fit$lambda.B, 2, sd)
    
    plot(bin.mids, 0*bin.mids, ylim = c(0, max(max(m.1+2*s.1), max(m.2+2*s.2))), ann = FALSE, ty = "n", bty = "n")
    polygon(bin.mids[c(1:nbins, nbins:1)], c(m.1 - 2*s.1, (m.1 + 2*s.1)[nbins:1]), col = tcol("orange", .5), border = tcol("orange", .5))
    polygon(bin.mids[c(1:nbins, nbins:1)], c(m.2 - 2*s.2, (m.2 + 2*s.2)[nbins:1]), col = tcol("cyan", .5), border = tcol("cyan", .5))
    abline(v = 0, lty = 2)
    title(xlab = "time (ms)", ylab = "firing rate (Hz)")
    title(main = "Single sound response")
    
    ## LS-Prob (to be removed by Azeem)
    #boxplot(lsprob.post, outline = FALSE, col = tcol(2, .3), border = tcol(2, .5), xlab = "length-scale (ms)", ylab = "probability", ylim = c(0,1))
    #if(add.prior) boxplot(lsprob.prior, outline = FALSE, col = tcol(4, .1), border = tcol(4, .5), add = TRUE)
    #title(main = "Length-scale distn")

    ## Minmax
    hist(minmax.post.pred, breaks = seq(0,1,.1), freq = FALSE, col = tcol(2, .3), border = tcol(2, .3), bty = "n", ann = FALSE)
    if(add.prior) hist(minmax.prior.pred, breaks = seq(0,1,.1), freq = FALSE, col = tcol(4, .1), border = tcol(4, .1), add = TRUE)
    title(xlab = "alpha minmax", ylab = "density", main = "Alpha MinMax")
    
    ## Average alpha.pred
    hist(rowMeans(alpha.post.pred), breaks = seq(0,1,.1), freq = FALSE, col = tcol(2, .3), border = tcol(2, .3), bty = "n", ann = FALSE)
    if(add.prior) hist(rowMeans(alpha.prior.pred), breaks = seq(0,1,.1), freq = FALSE, col = tcol(4, .1), border = tcol(4, .1), add = TRUE)
    title(xlab = "ave alpha", ylab = "density", main = "Ave Alpha")
    
    ## Swings
    swing.cuts <- c(0.2, 0.4, 0.6)
    nswing.prior.pred <- nswing.post.pred <- list()
    for(i in 1:length(swing.cuts)) nswing.post.pred[[i]] <- apply(alpha.post.pred, 1, swing.counter, width = swing.cuts[i])
    names(nswing.post.pred) <- swing.cuts
    boxplot(nswing.post.pred, xlab = "swing magnitude", ylab = "counts", col = tcol(2, .3), border = tcol(2, .5), pch = 19)

    if(add.prior){
        for(i in 1:length(swing.cuts)) nswing.prior.pred[[i]] <- apply(alpha.prior.pred, 1, swing.counter, width = swing.cuts[i])
        names(nswing.prior.pred) <- swing.cuts
        boxplot(nswing.prior.pred, col = tcol(4, .1), border = tcol(4, .5), pch = 19, add = TRUE)
    }
    title(main = "Swings")
    
}

dynamic.model.simu <- function(bin.mids = seq(12.5,1000,25), bin.width = 25, lengthScale = c(75, 125, 200, 300, 500), lsPrior = rep(1/length(lengthScale),length(lengthScale)), hyper = list(prec = c(1,1), sig0 = 1.87), nsamp = 1e3){
    
    sig0 <- hyper$sig0
    sigSq0 <- sig0^2
    prec.a <- hyper$prec[1]
    prec.b <- hyper$prec[2]
    
    nbins <- length(bin.mids)
    
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
        PREC[samp.count] <- prec <- rgamma(1, prec.a, prec.b)
        w1 <- prec; w2 <- 1
        LSPROB[samp.count,] <- ls.probs <- rDirichlet(1, lsPrior)
        psi <- rbeta(1, w1, w2)
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

rsynth <- function(ntrials = c(10, 10, 10), time.bins = 0:1000, lambda.A = 400, lambda.B = 100, pr.flat = 0.5, intervals = list(c(0,1)), wts = 1, span = c(0,1), period.range = c(400, 1000)){
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


bin.counter <- function(x, b) return(diff(sapply(b, function(a) sum(x <= a))))

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

rDirichlet <- function(n, alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- rowSums(x)
    return(x/sm)
}

rinvgamma <-function(n,shape, scale = 1) return(1/rgamma(n = n, shape = shape, rate = scale))

rtruncgamma <- function(n, low = 0, up = Inf, ...){
    unifs <- runif(n)
    plo <- pgamma(low, ...)
    pup <- pgamma(up, ...)
    return(qgamma(plo + unifs * (pup - plo), ...))
    #lplo <- pgamma(low, ..., log = TRUE)
    #lpup <- pgamma(up, ..., log = TRUE)
    #return(qgamma(lplo + log1p(unifs*(expm1(lpup - lplo))), ..., log.p = TRUE))
}

smoogam <- function(x){
    require(mgcv)
    T <- nrow(x)
    n <- ncol(x)
    if(T > 1){
        x.dd <- data.frame(cts = c(x), time = rep(1:T, n))
        x.gam <- gam(cts ~ s(time, bs = "ad"), data = x.dd, family = poisson(link = "log"))
        x.pred <- predict(x.gam, data.frame(cts = NA, time = 1:T), se.fit = TRUE)
        mu <- x.pred$fit; sig <- x.pred$se.fit
    } else {
        mu <- log(mean(c(x))); sig <- log(sd(c(x)))
    }
    firingrate.mean <- exp(mu)
    firingrate.vari <- expm1(sig^2/2) * firingrate.mean^2
    return(list(a = 1 / expm1(sig^2/2), b = 1/(expm1(sig^2/2) * firingrate.mean),
    mean = firingrate.mean, vari = firingrate.vari))
}

drop.item <- function(x, j) return(x[-match(j, x, nomatch = length(x) + 1)])


rank1up <- function(R, u){ # return chol(crossprod(R) + tcrossprod(u))
    n <- length(u)
    for(k in 1:n){
        r.old <- R[k,k]
        r.new <- sqrt(r.old^2 + u[k]^2)
        R[k, k] <- r.new
        if(k < n){
            v <- r.new / r.old
            s <- u[k] / r.old
            rest <- (k+1):n
            R[k,rest] <- (R[k,rest] + s * u[rest]) / v
            u[rest] <- v * u[rest] - s * R[k,rest]
        }
    }
    return(R)
}


swing.counter <- function(alpha, width){
    unique.vals <- sort(unique(alpha))
    n.unique <- length(unique.vals)
    swings <- rep(0, n.unique)

    ## first scan bottom up until out of width
    jj <- 1
    lo <- unique.vals[jj]
    up <- lo + width
    counter <- 0
    while(up <= unique.vals[n.unique]){
        counter <- counter + 1
        keys <- rep(0, length(alpha))
        keys[alpha >= up] <- 0.5
        keys[alpha <= lo] <- -0.5
        key.string <- keys[keys != 0]
        jumps <- diff(key.string)
        swings[counter] <- sum(abs(jumps))
        
        jj <- jj + 1
        lo <- unique.vals[jj]
        up <- lo + width
    }

    ## then scan top down until out of width
    jj <- n.unique
    up <- unique.vals[jj]
    lo <- up - width
    while(lo >= unique.vals[1]){
        counter <- counter + 1
        keys <- rep(0, length(alpha))
        keys[alpha >= up] <- 0.5
        keys[alpha <= lo] <- -0.5
        key.string <- keys[keys != 0]
        jumps <- diff(key.string)
        swings[counter] <- sum(abs(jumps))
        
        jj <- jj - 1
        up <- unique.vals[jj]
        lo <- up - width
    }
    return(max(swings))
}

tcol <- Vectorize(function(col, alpha = 1) {x <- col2rgb(col)/255; return(rgb(x[1],x[2],x[3],alpha = alpha))}, "col")

#### Wrapper function for real data runs ####

fitter.fn <- function(triplet, triplet.meta, start.time = 0, end.time = 1000, bw = 25, on.reward = TRUE, go.by.soff = TRUE, save.figure = FALSE, save.out = FALSE, tilt = TRUE, data.path = "http://www2.stat.duke.edu/~st118/Jenni/STCodes/Data", local.pull = FALSE, save.path = "./", ...){
    
    if(save.out | save.figure) if(!dir.exists(save.path)) stop("Invalid directory for storing summary and/or figures")
    
    freq <- triplet.meta[triplet, "AltFreq"]
    angl <- triplet.meta[triplet, "AltPos"]
    cell <- triplet.meta[triplet, "PairId"]
    fname <- triplet.meta[triplet, "CellId"]
    file.tag <- paste0(fname, "_cell", cell, "_freq", freq, "_pos", angl)
    
    is.JA <- is.na(cell)
    
    encase <- ifelse(local.pull, file, url)
    sub.dir <- ifelse(is.JA, "JA", "VC")
    sub.exp <- ifelse(is.JA, "", paste0("_cell", cell))
    
    trials <- read.table(encase(paste0(data.path, "/", sub.dir, "/", fname, ".txt")))
    spiketimes <- read.table(encase(paste0(data.path, "/", sub.dir, "/", fname, sub.exp, "_spiketimes.txt")))
    
    if(is.JA) trials <- cbind(trials, 600)
    colnames(trials) <- c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF")
    colnames(spiketimes) = c("TRIAL2", "TIMES")
    
    attach(trials)
    attach(spiketimes)
    
    frq <- c(freq, 742); pos <- c(angl, -144/angl)
    
    timestamps <- split(TIMES, TRIAL2)
    ntrials <- length(timestamps)
    trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
    
    ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
    ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
    ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
    
    if(on.reward){
        ix1 <- ix1 & REWARD == 1
        ix2 <- ix2 & REWARD == 1
        ix3 <- ix3 & REWARD == 1
    }
    
    sing1 <- trials[ix1, 1]
    sing2 <- trials[ix2, 1]
    success <- REWARD[ix3]
    
    if(go.by.soff) end.time <- min(SOFF[ix1 | ix2 | ix3])
    
    if(is.nan(end.time)){
        detach(trials)
        detach(spiketimes)
        stop("SOFF is NaN")
    }
    
    brk <- seq(start.time, end.time, bw)
    mpt <- (brk[-1] + brk[-length(brk)]) / 2
    
    spike.bincounter <- function(jj, brk){
        jj1 <- match(jj, trial.id)
        spks <- timestamps[[jj1]]
        return(bin.counter(spks, brk))
    }
    
    Acounts <- matrix(sapply(sing1, spike.bincounter, brk = brk), nrow = length(mpt))
    Bcounts <- matrix(sapply(sing2, spike.bincounter, brk = brk), nrow = length(mpt))
    nA <- ncol(Acounts)
    nB <- ncol(Bcounts)
    min.samp.size <- min(nA, nB)
    
    duplx <- trials[ix3, 1]
    ABcounts <- matrix(sapply(duplx, spike.bincounter, brk = brk), nrow = length(mpt))
    nAB <- ncol(ABcounts)
    
    detach(trials)
    detach(spiketimes)
    
    spike.counts <- list(Acounts = Acounts, Bcounts = Bcounts, ABcounts = ABcounts, bin.mids = mpt, bin.width = bw)
    fit.post <- dynamic.model.fit(spike.counts, ...)
    if(save.out){
        out.path <- paste0(save.path, "/Summaries/")
        if(!dir.exists(out.path)) dir.create(out.path)
        save(fit.post, file = paste0(out.path, file.tag, ".Rd"))
    }
    
    if(save.figure){
        figure.path <- paste0(save.path, "/Figures/")
        if(!dir.exists(figure.path)) dir.create(figure.path)
        pdf(file = paste0(figure.path, file.tag, ".pdf"), height = 9, width = 7)
    }
    plot(fit.post, tilt = tilt, nprior = 1e4)
    if(save.figure) dev.off()
    invisible(fit.post)
    
}
