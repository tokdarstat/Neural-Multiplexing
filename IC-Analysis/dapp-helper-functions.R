#### Helper functions to run DAPP analysis on IC data ####

## Wrapper function for real data runs ##
## Allows data pulling from a local directory named 'Data' ##
## Or data can be pulled via the internet ##

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
    fit.post <- dapp(spike.counts, ...)
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

## Additional function for quick visualization of the data ##
## Creates three plots for the three conditions (A, B, and, AB)
## For each conditions, every trial associated with the condition
## are shown by a curve over time giving a smoothed spike train.
## A and B smoothed spike train bands are shown as benchmark for
## AB spike trains.
## multi.color = TRUE will show each AB trial with a different color
## group = TRUE will do a quick classification of each AB trial as
## 'A-like' or 'B-like' based on aggregate spike count, and,
## distinguish the groups by color

visual.fn <- function(triplet, triplet.meta, start.time = 0, end.time = 1000, bw = 25, on.reward = TRUE, go.by.soff = TRUE, data.path = "http://www2.stat.duke.edu/~st118/Jenni/STCodes/Data", local.pull = FALSE, retain.plot.par = FALSE, remove.zeros = FALSE, multi.color = TRUE, group = FALSE, ...){
    
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
    
    x1 <- spike.counts$Acounts
    x2 <- spike.counts$Bcounts
    x3 <- spike.counts$ABcounts
    bin.mids <- spike.counts$bin.mids
    bin.width <- spike.counts$bin.width
    
    if(remove.zeros){
        x1 <- x1[,colSums(x1) > 0, drop = FALSE]
        x2 <- x2[,colSums(x2) > 0, drop = FALSE]
        x3 <- x3[,colSums(x3) > 0, drop = FALSE]
    }
    
    nbins <- length(bin.mids)
    if(nrow(x3) != nbins) stop("dimension mismatch between spike counts and bins")
    
    n1 <- ncol(x1)
    n2 <- ncol(x2)
    n3 <- ncol(x3)
    
    #get.hyper1 <- smoogam(x1);
    #get.hyper2 <- smoogam(x2);
    #m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, nbins); s.1 <- sqrt(am.1)/bm.1
    #m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, nbins); s.2 <- sqrt(am.2)/bm.2
    
    count.to.Hz.factor <- 1000/bw
    
    x1.smu <- sapply(1:n1, function(i) supsmu(bin.mids, x1[,i])$y)
    x2.smu <- sapply(1:n2, function(i) supsmu(bin.mids, x2[,i])$y)
    x3.smu <- sapply(1:n3, function(i) supsmu(bin.mids, x3[,i])$y)
    
    x.max <- max(max(x1.smu), max(x2.smu), max(x3.smu))
    x.min <- min(min(x1.smu), min(x2.smu), min(x3.smu))
    
    m.1 <- apply(x1.smu, 1, mean); s.1 <- apply(x1.smu, 1, sd)/sqrt(n1)
    m.2 <- apply(x2.smu, 1, mean); s.2 <- apply(x2.smu, 1, sd)/sqrt(n2)
    
    if(multi.color) group <- FALSE
    
    resp.type <- rep(0, n3)
    p1 <- rep(0, n3)
    if(group){
        mu.1 <- sum(m.1)
        mu.2 <- sum(m.2)
        x3.tot <- colSums(x3)
        resp.type <- 1 + (dpois(x3.tot, mu.1, log = TRUE) < dpois(x3.tot, mu.2, log = TRUE))
        p1 <- 1/ (1 + dpois(x3.tot, mu.1, log = TRUE) / dpois(x3.tot, mu.2, log = TRUE))
    }
    
    if(!retain.plot.par) par(mfrow = c(1,3))
    plot(bin.mids, 0*bin.mids, ylim = count.to.Hz.factor * c(x.min, x.max), ann = FALSE, ty = "n", bty = "n")
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.1 - 2*s.1, (m.1 + 2*s.1)[nbins:1]), col = tcol("orange", .5), border = tcol("orange", .5))
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.2 - 2*s.2, (m.2 + 2*s.2)[nbins:1]), col = tcol("cyan", .5), border = tcol("cyan", .5))
    abline(v = 0, lty = 2)
    title(xlab = "time (ms)", ylab = "firing rate (Hz)", ...)
    title(main = "Stimulus A", line = 0, font.main = 1)
    #for(j in 1:ncol(x1)) lines(bin.mids, count.to.Hz.factor * smooth(x1[,j]), col = j)
    for(j in 1:ncol(x1)) lines(bin.mids, count.to.Hz.factor * x1.smu[,j], col = tcol(multi.color*(j-1) + 1 + group, 0.5), lwd = 1 + 3 * group)
    
    plot(bin.mids, 0*bin.mids, ylim = count.to.Hz.factor * c(x.min, x.max), ann = FALSE, ty = "n", bty = "n")
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.1 - 2*s.1, (m.1 + 2*s.1)[nbins:1]), col = tcol("orange", .5), border = tcol("orange", .5))
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.2 - 2*s.2, (m.2 + 2*s.2)[nbins:1]), col = tcol("cyan", .5), border = tcol("cyan", .5))
    abline(v = 0, lty = 2)
    title(xlab = "time (ms)", ylab = "firing rate (Hz)", ...)
    title(main = "Stimulus B", line = 0, font.main = 1)
    #  for(j in 1:ncol(x2)) lines(bin.mids, count.to.Hz.factor * smooth(x2[,j]), col = j)
    for(j in 1:ncol(x2)) lines(bin.mids, count.to.Hz.factor * x2.smu[,j], col = tcol(multi.color*(j-1) + 1 + 3*group, 0.5), lwd = 1 + 3 * group)
    
    plot(bin.mids, 0*bin.mids, ylim = count.to.Hz.factor * c(x.min, x.max), ann = FALSE, ty = "n", bty = "n")
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.1 - 2*s.1, (m.1 + 2*s.1)[nbins:1]), col = tcol("orange", .5), border = tcol("orange", .5))
    polygon(bin.mids[c(1:nbins, nbins:1)], count.to.Hz.factor * c(m.2 - 2*s.2, (m.2 + 2*s.2)[nbins:1]), col = tcol("cyan", .5), border = tcol("cyan", .5))
    abline(v = 0, lty = 2)
    title(xlab = "time (ms)", ylab = "firing rate (Hz)", ...)
    title(main = "Stimuli A+B", line = 0, font.main = 1)
    #for(j in 1:ncol(x3)) lines(bin.mids, count.to.Hz.factor * smooth(x3[,j]), col = j)
    #for(j in 1:ncol(x3)) lines(bin.mids, count.to.Hz.factor * x3.smu[,j], col = tcol(multi.color*(j-1) + 1 + 2*resp.type[j] - group, .5))
    for(j in 1:ncol(x3)) lines(bin.mids, count.to.Hz.factor * x3.smu[,j], col = tcol(rgb(group*p1[j],0,group*(1-p1[j])), 0.5), lwd = 1 + 3 * group)
    
    invisible(list(spikes = spike.counts, cell = file.tag))
}



