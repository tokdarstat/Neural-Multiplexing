## Display raw spiketrains from JA data

display.raw.JA <- function(fname, path = "./"){
  data.path <- paste0(path, "JA/")
    
  infile1 <- paste(data.path, fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL"))
  
  infile2 <- paste(data.path, fname, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps <- split(TIMES, TRIAL2)
  ntrials <- length(timestamps)
  
  trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
  
  layout(matrix(c(rep(1,12),rep(2,4)), 4, 4, byrow = FALSE))
  
  plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(TRIAL), bty = "n", xlab = "time", ylab = "trial")  
  reward.locator <- rep(NA, ntrials)
  for(i in 1:ntrials){
    jj <- trial.id[i]
    jj2 <- which(TRIAL == jj)
    spks <- timestamps[[i]]
    nspk <- length(spks)
    reward.locator[i] <- REWARD[jj2]
    if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "darkgreen", "red")[1 + (spks > 0)])  
  }
  
  spk.counts <- sapply(timestamps, length)
  plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
  axis(1)
  segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "royalblue")[1 + reward.locator])
  
  title(main = fname, out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}


## Display raw spike trains from VC data

display.raw.VC <- function(fname, cell, path = "./"){
  data.path <- paste0(path, "VC/")
  
  infile1 <- paste(data.path, fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste(data.path, fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps <- split(TIMES, TRIAL2)
  ntrials <- length(timestamps)  
  trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
  
  layout(matrix(c(rep(1,12),rep(2,4)), 4, 4, byrow = FALSE))
  
  plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(TRIAL), bty = "n", xlab = "time", ylab = "trial")  
  reward.locator <- rep(NA, ntrials)
  for(i in 1:ntrials){
    jj <- trial.id[i]
    jj2 <- which(TRIAL == jj)
    spks <- timestamps[[i]]
    nspk <- length(spks)
    reward.locator[i] <- REWARD[jj2]
    if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "red", "orange")[1 + (spks > 0) + (spks > SOFF[jj2])])  
  }

  spk.counts <- sapply(timestamps, length)
  plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
  axis(1)
  segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "magenta")[1 + reward.locator])
  
  title(main = paste(fname, " (cell ", cell, ")", sep = ""), out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}

## display raw spikes from a triplet
display.triplet.JA <- function(fname, Afrq, Apos, on.reward = FALSE, path = "./"){
  data.path <- paste0(path, "JA/")
  
  infile1 <- paste(data.path, fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL"))
  trials$SOFF <- 600
  
  infile2 <- paste(data.path, fname, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps.all <- split(TIMES, TRIAL2)
  
  frq <- c(Afrq, 742); pos <- c(Apos, -144/Apos)
  ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
  ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
  ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
  
  if(on.reward){
    ix1 <- ix1 & REWARD == 1
    ix2 <- ix2 & REWARD == 1
    ix3 <- ix3 & REWARD == 1
  } 
  groups <- list(trials[ix1, 1], trials[ix2, 1], trials[ix3, 1])
  
  layout(matrix(rep(1:6, c(2,1,2,1,2,1)), 3,3, byrow = TRUE))
  
  group.names <- c(paste("A @", frq[1], "Hz", pos[1], "deg", sep = ""),
                   paste("B @", frq[2], "Hz", pos[2], "deg", sep = ""),
                   "AB")
  for(gg in 1:3){
    gg.sel <- na.omit(match(groups[[gg]], names(timestamps.all)))
    timestamps <- timestamps.all[gg.sel]
    ntrials <- length(timestamps)  
    trial.id <- as.numeric(names(timestamps.all))[gg.sel] ## same as unique(TRIAL2)
    
    plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(trial.id), bty = "n", xlab = "time", ylab = "trial")  
    reward.locator <- rep(NA, ntrials)
    for(i in 1:ntrials){
      jj <- trial.id[i]
      jj2 <- which(TRIAL == jj)
      spks <- timestamps[[i]]
      nspk <- length(spks)
      reward.locator[i] <- REWARD[jj2]
      if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "darkgreen", "red")[1 + (spks > 0) + (spks > SOFF[jj2])])  
    }
    title(main = group.names[gg], font = 1, line = 0)    
    
    spk.counts <- sapply(timestamps, length)
    plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
    axis(1)
    segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "royalblue")[1 + reward.locator])
  }
  title(main = paste(fname, " (cell ", cell, ")", sep = ""), out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}


## display raw spikes from a triplet
display.triplet.VC <- function(fname, cell, Afrq, Apos, on.reward = FALSE, path = "./"){
  data.path <- paste0(path, "VC/")
  
  infile1 <- paste(data.path, fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste(data.path, fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps.all <- split(TIMES, TRIAL2)
  
  frq <- c(Afrq, 742); pos <- c(Apos, -144/Apos)
  ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
  ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
  ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
  
  if(on.reward){
    ix1 <- ix1 & REWARD == 1
    ix2 <- ix2 & REWARD == 1
    ix3 <- ix3 & REWARD == 1
  } 
  groups <- list(trials[ix1, 1], trials[ix2, 1], trials[ix3, 1])
  
  layout(matrix(rep(1:6, c(2,1,2,1,2,1)), 3,3, byrow = TRUE))
  
  group.names <- c(paste("A @", frq[1], "Hz", pos[1], "deg", sep = ""),
                   paste("B @", frq[2], "Hz", pos[2], "deg", sep = ""),
                   "AB")
  for(gg in 1:3){
    gg.sel <- na.omit(match(groups[[gg]], names(timestamps.all)))
    timestamps <- timestamps.all[gg.sel]
    ntrials <- length(timestamps)  
    trial.id <- as.numeric(names(timestamps.all))[gg.sel] ## same as unique(TRIAL2)
    
    plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(trial.id), bty = "n", xlab = "time", ylab = "trial")  
    reward.locator <- rep(NA, ntrials)
    for(i in 1:ntrials){
      jj <- trial.id[i]
      jj2 <- which(TRIAL == jj)
      spks <- timestamps[[i]]
      nspk <- length(spks)
      reward.locator[i] <- REWARD[jj2]
      if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "darkgreen", "red")[1 + (spks > 0) + (spks > SOFF[jj2])])  
    }
    title(main = group.names[gg], font = 1, line = 0)    
    
    spk.counts <- sapply(timestamps, length)
    plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
    axis(1)
    segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "royalblue")[1 + reward.locator])
  }
  title(main = paste(fname, " (cell ", cell, ")", sep = ""), out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}

## Poisson analysis customized for spike train triplets

esti.Poi <- function(trials, spiketimes, frq = c(1100, 742), pos = c(24, -6), on.reward = TRUE, start.time = 0, end.time = 600, match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = FALSE, ...){
  
  source("poisson_analysis.R")

  attach(trials)
  attach(spiketimes)
  
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
  
  blev <- sort(unique(B_LEVEL[ix3]))
  targ.lev <- blev[blev > 0]
  lev <- "*"
  if(match.level) {
    if(length(targ.lev) > 1) {
      targ.lev <- max(targ.lev)
      warning("Multiple single sound levels, choosing largest one")
    }
    ix1 <- ix1 & A_LEVEL == targ.lev
    ix2 <- ix2 & A_LEVEL == targ.lev
    lev <- as.character(targ.lev)
  }

  if(AB.eqlevel) ix3 <- ix3 & (A_LEVEL == B_LEVEL)
  
  sing1 <- trials[ix1, 1]
  sing2 <- trials[ix2, 1]
  duplx <- trials[ix3, 1]    
  success <- REWARD[ix3]
  
  if(go.by.soff) end.time <- min(SOFF[ix1 | ix2 | ix3])
  
  spike.counter <- function(jj){
    jj1 <- match(jj, trial.id)
    spks <- timestamps[[jj1]]
    return(sum(spks > start.time & spks < end.time))
  }
  
  Acounts <- sapply(sing1, spike.counter)
  Bcounts <- sapply(sing2, spike.counter)
  ABcounts <- sapply(duplx, spike.counter)
  
  detach(trials)
  detach(spiketimes)
  
  
  s1 <- paste(frq[1],"Hz ",pos[1],"deg ",lev,"Db ", sep = "")
  s2 <- paste(frq[2],"Hz ",pos[2],"deg ",lev,"Db ", sep = "")
  dp <- paste("Duplex: ", lev, "Db ", sep = "")
  
  return(get.bayes.factors(Acounts, Bcounts, ABcounts, labels = c(s1, s2, dp), ...))   
}


poi.VC <- function(fname, cell, on.reward = TRUE, match.level = FALSE, AB.eqlevel = FALSE, outfile = "", start = 0, end = 600, go.by.soff = TRUE, remove.zeros = FALSE, plot = TRUE){
  
  infile1 <- paste(data.path, "/VC/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste(data.path, "/VC/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])  
  alt.pos <- c(-24, -6, 6, 24)
  
  par(mfcol = c(length(alt.freq),2), mar = c(2,3,3,0) + .1)
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({lbf <- round(esti.Poi(trials, spiketimes, c(fr, 742), c(po, -144/po), on.reward, start, end, match.level, AB.eqlevel, go.by.soff, remove.zeros), 4);
        cat(fname, cell, c(fr, po, lbf), "\n", file = outfile, append = TRUE)});
        if(plot) text(0, 0.09, paste0(fname, "_cell", cell), cex = 0.7, pos = 4)    
    }
  }
}


poi.JA <- function(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = FALSE, outfile = "", start = 0, end = 600, remove.zeros = FALSE, plot = TRUE){
  
  infile1 <- paste(data.path, "/JA/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL"))
  
  infile2 <- paste(data.path, "/JA/", fname, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])
  alt.pos <- c(-24, -6, 6, 24)
  
  par(mfcol = c(length(alt.freq),2), mar = c(2,3,3,0) + .1)
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({lbf <- round(esti.Poi(trials, spiketimes, c(fr, 742), c(po, -144/po), on.reward, start, end, match.level, AB.eqlevel, FALSE, remove.zeros, plot), 4);
        cat(fname, c(fr, po, lbf), "\n", file = outfile, append = TRUE);
        if(plot) text(0, 0.09, fname, cex = 0.7, pos = 4)    
      })
    }  
  }
}

## HMM analysis ###

esti.hmm <- function(trials, spiketimes, frq = c(1100, 742), pos = c(24, -6), on.reward = TRUE, start.time = 0, end.time = 600, bw = 25, target = c(25, 150, 800), match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = FALSE, n.iter = 1e3, plot = FALSE, faux.dual.mix = FALSE, faux.dual.int = FALSE, faux.alpha = 0.5, faux.dual.swi = FALSE, faux.target = c(100, 100), nAB = "match.realABcount", ...){
  attach(trials)
  attach(spiketimes)
  
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
  
  if(match.level) {
    ix1 <- ix1 & A_LEVEL == lev[1]
    ix2 <- ix2 & A_LEVEL == lev[1]
    ix3 <- ix3 & A_LEVEL == lev[2]
  }
  
  if(AB.eqlevel) ix3 <- ix3 & (A_LEVEL == B_LEVEL)

  if(min(sum(ix1), sum(ix2), sum(ix3)) == 0){
    detach(trials)
    detach(spiketimes)
    stop("Not enought data")
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
  
  Abincounts <- matrix(sapply(sing1, spike.bincounter, brk = brk), nrow = length(mpt))
  Bbincounts <- matrix(sapply(sing2, spike.bincounter, brk = brk), nrow = length(mpt))
  nA <- ncol(Abincounts)
  nB <- ncol(Bbincounts)
  min.samp.size <- min(nA, nB)

  nAB.target <- ifelse(nAB == "match.realABcount", sum(ix3), nAB)
  duplx <- trials[ix3, 1]
  dualbincounts <- matrix(sapply(duplx, spike.bincounter, brk = brk), nrow = length(mpt))
  if(faux.dual.mix){
    A.subsamp <- sample(nA, min(ceiling(nAB.target*faux.alpha), nA))
    B.subsamp <- sample(nB, min(ceiling(nAB.target*(1-faux.alpha)), nB))
    ABbincounts <- cbind(Abincounts[,A.subsamp], Bbincounts[,B.subsamp])
  } else if (faux.dual.int) {
    ABbincounts <- replicate(nAB.target, thincombo(Abincounts[,sample(nA,1)], Bbincounts[,sample(nB,1)], faux.alpha))
  } else if (faux.dual.swi) {
    #faux.theta <- bw / pmax(bw, faux.target)
    ABbincounts <- replicate(nAB.target, swicombo(Abincounts[,sample(nA,1)], Bbincounts[,sample(nB,1)], bw = bw, p1 = faux.alpha, mean.stay = faux.target, sd.stay = 0.25 * faux.target))    
  } else {
    ABbincounts <- dualbincounts
  }
  
  detach(trials)
  detach(spiketimes)
  
  nAB <- ncol(ABbincounts)
  
  source("multiplex-new.R")
  theta.hyper <- match.target(bw, target)
  #theta.hyper <- c(1/2, 1/2)
  #print(theta.hyper)
  a3.1 <- a3.2 <- theta.hyper[1]
  b3.1 <- b3.2 <- theta.hyper[2]
  
  oo <- main.fn(Abincounts, Bbincounts, ABbincounts, n.iter = n.iter, mpt = mpt, a3.1=a3.1,b3.1=b3.1,a3.2=a3.2,b3.2=b3.2,...)
  ss <- unique(ceiling(n.iter * seq(0.1, 1, length = 1000)))
  p.swi <- colMeans(oo$ptrial.swi.samp[ss,])
  p.div <- 1 - p.swi
  p.type <- cbind(p.swi, p.div)
  trial.type.winr <- apply(p.type, 1, which.max)
  p.max <- apply(p.type, 1, max)
  p.s1 <- colMeans(oo$state1.samp[ss,])
  
  groups <- list(switch.red = which(p.swi == p.max & p.s1 > 0.5),
                 switch.blue = which(p.swi == p.max & p.s1 <= 0.5),
                 divide = which(p.div == p.max))
  
  prob.stateA <- apply(oo$states.samp[ss,,,drop=FALSE] == 1, c(2,3), mean)
  p.A <- apply(prob.stateA, 2, median)
  
  prob.1state <- mean(oo$p1.samp[ss])
  prob.statechange <- 1 - colMeans(oo$theta.samp)
  
  pcell.swi <- oo$p.swi.samp[ss]; pcell.div <- 1 - pcell.swi
  pcell.type <- cbind(pcell.swi, pcell.div)
  prob.half <- colMeans(pcell.type > 0.5)
  prob.winr <- count(apply(pcell.type, 1, which.max), 1:2) / length(ss)
  swi.rate <- mean(pcell.swi); div.rate <- mean(pcell.div)
  
  pval.totfire <- min(apply(cbind(Abincounts, Bbincounts), 1, function(z) wilcox.test(z[1:nA], z[nA+1:nB], exact = FALSE)$p.val))
  pbin.sep <- eqpoi.test(matrix(Abincounts, ncol = nA), matrix(Bbincounts, ncol = nB))
  pw.est <- colMeans(oo$pw.samp[ss,])
  

  if(plot){
    
    m.1 <- colMeans(oo$mu1.samp); s.1 <- apply(oo$mu1.samp, 2, sd)
    m.2 <- colMeans(oo$mu2.samp); s.2 <- apply(oo$mu2.samp, 2, sd)
    
    fct <- 1e3 / bw ## count to rate factor
    T <- nrow(ABbincounts)
    require(KernSmooth)

    Tlong <- (end.time - start.time)
    tptlong <- start.time + 1:Tlong - 0.5

    Asmoothed <- apply(Abincounts, 2, function(z) locpoly(mpt, z, bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y)
    Bsmoothed <- apply(Bbincounts, 2, function(z) locpoly(mpt, z, bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y)
    ABsmoothed <- apply(ABbincounts, 2, function(z) locpoly(mpt, z, bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y)
    dualsmoothed <- apply(dualbincounts, 2, function(z) locpoly(mpt, z, bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y)
    rate.max <- max(c(max(Asmoothed), max(Bsmoothed), max(dualsmoothed)))  
    
    
    par(mfrow = c(3,2))  
    for(j in 1:3){
      plot(mpt, fct*ABbincounts[,1], ylim = fct*c(0,1.0 * rate.max), ann = FALSE, ty = "n", bty = "n")
      polygon(mpt[c(1:T, T:1)], fct * c(m.1 - 2*s.1, (m.1 + 2*s.1)[T:1]), col = '#AA000066', border = '#AA000066')
      polygon(mpt[c(1:T, T:1)], fct * c(m.2 - 2*s.2, (m.2 + 2*s.2)[T:1]), col = '#0000AA66', border = '#0000AA66')      
      for(k in groups[[j]]) lines(tptlong, fct*ABsmoothed[,k], lwd=.33, col = k)
      abline(v = 0, lty = 2)
      title(xlab = "Time (ms)", ylab = "Firing rate (Hz)")
      title(main = paste(names(groups)[j], " (", length(groups[[j]]), ")", sep = ""))
    }
    plot(mpt, fct*ABbincounts[,1], ylim = fct*c(0,1.0 * rate.max), ann = FALSE, ty = "n", bty = "n")
    polygon(mpt[c(1:T, T:1)], fct * c(m.1 - 2*s.1, (m.1 + 2*s.1)[T:1]), col = '#AA000066', border = '#AA000066')
    polygon(mpt[c(1:T, T:1)], fct * c(m.2 - 2*s.2, (m.2 + 2*s.2)[T:1]), col = '#0000AA66', border = '#0000AA66')      
    for(k in 1:nA) lines(tptlong, fct*Asmoothed[,k], lwd=.33, col = 2)
    for(k in 1:nB) lines(tptlong, fct*Bsmoothed[,k], lwd=.33, col = 4)
    abline(v = 0, lty = 2)
    title(xlab = "Time (ms)", ylab = "Firing rate (Hz)")
    title(main = paste("singles (nA = ", nA, ",nB = ", nB, ")", sep = ""))
    #hh <- hist(c(oo$ptrial.w1.samp), breaks = 0:11/10 - 0.05, plot = FALSE)
    #plot(hh$mids, 10 * hh$dens, ty = "h", ylim = c(0, 10*max(hh$dens)), bty = "n", ann = FALSE, col = "darkgreen", lwd = 4)
    #abline(h = seq(0,100,10), col = "lightblue", lty = 1)
    #points(hh$mids, 10 * hh$dens, pch = 15, col = "darkgreen")
    plot(oo$w.grid, 100*pw.est, ty = "h", bty = "n", ann = FALSE, col = "darkgreen", lwd = 4)
    points(oo$w.grid, 100*pw.est, pch = 15, col = "darkgreen")
    title(xlab = "alpha", ylab = "% Trials", main = "attention division")
    #cat("fine so far")
    
    try({
      fano.A <- locpoly(mpt, apply(Abincounts, 1, fano), bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y;
      fano.B <- locpoly(mpt, apply(Bbincounts, 1, fano), bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y;
      fano.AB <- locpoly(mpt, apply(ABbincounts, 1, fano), bandwidth=max(bw,10), gridsize=Tlong,range.x=range(tptlong))$y;
      plot(tptlong, fano.A, ann = FALSE, col = 2, ty = "l", bty = "n", ylim = range(c(range(fano.A), range(fano.B), range(fano.AB))), lwd = 2);
      lines(tptlong, fano.B, col = 4, lwd = 2);
      lines(tptlong, fano.AB, col = "darkgreen", lwd = 2);
      title(xlab = "Time (ms)", ylab = "Fano factor", main = "Fano")
    })
    title(paste(paste(frq, pos, sep = "@", collapse = " "), ": P(switch / divide) =", paste(round(c(swi.rate, div.rate), 2), collapse = " / ")), outer = TRUE, line = -1)    
    
  }
  return(list(probs = c(swi.rate, div.rate, prob.1state, prob.statechange, prob.winr, pval.totfire, pbin.sep, min.samp.size, round(pw.est,2)), best.trial.label = trial.type.winr, prob.best.label = p.max, pstateA = prob.stateA))
}


hmm.JA <- function(fname, outfile = "", ...){
  
  infile1 <- paste("JA/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL"))
  
  infile2 <- paste("JA/", fname, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])
  alt.pos <- c(-24, -6, 6, 24)
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({hrun <- esti.hmm(trials, spiketimes, c(fr, 742), c(po, -144/po), ...);
           probs <- round(hrun$probs, 2); 
           cts <- sapply(1:3, function(j) bin.counter(hrun$prob.best.label[hrun$best.trial.label == j], 0:5/5));
           cat(fname, c(fr, po, probs, cts), "\n", file = outfile, append = TRUE)})
    }
  }
}


hmm.VC <- function(fname, cell, outfile = "", ...){
  infile1 <- paste("VC/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste("VC/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))

  FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
  alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])
  alt.pos <- c(-24, -6, 6, 24)
  
  for(fr in alt.freq){
    for(po in alt.pos){
      try({hrun <- esti.hmm(trials, spiketimes, c(fr, 742), c(po, -144/po), ...);
           probs <- round(hrun$probs, 2); 
           cts <- sapply(1:3, function(j) bin.counter(hrun$prob.best.label[hrun$best.trial.label == j], 0:5/5));
           cat(fname, cell, c(fr, po, probs, cts), "\n", file = outfile, append = TRUE)
           })
    }
  }
}


### AUX functions ###
log.pm <- function(x, a, b){
  a.x <- a + sum(x)
  b.x <- b + length(x)
  mu.map <- a.x / b.x
  return(sum(dpois(x, mu.map, log = TRUE)) - diff(dgamma(mu.map, c(a, a.x), c(b, b.x), log = TRUE)))
}


eqpoi.test <- function(x1, x2){
  T <- nrow(x1)
  n1 <- ncol(x1)
  n2 <- ncol(x2)
  
  a1 <- mean(c(x1)); b1 <- 1
  a2 <- mean(c(x2)); b2 <- 1
  a3 <- mean(c(c(x1), c(x2))); b3 <- 1
  lps <- sapply(1:T, function(j) return(log.pm(x1[j,],a1,b1) + log.pm(x2[j,],a2,b2) - log.pm(c(x1[j,],x2[j,]),a3,b3)))
  return(prod(1 / (1 + exp(lps)/9)))
}

bin.counter <- function(x, b) return(diff(sapply(b, function(a) sum(x <= a))))

thincombo <- function(x, y, alpha){
  n <- length(x)
  if(length(y) != n) stop("x and y lengths must be same")
  xthin <- rbinom(n,x,alpha)
  ythin <- rbinom(n,y,1-alpha)
  return(xthin + ythin)
}
swicombo <- function(x, y, p1, bw = 25, mean.stay = c(100, 100), sd.stay = c(25,25)){
  n <- length(x)
  mean.stay <- mean.stay/bw
  sd.stay <- sd.stay / bw
  prob <- pmin(mean.stay/sd.stay^2, 0.99)
  size <- mean.stay * prob / (1 - prob)
  if(length(y) != n) stop("x and y lengths must be same")
  start.is.2 <- (runif(1) > p1)
  stays <- 1 + rnbinom(n + 1, size, prob)[start.is.2 + 1:n]
  state.order <- rep(c(1,2),n)[start.is.2 + 1:n]
  states <- rep(state.order, stays)[1:n]
  z <- x
  z[states == 2] <- y[states == 2]
  return(z)
}

match.target <- function(bw, target, maxit = 100){
  ix <- (target > bw)
  f <- function(par) return(sum((1 - bw/(target*qbeta(c(.975, .5, .025), exp(par[1]), exp(par[2])))[ix])^2))
  f.op <- optim(log(c(2.5, 7.5)), f, control=list(maxit = maxit))
  ab <- exp(f.op$par)
  attr(ab, "conv:code") <- f.op$convergence
  return(ab)
}

fano <- function(x) return(var(x) / mean(x))