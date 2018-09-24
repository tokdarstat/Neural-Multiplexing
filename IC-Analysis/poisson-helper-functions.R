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
  
  require(neuromplex)

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
  
  return(poisson.tests(Acounts, Bcounts, ABcounts, labels = c(s1, s2, dp), ...))
}


poi.VC <- function(fname, cell, data.path = "Data/", on.reward = TRUE,
                   match.level = FALSE, AB.eqlevel = FALSE, outfile = "",
                   start = 0, end = 600, go.by.soff = TRUE, remove.zeros = FALSE,
                   plot = TRUE){
  
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


poi.JA <- function(fname, on.reward = TRUE, data.path = "Data/",
                   match.level = FALSE, AB.eqlevel = FALSE, outfile = "",
                   start = 0, end = 600, remove.zeros = FALSE,
                   plot = TRUE){
  
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




