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



source("dynamic_neural_model.R")
triplet.meta <- read.csv("~/Desktop/Research/Neuro/363triplets.csv")
local.path <- "~/Desktop/Research/Neuro/"
web.url <- "http://www2.stat.duke.edu/~st118/Jenni/STCodes/"

require(parallel)
trip.set <- 1:363

BW <- 50
data.path <- paste0(local.path, "Data")
save.path <- paste0(local.path, "ResultsBW", BW)
if(!dir.exists(save.path)) dir.create(save.path)
if(!dir.exists(paste0(save.path, "/Figures"))) dir.create(paste0(save.path, "/Figures"))
if(!dir.exists(paste0(save.path, "/Summaries"))) dir.create(paste0(save.path, "/Summaries"))

all.set <- mclapply(trip.set, function(jj) try(fitter.fn(jj, triplet.meta, plot = FALSE, verbose = FALSE, bw = BW, save.figure = TRUE, save.out = TRUE, data.path = data.path, local.pull = TRUE, save.path = save.path)), mc.cores = 4)
trip.err <- (sapply(all.set, class) == "try-error")


### Post processing for summary and plotting

file.tag <- rep(NA, length(trip.set))
for(jj in 1:length(trip.set)){
  triplet <- trip.set[jj]
  freq <- triplet.meta[triplet, "AltFreq"]
  angl <- triplet.meta[triplet, "AltPos"]
  cell <- triplet.meta[triplet, "PairId"]
  fname <- triplet.meta[triplet, "CellId"]
  file.tag[jj] <- paste0(fname, "_cell", cell, "_freq", freq, "_pos", angl)
}

triplet.combo <- cbind(triplet.meta, tag = file.tag)

summary.files <- list.files(paste0(save.path, "/Summaries"))
nfiles <- length(summary.files)
ncuts <- 4
summary.meta <- as.data.frame(matrix(nrow = nfiles, ncol = 1+4*ncuts))
names(summary.meta) <- c("tag", paste0("mmPriorQ", 1:ncuts), paste0("mmPostQ", 1:ncuts), paste0("avPriorQ", 1:ncuts), paste0("avPostQ", 1:ncuts))
for(jj in 1:nfiles){
  ff <- summary.files[jj]
  load(paste0(save.path, "/Summaries/", ff))
  ff.tag <- substring(ff, 1, nchar(ff)-3)
  sm <- summary(fit.post, tilt = TRUE, cut = 1/ncuts, mesh = 1/ncuts, nprior = 1e4)
  summary.meta[jj,1] <- ff.tag
  summary.meta[jj,-1] <- c(sm)
  rm(fit.post)
}

select.tags <- c("YKIC092311_1Loc_DoubleSound_cellNA_freq609_pos24",
                 "YKIC140121Loc_DoubleSound_cell1_freq500_pos24",
                 "YKIC141016Loc_DoubleSound_cell1_freq500_pos-24",
                 "YKIC080311Loc_DoubleSound_cellNA_freq1100_pos24",
                 "YKIC140501Loc_DoubleSound_cell1_freq500_pos24",
                 "YKIC140226Loc_DoubleSound_cell1_freq1988_pos-6",
                 "YKIC102711Loc_DoubleSound_cellNA_freq1100_pos24"
)

triplet.final <- merge(triplet.combo, summary.meta, by = "tag", all = TRUE)
triplet.final$wavieness <- with(triplet.final, (mmPostQ3 + mmPostQ4)/(mmPostQ1 + mmPostQ2))
triplet.final$centrality <- with(triplet.final, (avPostQ2 + avPostQ3)/(avPostQ1 + avPostQ4))
triplet.final$skewness <- with(triplet.final, pmax((avPostQ3 + avPostQ4)/(1e-3 + avPostQ1 + avPostQ2), (avPostQ1 + avPostQ2)/(1e-3 + avPostQ3 + avPostQ4)))
triplet.final$scaledskew <- log(pmax(1, pmin(triplet.final$skewness, 100))) / log(1e2)
triplet.final$chosen.tag <- ""
select.triplets <- sapply(select.tags, match, table = triplet.final$tag)
triplet.final$chosen.tag[select.triplets] <- as.character(triplet.final$tag[select.triplets])
ggplot(triplet.final, aes(centrality, wavieness)) + geom_point(aes(color = rank(skewness), shape = WinModels), cex = 4, alpha = 0.7) + scale_colour_gradientn(colours=rainbow(2)) + geom_hline(yintercept = 1, lty = 2, col = "gray") + geom_vline(xintercept = 7/3, lty = 2, col = "gray") + geom_text(aes(label = chosen.tag), size = 2, vjust = -.01, hjust = -.01)


three.tags <- list(wavieness = c("flat", "", "wavy"),
                   centrality = c("extreme", "", "central"),
                   skewness = c("symmetric", "", "skewed"))
pr.cuts <- list(wavieness = c(0, .78, 1.3, Inf), centrality = c(0,1.68,3.24,Inf), skewness = c(0,2,4,Inf))
#Prior wavieness = 1, centrality = 2.33 and skewness = 1.
# The cut-offs used above are resonable extremeties -- but could be altered.

three.summaries <- data.frame(matrix(NA, nrow = nrow(triplet.final), ncol = 3))
names(three.summaries) <- names(three.tags)

for(sm.name in c("centrality", "wavieness", "skewness")){
  groups <- cut(triplet.final[,sm.name], pr.cuts[[sm.name]], labels = FALSE, include.lowest = TRUE)
  three.summaries[,sm.name] <- three.tags[[sm.name]][groups]
}
triplet.final$dynamic.analysis.tag <- apply(three.summaries, 1, paste, collapse = "-")

whole3summ <- with(triplet.final, table(WinModels, dynamic.analysis.tag))
whole3summ <- whole3summ / rowSums(whole3summ)
whole3summ <- whole3summ[,apply(whole3summ, 2, max) >= 0.1]
whole3summ <- cbind(whole3summ, other = 1 - rowSums(whole3summ))
class(whole3summ) <- "table"
whole3summ.df <- as.data.frame(whole3summ)
names(whole3summ.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = "All 363 triplets, BW = 50")

top.cut <- 0.95
whole3summ.top <- with(subset(triplet.final, WinPr > top.cut), table(WinModels, dynamic.analysis.tag))
whole3summ.top <- whole3summ.top / rowSums(whole3summ.top)
whole3summ.top <- whole3summ.top[,apply(whole3summ.top, 2, max, na.rm = TRUE) >= 0.1]
whole3summ.top <- cbind(whole3summ.top, other = 1 - rowSums(whole3summ.top))
class(whole3summ.top) <- "table"
whole3summ.top.df <- as.data.frame(whole3summ.top)
names(whole3summ.top.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.top.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = paste0("WinPr > ", top.cut, ", BW = 50"))

## whole trial summary
qplot(factor(WinModels), data = triplet.final, fill = factor(WinPr < 0.95))
qplot(factor(WinModels), data = triplet.final, fill = factor(SepBF < 7))


## dynamic analysis
ggplot(triplet.final, aes(dynamic.analysis.tag, fill=(WinPr < 0.95))) + geom_bar() + coord_flip()
ggplot(triplet.final, aes(dynamic.analysis.tag, fill=(SepBF > 7.0))) + geom_bar() + coord_flip()


sep.cut <- 14.0
whole3summ.top <- with(subset(triplet.final, SepBF > sep.cut), table(WinModels, dynamic.analysis.tag))
whole3summ.top <- whole3summ.top / rowSums(whole3summ.top)
whole3summ.top <- whole3summ.top[,apply(whole3summ.top, 2, max, na.rm = TRUE) >= 0.1]
whole3summ.top <- cbind(whole3summ.top, other = 1 - rowSums(whole3summ.top))
class(whole3summ.top) <- "table"
whole3summ.top.df <- as.data.frame(whole3summ.top)
names(whole3summ.top.df)[1:2] <- c("WinModels", "dynamic.analysis.tag")
ggplot(subset(whole3summ.top.df, Freq > 0.1), aes(WinModels, dynamic.analysis.tag)) + geom_tile(aes(fill = Freq)) + scale_fill_gradientn(colors = rainbow(4)) +  labs(title = paste0("WinPr > ", top.cut, ", BW = 50"))
