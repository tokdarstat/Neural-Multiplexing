pkgname <- "neuromplex"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "neuromplex-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('neuromplex')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bin.counter")
### * bin.counter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bin.counter
### Title: Bin Counting
### Aliases: bin.counter
### Keywords: programming

### ** Examples

## generate 20 AB trials, roughl half with flat weight curves
## with a constant intensity either in (0,.1) or in (0.9, 1)
## (equally likely). The remaining curves are sinusoidal
## that snake between 0.1 and 0.9 with a period randomly
## drawn between 500 and 1500

synth.data <- synthesis.dapp(ntrials = c(15, 20, 20), pr.flat = 1,
                             intervals = list(c(0,.1), c(.45,.55), c(.9,1)),
                             wts = c(1/3, 1/3, 1/3), span = c(.1,.9),
                             period = c(500, 1500))

spike.counts <- list()
breaks <- seq(0, 1e3, 25)
spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bin.counter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dapp")
### * dapp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dapp
### Title: Dynamic Admixture of Poisson Process
### Aliases: dapp
### Keywords: programming

### ** Examples

## Not run: 
##D synth.data <- synthesis.dapp(ntrials = c(15, 20, 20), pr.flat = 1,
##D                              intervals = list(c(0,.1), c(.45,.55), c(.9,1)),
##D                              wts = c(1/3, 1/3, 1/3), span = c(.1,.9),
##D                              period = c(500, 1500))
##D 
##D spike.counts <- list()
##D breaks <- seq(0, 1e3, 25)
##D spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
##D spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
##D spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)
##D spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
##D spike.counts$bin.width <- diff(breaks)[1]
##D 
##D fit.post <- dapp(spike.counts)
##D plot(fit.post, synth.data = synth.data)
##D ## reanalyze forcing a uniform prior on range(alpha)
##D plot(fit.post, synth.data = synth.data, tilt = TRUE)
##D print(summary(fit.post, tilt = TRUE))
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dapp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dapp.simulate")
### * dapp.simulate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dapp.simulate
### Title: Simulate from Dynamic Admixture of Poisson Process
### Aliases: dapp.simulate
### Keywords: programming

### ** Examples

## Not run: 
##D prior <- dapp.simulate(1000, 25)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dapp.simulate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.dapp")
### * plot.dapp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.dapp
### Title: Plotting Method for Dynamic Admixture of Poisson Process
### Aliases: plot.dapp
### Keywords: programming

### ** Examples

## Not run: 
##D synth.data <- synthesis.dapp(ntrials = c(15, 20, 20), pr.flat = 1,
##D                              intervals = list(c(0,.1), c(.45,.55), c(.9,1)),
##D                              wts = c(1/3, 1/3, 1/3), span = c(.1,.9),
##D                             period = c(500, 1500))
##D 
##D spike.counts <- list()
##D breaks <- seq(0, 1e3, 25)
##D spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
##D spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
##D spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)
##D spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
##D spike.counts$bin.width <- diff(breaks)[1]
##D 
##D fit.post <- dapp(spike.counts)
##D plot(fit.post, synth.data = synth.data)
##D ## reanalyze forcing a uniform prior on range(alpha)
##D plot(fit.post, synth.data = synth.data, tilt = TRUE)
##D print(summary(fit.post, tilt = TRUE))
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.dapp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("poisson.tests")
### * poisson.tests

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: poisson.tests
### Title: Poisson Tests for Whole Trial Spike Counts
### Aliases: poisson.tests
### Keywords: programming

### ** Examples

## Not run: 
##D nA <- 20; nB <- 15; nAB <- 25
##D muA <- 25; muB <- 40
##D Acounts <- rpois(nA, muA)
##D Bcounts <- rpois(nB, muB)
##D ABcounts <- rpois(nAB, sample(c(muA, muB), nAB, replace = TRUE))
##D poisson.tests(Acounts, Bcounts, ABcounts)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("poisson.tests", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.dapp")
### * summary.dapp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.dapp
### Title: Summary Method for Dynamic Admixture of Poisson Process
### Aliases: summary.dapp
### Keywords: programming

### ** Examples

## Not run: 
##D synth.data <- synthesis.dapp(ntrials = c(15, 20, 20), pr.flat = 1,
##D                              intervals = list(c(0,.1), c(.45,.55), c(.9,1)),
##D                              wts = c(1/3, 1/3, 1/3), span = c(.1,.9),
##D                              period = c(500, 1500))
##D 
##D spike.counts <- list()
##D breaks <- seq(0, 1e3, 25)
##D spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
##D spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
##D spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)
##D spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
##D spike.counts$bin.width <- diff(breaks)[1]
##D 
##D fit.post <- dapp(spike.counts)
##D plot(fit.post, synth.data = synth.data)
##D ## reanalyze forcing a uniform prior on range(alpha)
##D plot(fit.post, synth.data = synth.data, tilt = TRUE)
##D print(summary(fit.post, tilt = TRUE))
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.dapp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("synthesis.dapp")
### * synthesis.dapp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: synthesis.dapp
### Title: Simulate Multiplexing Data for DAPP Analysis
### Aliases: synthesis.dapp
### Keywords: programming

### ** Examples

## generate 20 AB trials, roughl half with flat weight curves
## with a constant intensity either in (0,.1) or in (0.9, 1)
## (equally likely). The remaining curves are sinusoidal
## that snake between 0.1 and 0.9 with a period randomly
## drawn between 500 and 1500

synth.data <- synthesis.dapp(ntrials = c(15, 20, 20), pr.flat = 1,
                             intervals = list(c(0,.1), c(.45,.55), c(.9,1)),
                             wts = c(1/3, 1/3, 1/3), span = c(.1,.9),
                             period = c(500, 1500))

spike.counts <- list()
breaks <- seq(0, 1e3, 25)
spike.counts$Acounts <- sapply(synth.data$spiketimes$A, bin.counter, b = breaks)
spike.counts$Bcounts <- sapply(synth.data$spiketimes$B, bin.counter, b = breaks)
spike.counts$ABcounts <- sapply(synth.data$spiketimes$AB, bin.counter, b = breaks)
spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
spike.counts$bin.width <- diff(breaks)[1]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("synthesis.dapp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
