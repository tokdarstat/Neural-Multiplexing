\name{summary.dapp}
\Rdversion{0.1}
\alias{summary.dapp}
\title{Summary Method for Dynamic Admixture of Poisson Process}
\description{Summarizes model fit of the DAPP model to binned spiking data}
\usage{ 
\method{summary}{dapp}(object, cut.width = 0.1, tilt.prior = FALSE,
     mesh.tilt = 0.1, nprior = object$mcmc["nsamp"], ...)
} 
\arguments{
\item{object}{a fitted model of the class 'dapp'}
\item{cut.width}{width of the bins to be used to cut range(alpha) and average(alpha) characteristics}
\item{tilt.prior}{lofical giving whether the prior should be tilted to mimic an analysis done with a uniform prior on the range(alpha)}
\item{mesh.tilt}{a tuning parameter that controls how exactly tilting is done. Shorter mesh value gives tighter match but will require more Monte Carlo simulations}
\item{nprior}{number of prior draws to be used for display}
\item{...}{no addiitonal parameters used at this point}
}
\value{
Gives prior and posterior summaries of the range and average predicted alpha curves
}

\details{To be added...}


\examples{
\dontrun{
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

fit.post <- dapp(spike.counts)
plot(fit.post, synth.data = synth.data)
## reanalyze forcing a uniform prior on range(alpha)
plot(fit.post, synth.data = synth.data, tilt = TRUE)
print(summary(fit.post, tilt = TRUE))

}
}
\keyword{programming}
