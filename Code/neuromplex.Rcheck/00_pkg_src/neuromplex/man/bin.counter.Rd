\name{bin.counter}
\Rdversion{0.1}
\alias{bin.counter}
\title{Bin Counting}
\description{Fast bin counts of spike times}
\usage{ 
bin.counter(x, b)
} 
\arguments{
\item{x}{spike times}
\item{b}{break points defining time bins. Must be an ordered vector with no duplications. Allowed to not cover the entire span of spike times}
}
\value{
Returns a vector giving the bin counts.
}

\examples{
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
}
\keyword{programming}
