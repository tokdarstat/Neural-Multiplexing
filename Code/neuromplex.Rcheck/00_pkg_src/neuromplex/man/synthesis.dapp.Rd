\name{synthesis.dapp}
\Rdversion{0.1}
\alias{synthesis.dapp}
\title{Simulate Multiplexing Data for DAPP Analysis}
\description{Simulate spike trains from controlled DAPP setting with flat and sinusoidal weight curves}
\usage{ 
synthesis.dapp(ntrials = c(10, 10, 10), time.bins = 0:1000, lambda.A = 400,
     lambda.B = 100, pr.flat = 0.5, intervals = list(c(0,1)),
     wts = 1, span = c(0,1), period.range = c(400, 1000))
} 
\arguments{
\item{ntrials}{an array with 3 elements giving the trial counts for conditions A, B and AB}
\item{time.bins}{time bins (in ms) giving the break points of the time bins in which Poisson draws should be made to mimic a Poisson process generation}
\item{lambda.A}{a flat intensity (in Hz) for condition A}
\item{lambda.B}{a flat intensity (in Hz) for condition B}
\item{pr.flat}{proportion of flat weight curves to be generated}
\item{intervals}{a list of sub-intervals (each represented by the 2-vector giving the sub-interval end-points) which determine the ranges of the flat weight curves}
\item{wts}{the relative weights of the sub-intervals above}
\item{span}{a two-vector giving the range of the sinusoidal weight curves}
\item{period.range}{the range from which the sinusoidal periods are drawn randomly (and uniformly)}
}
\value{
Returns a list containting the following items.

\item{spiketimes}{a list of spiketimes for each AB trial}
\item{alphas}{true underlying weight curves for each AB trial}
\item{lambdas}{corresponding intensity curves for each AB trial}
\item{time.pts}{time points associated with alphas and lambdas}
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
spike.counts$bin.mids <- breaks[-1] - mean(diff(breaks))/2
spike.counts$bin.width <- diff(breaks)[1]

}
\keyword{programming}
