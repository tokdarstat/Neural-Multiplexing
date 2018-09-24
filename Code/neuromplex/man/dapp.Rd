\name{dapp}
\Rdversion{0.1}
\alias{dapp}
\title{Dynamic Admixture of Poisson Process}
\description{Fits the DAPP model to binned spiking data}
\usage{ 
dapp(spike.counts, lengthScale = NULL, lsPrior = NULL,
     hyper = list(prec = c(1,1), sig0 = 1.87), burnIn = 1e3,
     nsamp = 1e3, thin = 4, plot = FALSE, verbose = TRUE,
     remove.zeros = FALSE)
} 
\arguments{
\item{spike.counts}{A list with the following items. 'Acounts': binned spike counts under condition A presented as a matrix. Rows are bins, columns are replicates (trials). 'Bcount': binned spike counts under condition B. 'ABcounts': binned spike counts under condition AB. 'bin.mids': an array giving the mid-points of the time bins. 'bin.width': a scalar giving the bin width.}
\item{lengthScale}{an array giving the length scale parameter values to be used for Gaussian process prior. Defaults to \code{sort(0.16 * resp.horiz / c(4, 3, 2, 1, 0.5, 0.1))} where \code{resp.horiz} is the time horizon of the response period.}
\item{lsPrior}{an array of the same length as \code{lengthScale} giving the prior probabilities of the length scale values.}
\item{hyper}{a list of hyper parameters with the following iterms. 'prec': a 2-vector giving the shape and rate parameters of the gamma distribution on the Dirichlet precision parameter. 'sig0': a scalaer giving the scale of the (centered) logistic distribution used in transforming the Gaussian random curves into curves restricted between 0 and 1.}
\item{burnIn}{number of MCMC iterations to discard as burn-in.}
\item{nsamp}{number of MCMC draws to be saved for posterior inference.}
\item{thin}{the thinning rate at which MCMC draws are to be saved. The total number of iterations equals \code{burnIn + nsamp * thin}}
\item{plot}{logical indicating if a graphical update should be plotted during the course of the MCMC}
\item{verbose}{logical indicating if some fit details should be printed during the course of the MCMC}
\item{remove.zeros}{logical indicating if trials with zero spike count shuold be removed from the analysis}
}
\value{
Returns a list of class "dapp" containting the following items.

\item{lsProb}{posterior preditctive draws of length scale}
\item{lambda.A}{posterior draws of lambda.A at bin mid-points}
\item{lambda.B}{posterior draws of lambda.B at bin mid-points}
\item{alpha}{posterior draws of the alpha curves at bin mid-points}
\item{A}{posterior draws of the latent variable A which gives the AB spike counts (by bin) that are to be attributed to signal A (the remaining are attributed to signal B)}
\item{prec}{posterior draws of precision}
\item{alpha.pred}{posterior predictive draws of alpha (of a future trial)}
\item{psl.pred}{posterior predictive draw of the feature parameters (phi, psi, ell) (of a future trial)}
\item{details}{mcmc details given as an array of \code{c(niter, nsamp, burnIn, thin, MH acceptance rate)}}
\item{hyper}{hyper parameters used in model fitting}
\item{lengthScale}{length scale set used in model fitting}
\item{lsPrior}{length scale prior}
\item{bin.mids}{bin mid-points}
\item{bin.width}{bin width}
\item{mcmc}{mcmc controls (burn-in length, thinning rate and number of saved draws)}
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
