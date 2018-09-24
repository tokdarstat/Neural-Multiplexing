\name{dapp.simulate}
\Rdversion{0.1}
\alias{dapp.simulate}
\title{Simulate from Dynamic Admixture of Poisson Process}
\description{Simulate spike trains from DAPP model to binned spiking data}
\usage{ 
dapp.simulate(horizon = 1000, bin.width = 25, lengthScale,
     lsPrior = rep(1/length(lengthScale),length(lengthScale)),
     hyper = list(prec = c(1,1), sig0 = 1.87), nsamp = 1e3)
} 
\arguments{
\item{horizon}{time horizon of the response period (in ms)}
\item{bin.width}{width of the time bins (in ms) to be used to aggregate spike counts}
\item{lengthScale}{an array giving the length scale parameter values to be used for Gaussian process prior. Defaults to \code{sort(0.16 * resp.horiz / c(4, 3, 2, 1, 0.5, 0.1))} where \code{resp.horiz} is the time horizon of the response period.}
\item{lsPrior}{an array of the same length as \code{lengthScale} giving the prior probabilities of the length scale values.}
\item{hyper}{a list of hyper parameters with the following iterms. 'prec': a 2-vector giving the shape and rate parameters of the gamma distribution on the Dirichlet precision parameter. 'sig0': a scalaer giving the scale of the (centered) logistic distribution used in transforming the Gaussian random curves into curves restricted between 0 and 1.}
\item{nsamp}{number of priors draws to be made}
}
\value{
Returns a list of class "dapp" containting the following items.

\item{lsProb}{draws of length scale}
\item{alpha.pred}{prior predictive draws of alpha}
\item{prec}{draws of precision}
}

\details{Primarily intended to be used internally by the \code{\link{summary.dapp}} and \code{\link{plot.dapp}} functions. Could also be use to draw directly from the model.}

\examples{
\dontrun{
prior <- dapp.simulate(1000, 25)
}
}
\keyword{programming}
