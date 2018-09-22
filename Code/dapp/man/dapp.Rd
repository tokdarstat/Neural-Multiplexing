\name{dapp}
\Rdversion{0.1}
\alias{dapp}
\title{Dynamic Admixture of Poisson Process}
\description{Fits the DAPP model to binned spiking data}
\usage{ 
dapp(spike.counts, lengthScale = NULL, lsPrior = NULL, hyper = list(prec = c(1,1), sig0 = 1.87), burnIn = 1e3, nsamp = 1e3, thin = 4, plot = FALSE, verbose = TRUE, remove.zeros = FALSE)
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

\item{lsProb}{}
\item{lambda.A}}
\item{lambda.B}{}
\item{alpha}{}
\item{A}{}
\item{prec}{}
\item{alpha.pred}{}
\item{psl.pred}{}
\item{details}{}
\item{hyper}{}
\item{lengthScale}{}
\item{lsPrior}{}
\item{bin.mids}{}
\item{bin.width}{}
\item{mcmc}{}
}

\details{}


\examples{
\dontrun{
a <- 1:3
}
}
\keyword{programming}
