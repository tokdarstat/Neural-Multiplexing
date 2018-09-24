\name{poisson.tests}
\Rdversion{0.1}
\alias{poisson.tests}
\title{Poisson Tests for Whole Trial Spike Counts}
\description{Carries out various Poisson related tests for double-stimuli spike count distribution.}
\usage{ 
poisson.tests(xA, xB, xAB, labels = c("A", "B", "AB"),
     remove.zeros = FALSE, plot = FALSE,
     gamma.pars = c(0.5, 2e-10), beta.pars = c(0.5, 0.5))
} 
\arguments{
\item{xA}{an array of whole-trial spike counts under stimulus 1}
\item{xB}{an array of whole-trial spike counts under stimulus 2}
\item{xAB}{an array of whole-trial spike counts when both stimuli are present together}
\item{labels}{labels for stimlus conditions}
\item{remove.zeros}{whether to remove trials with zero spike counts}
\item{plot}{logical indicating if a visualization plot should be made}
\item{gamma.pars}{shape and rate parameters of the gamma prior on Poisson mean}
\item{beta.pars}{shape parameters of the beta prior for the mixture/intermediate parameter}
}
\value{
Returns a list with the following items:
\item{separation.logBF}{the (log) Bayes factor for testing that that two single stimulus distributions are different}
\item{post.prob}{posterior probabilities of the four hypotheses (Mixture, Intermediate, Outside, Single) under equal prior probabilities}
\item{pois.pvalue}{minimum of the two p-values checking for Poisson-ness of each single stimulus distribution}
\item{sample.sizes}{three trial counts for A, B and AB conditions}
}

\details{To be added...}
\examples{
\dontrun{
nA <- 20; nB <- 15; nAB <- 25
muA <- 25; muB <- 40
Acounts <- rpois(nA, muA)
Bcounts <- rpois(nB, muB)
ABcounts <- rpois(nAB, sample(c(muA, muB), nAB, replace = TRUE))
poisson.tests(Acounts, Bcounts, ABcounts)
}
}
\keyword{programming}
