\name{ffobi}
\alias{ffobi}
\title{Smoothed functional ICA in terms of basis functions coordinates}
\description{
This function computes the ordinary ICA procedure from a sample represented by basis functions (Fourier, B-splines...). The estimation method is based on the use of fourth moments (kurtosis), in which it is assumed that the independent components have different kurtosis values. The proposed algorithm can be considered an extension of the implementation of the kurtosis operator introduced in Peña et. al (2014), whose decomposition is used to identify cluster structures and outliers.}
\usage{
ffobi(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
      pr = c("fdx", "fdx.st"), shrinkage = FALSE,
      center = FALSE, plotfd = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{eigenfPar}{a functional parameter object, obtained from the \strong{fda} package, that defines the independent component functions to be estimated.}
    \item{pr}{the functional data object to project  into the space spanned by the eigenfunctions of the kurtosis operator.  To compute the independent components, the usual procedure is to use \code{fdx.st}, the standardized basis expansion. Thus, if \code{pr} is not supplied, \code{fdx.st} is used.}
  \item{shrinkage}{uses shrinkage estimators to compute the covariance matrix of the coordinate vectors.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
  \item{plotfd}{a logical value indicating whether to plot the eigenfunctions \cr of the kurtosis operator.}
}
\details{This IC model for functional data consists in performing the multivariate ICA of a transformation of the coordinate vectors associated to a basis of functions. The algorithm also incorporates a continuous penalty in the orthonormality constraint.}
\value{a list with the following named entries:
\item{eigenbasis}{a functional data object for the eigenfunctions or independent factors.}
\item{kurtosis}{a numeric vector giving the proportion of variance kurtosis explained by each eigenfunction.}
\item{scores}{a matrix whose column vectors are the independent components.}}
\references{
Peña, C., J. Prieto, and C. Rendón (2014). \emph{Independent components techniques based
on kurtosis for functional data analysis}. Working paper 14-10. Universidad Carlos
III de Madrid.

Ramsay, J. and B. Silverman (2005). \emph{Functional Data Analysis}. Springer.

Schafer, J. and K. Strimmer (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical Applications in Genetics and Molecular Biology, 4.32(1).

Vidal, M. (2020). \emph{Functional Independent Component Analysis in Bioelectrical Signal Processing}. MA thesis. Universidad de Granada.
}
\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{kffobi}}}
\examples{
## Canadian Weather data
library(fda)
arg <- 1:365
Temp <- CanadianWeather$dailyAv[,,1]
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=16)
x <- Data2fd(Temp, argvals = arg, B)
Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf <- fdPar(B, Lfdobj, lambda=10^4)
ica.fd <- ffobi(x, 16, penf)

## Plot by region in order of maximum kurtosis (outliers)
sc <- ica.fd$scores
plot(sc[,1], sc[,2], ylab = "", xlab = "")
text(sc[,1], sc[,2], CanadianWeather$region, pch=0.5, cex=0.6)
}
\keyword{functional ICA}
