\name{ffobi}
\alias{ffobi}
\title{Smoothed functional ICA}
\description{
This function computes the ordinary ICA procedure from a sample represented by basis functions (Fourier, B-splines...). The estimation method is based on the use of fourth moments (kurtosis), in which it is assumed that the independent components have different kurtosis values. The proposed algorithm is an implementation of the kurtosis operator introduced in Peña et. al (2014).}
\usage{
ffobi(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
      w = c("PCA", "PCA-cor","ZCA", "ZCA-cor",
            "Varimax", "Varimax-cor", "Cholesky"),
      pr = c("fdx", "wfdx"), center = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{eigenfPar}{a functional parameter object obtained from the \strong{fda} package.}
  \item{w}{the whitening procedure. By default \code{ZCA} (Mahalanobis whitening) is used.}
    \item{pr}{the functional data object to project on to the space spanned by the eigenfunctions of the kurtosis kernel function.  To compute the independent components, the usual procedure is to use \code{wfdx}, the whitened basis expansion. Thus, if \code{pr} is not supplied, \code{wfdx} is used.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
}
\details{This functional ICA consists of performing the multivariate ICA of a transformation of the coordinate vectors associated with a basis of functions. The algorithm also incorporates a continuous penalty in the orthonormality constraint of the kurtosis eigenfunctions.}
\value{a list with the following named entries:
\item{ICA.eigv}{a numeric vector giving the eigenvalues of the kurtosis kernel function.}
\item{ICA.basis}{a functional data object for the kurtosis kernel eigenfunctions.}
\item{ICA.scores}{a matrix whose column vectors are the projection coefficients for \code{fdx} or \code{wfdx}.}
}
\references{
Peña, C., J. Prieto and C. Rendón (2014). \emph{Independent components techniques based on kurtosis for functional data analysis}. Working paper 14–10 Statistics and Econometric Series (06); Universidad Carlos III de Madrid, Madrid, 2014.}
\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{kffobi}}, \code{\link{whiten.fd}}}
\examples{
## Canadian Weather data
library(fda)
arg <- 1:365
Temp <- CanadianWeather$dailyAv[,,1]
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=16)
x <- Data2fd(Temp, argvals = arg, B)
Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf <- fdPar(B, Lfdobj, lambda=10^4) #lambda
ica.fd <- ffobi(x, 16, penf)
}
\keyword{functional ICA}
