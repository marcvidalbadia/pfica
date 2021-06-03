\name{pspline.kffobi}
\alias{pspline.kffobi}
\title{P-Spline smoothed functional ICA}
\description{
This function provides an alternative form of computing the smoothed functional ICA in terms of principal components (function \code{\link{kffobi}}). A discrete penalty  that measures the roughness of principal factors by summing squared \emph{r}-order differences between adjacent B-spline coefficients (P-spline penalty) is used; see Aguilera and Aguilera-Morillo (2013) for a detailed discussion.}
\usage{
pspline.kffobi(fdx, ncomp = fdx$basis$nbasis, pp = 0, r = 2,
               pr = c("fdx", "fdx.st", "KL", "KL.st"),
               shrinkage = FALSE, center = FALSE, plotfd = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{pp}{the penalty parameter. It can be estimated by \emph{leave-one-out} cross-validation or \emph{baseline} cross-validation (see \code{\link{bcv}}).}
  \item{r}{a number indicating the order of the penalty.}
  \item{pr}{the functional data object to project  into the space spanned by the eigenfunctions of the kurtosis operator.  To compute the independent components, the usual procedure is to use \code{KL.st}, the standardized principal component expansion. Thus, if \code{pr} is not supplied, then \code{KL.st} is used.}
  \item{shrinkage}{uses shrinkage estimators to compute the covariance matrix of the coordinate vectors.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
  \item{plotfd}{a logical value indicating whether to plot the eigenfunctions \cr of the kurtosis operator.}
}
\details{To compute the penalty matrix, the following code is used: \eqn{\Delta^{2}=\mathtt{diff(diag(nknots+2), differences = 2)}}, where \code{nknots} is the  number of basis knots. As in \code{\link{kffobi}}, the functional ICA of the principal component expansion is equivalent to the multivariate ICA of the principal coordinate vectors; see \emph{Details} in \code{\link{kffobi}}.}
\value{
a list with the following named entries:
\item{PCA.eigv}{a numeric vector giving the eigenvalues of the covariance operator.}
\item{PCA.basis}{a functional data object for the eigenfunctions of the covariance operator.}
\item{PCA.scores}{a matrix whose column vectors are the principal components.}
\item{ICA.eigv}{a numeric vector giving the eigenvalues of the kurtosis operator.}
\item{ICA.basis}{a functional data object for the eigenfunctions of the covariance operator.}
\item{ICA.scores}{a matrix whose column vectors are the projection coefficients for \code{fdx, fdx.st, KL} or \code{KL.st}.}
\item{ICA.kurtosis}{a numeric vector giving the kurtosis of each component vector.}
}
\references{
Aguilera, AM. and MC. Aguilera-Morillo (2013). “Penalized PCA approaches for B-spline expansions of smooth functional data”. In: \emph{Applied Mathematics and Computation} 219(14), pp. 7805–7819.

Vidal, M., M. Rosso and AM. Aguilera. (2021). Bi-Smoothed Functional Independent Component Analysis for EEG Artifact Removal. Mathematics 9(11) 1243.
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
ica.fd <- pspline.kffobi(x, 16, pp = 10)
## Plot by region in order of maximum kurtosis (outliers)
sc <- ica.fd$ICA.scores
plot(sc[,1], sc[,2], ylab = "", xlab = "")
text(sc[,1], sc[,2], CanadianWeather$region, pch=0.5, cex=0.6)}
\keyword{functional ICA}

