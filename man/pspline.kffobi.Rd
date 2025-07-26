\name{pspline.kffobi}
\alias{pspline.kffobi}
\title{P-Spline smoothed functional PCA/ICA}
\description{
This function is an alternative way of computing the smoothed functional ICA in terms of principal components (see \code{\link{kffobi}}). The function further provides smoothed (and non-smoothed) functional PCA estimates. A discrete penalty  that measures the roughness of principal factors by summing squared \emph{r}-order differences between adjacent B-spline coefficients (P-spline penalty) is used; see Aguilera and Aguilera-Morillo (2013) and Eliers and Marx (2021).}
\usage{
pspline.kffobi(fdx, ncomp = fdx$basis$nbasis, pp = 0, r = 2,
               w = c("PCA", "PCA-cor","ZCA", "ZCA-cor","Cholesky"),
               pr = c("fdx", "wfdx", "KL", "wKL"),
               center = TRUE)}
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute (must be > 1).}
  \item{pp}{the penalty parameter. It can be estimated using \emph{leave-one-out} cross-validation.}
  \item{r}{a number indicating the order of the penalty.}
  \item{w}{the whitening procedure. By default \code{ZCA} (Mahalanobis whitening) is used.}
  \item{pr}{the functional data object to project on to the space spanned by the eigenfunctions of the kurtosis kernel function.  To compute the independent components, the usual procedure is to use \code{wKL}, the whitened principal component expansion. If \code{pr} is not supplied, then \code{wKL} is used. Note that \code{w} stands for whitened, e.g., \code{wfdx} is the whitened original data.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
}
\value{
a list with the following named entries:
\item{PCA.eigv}{a numeric vector giving the eigenvalues of the covariance kernel function.}
\item{PCA.basis}{a functional data object for the eigenfunctions of the covariance kernel function.}
\item{PCA.scores}{a matrix whose column vectors are the principal components.}
\item{ICA.eigv}{a numeric vector giving the eigenvalues of the kurtosis kernel function.}
\item{ICA.basis}{a functional data object for the kurtosis kernel eigenfunctions.}
\item{ICA.scores}{a matrix whose column vectors are the projection coefficients for \code{fdx, wfdx, KL} or \code{wKL}.}
\item{wKL}{the whitened principal components expansion whith coefficients in terms of basis functions.}
}
\references{
Aguilera, A.M. and M.C. Aguilera-Morillo (2013). \emph{Penalized PCA approaches for B-spline expansions of smooth functional data.} Applied Mathematics and Computation 219(14), 7805–7819, <doi:10.1016/j.amc.2013.02.009>.

Eilers, P. and B. Marx (2021). \emph{Practical Smoothing: The Joys of P-splines.} Cambridge: Cambridge University Press, <doi:10.1017/9781108610247>.

Vidal, M., M. Rosso and A.M. Aguilera. (2021). \emph{Bi-Smoothed Functional Independent Component Analysis for EEG Artifact Removal}. Mathematics 9(11), 1243, <doi:10.3390/math9111243>.
}
\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{kffobi}}, \code{\link{whiten.fd}}}
\examples{
## Canadian Weather data
library(fda)
arg <- 1:365
Temp <- CanadianWeather$dailyAv[,,1]
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=14)
x <- Data2fd(Temp, argvals = arg, B)
#plot(x) #plot data
ica.fd <- pspline.kffobi(x, 2, pp = 10^4, w="ZCA-cor")

## Whitened data using the two first smoothed principal components
wKL <- ica.fd$wKL

## Plot by region
sc <- ica.fd$ICA.scores
plot(sc[,2], sc[,1], pch = 20, col = factor(CanadianWeather$region),
     ylab = "IC 1", xlab = "IC 2")
}
\keyword{functional ICA}

