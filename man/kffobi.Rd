\name{kffobi}
\alias{kffobi}
\title{Smoothed functional PCA/ICA}
\description{
This function computes the ordinary ICA procedure on a smoothed principal component expansion (also known as Karhunen-Loeve expansion) whose eigenbasis is expressed in terms of basis functions (Fourier, B-splines...). The estimation method is based on the use of fourth moments (kurtosis), in which it is assumed that the independent components have different kurtosis values.  The function further provides smoothed (and non-smoothed) functional PCA estimates based on Ocaña et al. (1999).}
\usage{
kffobi(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
       w = c("PCA", "PCA-cor","ZCA", "ZCA-cor","Cholesky"),
       pr = c("fdx", "wfdx", "KL", "wKL"),
       center = TRUE)}
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute (must be > 1).}
  \item{eigenfPar}{a functional parameter object obtained from the \strong{fda} package.}
  \item{w}{the whitening procedure. By default \code{ZCA} (Mahalanobis whitening) is used.}
  \item{pr}{the functional data object to project  on to the space spanned by the eigenfunctions of the kurtosis kernel function.  To compute the independent components, the usual procedure is to use \code{wKL}, the whitened principal component expansion. If \code{pr} is not supplied, then \code{wKL} is used. Note that \code{w} stands for whitened, e.g., \code{wfdx} is the whitened original data.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
}
\details{This functional ICA consists of performing the multivariate ICA of the PC coordinate vectors in terms of the principal component weight functions.}

\value{a list with the following named entries:
\item{PCA.eigv}{a numeric vector giving the eigenvalues of the covariance kernel function.}
\item{PCA.basis}{a functional data object for the eigenfunctions of the covariance kernel function.}
\item{PCA.scores}{a matrix whose column vectors are the principal components.}
\item{ICA.eigv}{a numeric vector giving the eigenvalues of the kurtosis kernel function.}
\item{ICA.eigv}{a numeric vector giving the eigenvalues of the kurtosis kernel function.}
\item{ICA.scores}{a matrix whose column vectors are the projection coefficients for \code{fdx, wfdx, KL} or \code{wKL}.}
\item{wKL}{the whitened principal components expansion whith coefficients in terms of basis functions.}
}
\references{
Li, B., G. Van Bever, H. Oja, R. Sabolova and F. Critchley (2019). \emph{Functional independent component analysis: an extension of the fourth-orderblind identification}. Technical Report, Univ. Namur.

Ocaña, F.A., A.M. Aguilera and M.J. Valderrama (1999). \emph{Functional Principal Components Analysis by Choice of Norm}. Journal of Multivariate Analysis, 71(2), <doi:10.1006/jmva.1999.1844>.
}
\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{ffobi}}, \code{\link{pspline.kffobi}}, \code{\link{whiten.fd}}}
\examples{
## Canadian Weather data
library(fda)
arg <- 1:365
Temp <- CanadianWeather$dailyAv[,,1]
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=14)
x <- Data2fd(Temp, argvals = arg, B)
#plot(x) #plot data
Lfdobj <- int2Lfd(max(0, norder(B)-2))
penf <- fdPar(B, Lfdobj, lambda=10^9) #penalty parameter (lambda)
ica.fd <- kffobi(x, 2, eigenfPar = penf, w="ZCA-cor")

## Whitened data using the two first smoothed principal components
wKL <- ica.fd$wKL

## Plot by region
sc <- ica.fd$ICA.scores
plot(sc[,1], sc[,2], pch = 20, col = factor(CanadianWeather$region),
     ylab = "IC 1", xlab = "IC 2")
}

\keyword{functional ICA}
