\name{whiten.fd}
\alias{whiten.fd}
\title{Whitening (or sphering) functional data}
\description{
This function whitens functional data in terms of basis functions expansions.

\strong{Functional whitening procedures available:}
 \itemize{
    \item \code{PCA}: uses a projection that combines the covariance function eigenfactors and an arbitrary orthonormal basis
    \item \code{PCA-cor}: PCA whitening based on the correlation  function
    \item \code{ZCA}: or Mahalanobis whitening, is a symmetric solution based on the spectral decomposition of the inverse square root of the covariance function
    \item \code{ZCA-cor}: Mahalanobis whitening based on the correlation function
    \item \code{Varimax}: combines the covariance function eigenfactors and its varimax rotation
    \item \code{Varimax-cor}: a Varimax solution based on the correlation function
    \item \code{Cholesky}: uses the Cholesky factorization of the inverse covariance function.
  }
\strong{Which functional whitening procedure to use?} \code{PCA} and \code{PCA-cor} allow for maximum compression of the functional observations and therefore the new data representations might be rather noisy. \code{ZCA} and \code{ZCA-cor} maximize the similarity with the original data, providing most robust whitening. The \code{Varimax}/\code{Varimax-cor} procedures, which are analogous to \code{ZCA}/\code{ZCA-cor}, can be used when the first component of the functional PCA accounts for a high percentage of the total variance. \code{Cholesky} takes an intermediate position in scoring for maximum compression and similarity indices. To explore low dimensional structures in the data, \code{\link{kffobi}} and \code{\link{pspline.kffobi}} allows to compute these functional whitening representations using a few principal components.

For further details, see Vidal and Aguilera (2022).
}
\usage{
whiten.fd(fdx, w = c("PCA", "PCA-cor","ZCA", "ZCA-cor",
                     "Varimax", "Varimax-cor", "Cholesky"),
          pre.center = TRUE, post.center = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{w}{the whitening procedure. By default \code{ZCA} (Mahalanobis whitening) is used.}
  \item{pre.center}{centering the functional observations before performing the transformation.}
  \item{post.center}{centering the functional observations after performing the transformation.}
}
\details{As in the \strong{whitening} package (see Kessy et al. (2018)), to solve the sign ambiguity in  \code{PCA}, \code{PCA-cor} we use eigenvector matrices with a positive diagonal. The \code{pre.centrer}/\code{post.centrer} options can be used interchangeably.}
\value{
\item{wfdx}{a functional data object with whitened coefficients.}
}
\references{
Acal, C., A.M.  Aguilera and M. Escabias (2020). \emph{New Modeling Approaches Based on Varimax Rotation of Functional Principal Components}. Mathematics, 8, issue 11, p. 1-15. <doi:10.3390/math8112085>

Kessy, A., A. Lewin and K. Strimmer (2018). \emph{Optimal Whitening and Decorrelation}. The American Statistician, 72:4, 309-314, <doi:10.1080/00031305.2016.1277159>.

Vidal, M. and and A.M. Aguilera (2022). \emph{Novel whitening approaches in functional settings}.  Stat, e516,
<doi:10.1002/sta4.516>.
}

\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{ffobi}}}
\examples{
## Canadian Weather data
library(fda)
arg <- 1:365
Temp <- CanadianWeather$dailyAv[,,1]
B <- create.bspline.basis(rangeval=c(min(arg),max(arg)), nbasis=16)
x <- Data2fd(Temp, argvals = arg, B)
wx <- whiten.fd(x)
}
\keyword{functional whitening}
