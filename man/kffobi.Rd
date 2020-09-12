\name{kffobi}
\alias{kffobi}
\title{Smoothed functional ICA in terms of principal components}
\description{
This function computes the ordinary ICA procedure from a penalized principal component expansion (also known as Karhunen-Loève expansion) whose basis, the eigenfunctions of the covariance matrix operator, is expressed in terms of basis functions (Fourier, B-splines...). The estimation method is based on the use of fourth moments (FOBI), in which it is assumed that the independent components have different kurtosis values. The proposed algorithm can be considered an extension of the IC model proposed in Li et al. (2015). This function provides more accurate estimates than \link{ffobi} and was used in Vidal (2020) to identify artifactual independent curves in bioelectrical signals.}
\usage{
kffobi(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
       pr = c("fdx", "fdx.st", "KL", "KL.st"),
       center = FALSE, plotfd = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{eigenfPar}{a functional parameter object, obtained from the \strong{fda} package, that defines the principal component functions to be estimated.}
  \item{pr}{the functional data object to project  into the space spanned by the eigenfunctions of the FOBI operator.  To compute the independent components, the usual procedure is to use \code{KL.st}, the standardized principal component expansion. Thus, if \code{pr} is not supplied, \code{KL.st} is used.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
  \item{plotfd}{a logical value indicating whether to plot the eigenfunctions \cr of the FOBI operator.}}
\details{Note that \code{kffobi} first computes the (penalized) functional PCA; see Aguilera and Aguilera-Morillo (2013) for a detailed discussion. Thus here,  the IC model for functional data consists in performing the multivariate ICA of a transformation of the principal component coordinate vectors.

Let \eqn{x=(x_{i}(t)\colon i=1,\ldots,n;t\in T)} be a set of functions in a finite-dimensional space spanned by a basis of functions \eqn{\{\phi_{1}(t),\ldots,\phi_{p}(t)\}}. Assume that the principal component expansion of \eqn{x} is  given by \deqn{x_{i}^{(q)}\left(t\right)=\sum_{j=1}^{q}z_{ij}f_{j}\left(t\right)} where \eqn{(z_{ij})_{n\times q}=Z} is the matrix of principal components and \eqn{f_j} the principal factors of \eqn{x}. Notice that \eqn{f_j} can be expressed in terms of basis functions, i.e. \eqn{f_{j}(t)=\phi(t)'b_{j}}, where the coordinates \eqn{b_j} are computed from the FPCA of \eqn{x}. The ICA procedure consists of the same steps as the ones described in \code{\link{ffobi}}, but replacing the matrix \eqn{A} by \eqn{Z}. As the  orthogonalization of the basis functions is implemented in the FPCA algorithm, now \eqn{\int_{T}f_{j}(t)f_{j'}(t)dt=I_{q};\forall j,j'=1,\ldots,q.}, where \eqn{I_{q}} denotes the identity matrix. As such, the ICA procedure from a principal component expansion is easy to compute.}
\value{a list with the following named entries:
\item{eigenbasis}{a functional data object for the eigenfunctions or independent factors.}
\item{kurtosis}{a numeric vector giving the kurtosis associated to each independent component vector.}
\item{scores}{a matrix whose column vectors are the independent components.}}
\references{
Aguilera, AM. and MC. Aguilera-Morillo (2013). “Penalized PCA approaches for B-spline expansions of smooth functional data”. In: \emph{Applied Mathematics and Computation} 219(14), pp. 7805–7819.

Li, B., G. Van Bever, H. Oja, R. Sabolová, and F. Critchley (2015). “Functional independent component analysis: an extension of the fourth-orderblind identification.” \emph{Submitted}.

Miettinen, J., K. Nordhausen, and S. Taskinen (2017). “Blind source separation based on joint diagonalization in R: The packages JADE and BSSasymp”. In: \emph{Journal of Statistical Software} 76.2, pp. 1–31.

Ramsay, J. and B. Silverman (2005). \emph{Functional Data Analysis}. Springer.

Vidal, M. (2020). \emph{Functional Independent Component Analysis in Bioelectrical Signal Processing}. MA thesis. Universidad de Granada.
}
\author{Marc Vidal, Ana Mª Aguilera}
\seealso{\code{\link{kd}}, \code{\link{ffobi}}}
\examples{
\donttest{
## foetal_ecg data
library(fda)
dataset <- matrix(
  scan("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v076i02/foetal_ecg.dat"),
  2500, 9, byrow = TRUE);
X <- dataset[1:1000, 2:9]
arg <- 1:1000
basis <- create.fourier.basis(rangeval=c(min(arg), max(arg)), nbasis=301, basisvalues=TRUE)
x <- Data2fd(X, argvals=arg, basis)
## Penalization can be considered:
#Lfdobj <- vec2Lfd(c(0,(2*pi/diff(4))^2,0), 4)
#hm <- fdPar(base, Lfdobj, lambda=2)
## Select the number of components with the kurtosis distance:
#kurt.dist <- kd(x, qmax = 8)
aci <- kffobi(x, 7, plotfd = TRUE)
}}

\keyword{functional ICA}
