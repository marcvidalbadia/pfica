\name{ffobi}
\alias{ffobi}
\title{Smoothed functional ICA in terms of basis functions coordinates}
\description{
This function computes the ordinary ICA procedure from a sample represented by basis functions (Fourier, B-splines...). The estimation method is based on the use of fourth moments (FOBI), in which it is assumed that the independent components have different kurtosis values. The proposed algorithm can be considered an extension of the implementation of the kurtosis operator introduced in Peña et. al (2014), whose decomposition is used to identify cluster structures and outliers.}
\usage{
ffobi(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
      center = FALSE, plotfd = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{eigenfPar}{a functional parameter object, obtained from the \strong{fda} package, that defines the independent component functions to be estimated.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
  \item{plotfd}{a logical value indicating whether to plot the eigenfunctions \cr of the FOBI operator.}
}
\details{This IC model for functional data consists in performing the multivariate ICA of a transformation of the coordinate vectors associated to a basis of functions. This algorithm also incorporates a continuous penalty in the orthonormality constraint.

Let \eqn{x=(x_{i}(t)\colon i=1,\ldots,n;t\in T)} be a set of functions in a finite-dimensional space spanned by a basis of functions \eqn{\{\phi_{1}(t),\ldots,\phi_{p}(t)\}}. Assume that each function of \eqn{x} and an independent factor \eqn{h} can be expanded as \deqn{x_{i}(t)=\sum_{j=1}^{p}a_{ij}\phi_{j}(t),\quad h(t)=\sum_{j=1}^{p}c_{j}\phi_{j}(t)=\phi(t)^{\prime}c.} Consider also the matrix  \eqn{G=(\int_{T}\phi_{j}(t)\phi_{j'}(t)dt;j,j'=1,\ldots,p)} and the continuous matrix of penalties \eqn{G_{r}=\int_{T}\left(\mathcal{D}^{r}\phi_{j}(t)\right)\left(\mathcal{D}^{r}\phi_{j'}(t)\right)'dt,} where \eqn{\mathcal{D}} denotes \emph{r}-order difference operator. \cr Then, the followed ICA procedure can be summarized in the following steps:

\enumerate{
   \item Standardization of the coordinates: \eqn{A_{st}=\mathrm{cov}(A)^{-\frac{1}{2}}\left(A-\mathrm{E}\left[A\right]\right),} where \eqn{A=(a_{ij})_{n\times p}}.
   \item Find solutions to the matrix generalized eigenproblem: \deqn{n^{-1}K^{\bullet}c=\lambda(G+\rho G_{r})c} with \eqn{K^{\bullet}} a fixed FOBI matrix of \eqn{A_{st}G} and \eqn{\rho} a smoothing parameter.
   \item Obtain the matrix of coordinates \eqn{C_{p\times p},} whose columns are the eigenvectors of the FOBI matrix \eqn{K^{\bullet},} such that \eqn{A_{st}C} are the independent components of \eqn{x}.}
Note that, by applying Cholesky factorization of the form \eqn{LL'=G+\rho G_{r},} the generalized eigenproblem in step 2 leads to an eigenproblem  of a symmetric FOBI matrix of \eqn{A_{st}GL^{-1'}}. Then, the basis coefficients of the independent factors are given by \eqn{c_{j}=L^{-1^{\prime}} v_{j}}, where \eqn{v_j} are the eigenvectors of the above FOBI matrix.}
\value{a list with the following named entries:
\item{eigenbasis}{a functional data object for the eigenfunctions or independent factors.}
\item{kurtosis}{a numeric vector giving the kurtosis associated to each independent component vector.}
\item{scores}{a matrix whose column vectors are the independent components.}}
\references{
Aguilera, AM. and MC. Aguilera-Morillo (2013). “Penalized PCA approaches for B-spline expansions of smooth functional data”. In: \emph{Applied Mathematics and Computation} 219(14), pp. 7805–7819.

Peña, C., J. Prieto, and C. Rendón (2014). \emph{Independent components techniques based
on kurtosis for functional data analysis}. Working paper 14-10. Universidad Carlos
III de Madrid.

Ramsay, J. and B. Silverman (2005). \emph{Functional Data Analysis}. Springer.

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
## Plot by region: classification in order of maximum kurtosis
k1 <- which.max(ica.fd$kurtosis)
k2 <- which.max(ica.fd$kurtosis[c(-k1)])+1
sc <- ica.fd$scores
plot(sc[,c(k1)], sc[,c(k2)], ylab = "", xlab = "")
text(sc[,c(k1)], sc[,c(k2)], CanadianWeather$region, pch=0.5, cex=0.6)
}
\keyword{functional ICA}
