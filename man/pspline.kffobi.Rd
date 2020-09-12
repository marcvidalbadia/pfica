\name{pspline.kffobi}
\alias{pspline.kffobi}
\title{P-Spline smoothed functional ICA}
\description{
This function provides an alternative form of computing the smoothed functional ICA in terms of principal components (function \code{\link{kffobi}}), using a discrete penalty  that measures the roughness of principal factors by summing squared \emph{r}-order differences between adjacent B-spline coefficients (P-spline penalty); see Aguilera and Aguilera-Morillo (2013) for a detailed discussion.}
\usage{
pspline.kffobi(fdx, ncomp = fdx$basis$nbasis, rho = 0, r = 2,
               pr = c("fdx", "fdx.st", "KL", "KL.st"),
               center = FALSE, plotfd = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdx}{a functional data object obtained from the \strong{fda} package.}
  \item{ncomp}{number of independent components to compute.}
  \item{rho}{the smoothing parameter. It can be estimated by \emph{leave-one-out} cross-validation.}
  \item{r}{a number indicating the order of the penalty.}
  \item{pr}{the functional data object to project  into the space spanned by the eigenfunctions of the FOBI operator.  To compute the independent components, the usual procedure is to use \code{KL.st}, the standardized principal component expansion. Thus, if \code{pr} is not supplied, then \code{KL.st} is used.}
  \item{center}{a logical value indicating whether the mean function has to be subtracted from each functional observation.}
  \item{plotfd}{a logical value indicating whether to plot the eigenfunctions \cr of the FOBI operator.}
}
\details{The smooth FPCA can be summarized in the following generalized matrix eigenproblem: \deqn{n^{-1}Vc=\lambda(G+\rho G_{r})c,} with \eqn{V} as the covariance matrix of \eqn{A}; see \emph{Details} in \code{\link{kffobi}}. Then, the matrix of penalties \eqn{G_{r}} is defined as \eqn{G_{r}=\left(\Delta^{r}\right)^{\prime} \Delta^{r}} with \eqn{\Delta^{r}} as the matrix representation of the \emph{r}-order difference operator. \cr For \eqn{r=2} this matrix is given by \deqn{\Delta^{2}=\left[\begin{array}{ccccc}
1 & -2 & 1 & 0 & \cdots \\
0 & 1 & -2 & 1 & \cdots \\
0 & 0 & 1 & -2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \ddots
\end{array}\right].} In \strong{R}, this can be translated into the following code: \eqn{\Delta^{2}=\mathtt{diff(diag(nknots+2), differences = 2)}}, where \code{nknots} is the  number of basis knots.

Once the P-spline approximation of eigenfunctions has been computed and the smooth principal component expansion of P-splines is obtained, then the ICA procedure is conducted. As in \code{\link{kffobi}}, the functional ICA of the principal component expansion is equivalent to the ICA of the matrix \eqn{Z}; see \emph{Details} in \code{\link{kffobi}}.}
\value{
a list with the following named entries:
\item{eigenbasis}{a functional data object for the eigenfunctions or independent factors.}
\item{kurtosis}{a numeric vector giving the kurtosis associated to each independent component vector.}
\item{scores}{a matrix whose column vectors are the independent components.}}
\references{
Aguilera, AM. and MC. Aguilera-Morillo (2013). “Penalized PCA approaches for B-spline expansions of smooth functional data”. In: \emph{Applied Mathematics and Computation} 219(14), pp. 7805–7819.

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
ica.fd <- pspline.kffobi(x, 16, rho = 10)
## Plot by region: classification in order of maximum kurtosis
k1 <- which.max(ica.fd$kurtosis)
k2 <- which.max(ica.fd$kurtosis[c(-k1)])+1
sc <- ica.fd$scores
plot(sc[,c(k1)], sc[,c(k2)], ylab = "", xlab = "")
text(sc[,c(k1)], sc[,c(k2)], CanadianWeather$region, pch=0.5, cex=0.6)}
\keyword{functional ICA}

