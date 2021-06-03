\name{bcv}
\alias{bcv}
\title{Compute baseline cross-validation criterion (P-spline version)}
\description{The generalized cross-validation method is adapted for B-spline representations with large number of knots to select the penalty parameter in \code{\link{pspline.kffobi}}. This function computes the log cross-validation score for a given penalty parameter.}
\usage{ bcv(a, phi, G, r = 2, ncomp, baseline, value = 0.1) }
\arguments{
  \item{a}{the coefficient matrix of a functional data object obtained from the \strong{fda} package.}
  \item{phi}{a functional basis object associated to the coefficient matrix \code{a}.}
  \item{G}{the Gramian matrix of inner products between pairs of basis elements in \code{phi}.}
  \item{r}{a number indicating the order of the penalty.}
  \item{ncomp}{number of principal components to compute.}
  \item{baseline}{a non-negative value (the penalty parameter).}
  \item{value}{a number that slightly  increases the \code{baseline} penalty parameter. This value may automatically change depending on the range of  \code{baseline}. The default is set to 0.1.}
}
%\details{ }
\value{
The baseline log cross-validation score.
}
\references{
Vidal, M., Rosso, M. and Aguilera, A. M. (2021). Bi-Smoothed Functional Independent Component Analysis for EEG Artifact Removal. Mathematics 9(11) 1243.
}
\author{ Marc Vidal, Ana Mª Aguilera}
%\note{}
%\examples{ }
\seealso{\code{\link{pspline.kffobi}}}
\keyword{utilities}

