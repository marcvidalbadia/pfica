ffobi <- function(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
                  center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (center) {
    fdx <- center.fd(fdx)}
  a <- fdx$coefs
  phi <- fdx$basis
  rphi <- eigenfPar$fd$basis
  Lfdobj <- eigenfPar$Lfd
  lambda <- eigenfPar$lambda
  svda <- La.svd(tcrossprod(a))
  wa <- svda$u %*% diag(c(1/sqrt(svda$d))) %*% t(svda$u)
  wa <- sqrt(ncol(a)-1)*wa
  ast <- wa %*% a
  L <- eval.penalty(rphi, 0)
  if (lambda > 0) {
    R <- eval.penalty(rphi, Lfdobj)
    L <- L + lambda * R
  }; L <- (L + t(L))/2
  W <- solve(chol(L))
  J <- inprod(rphi, phi)
  rGram <- crossprod(W, J)
  nr <- sqrt(rowSums(ast^2))
  ast <- nr * ast
  kurt <- tcrossprod(ast)/(ncol(ast) * (ncomp + 2))
  C4 <- rGram %*% kurt %*% t(rGram); C4  <- (C4 + t(C4))/2
  svdk <- La.svd(C4)
  u <- as.matrix(svdk$u[, 1:ncomp])
  b <- W %*% u
  h <- fd(b, rphi)
  xst <- fd(ast, phi)
  zi <- inprod(xst, h)
  kurt <- vector()
  for (i in 1:ncomp) {
    K <- kurtosis(zi[,i])
    kurt[i]=K}
  if (plotfd) {
    oldpar <- par(mfrow=c(2,3), no.readonly = TRUE)
    on.exit(par(oldpar))
    for (j in 1:ncomp){
      txt <- paste('IC', j, "Kurt:", round(kurt[j],3))
      plot(h[j], col="#6495ED", bty="n", main=txt)
      par(ask=T)}; par(ask=F)}
  colnames(h$coefs) <- paste("eigenf.", c(1:ncomp), sep = "-")
  rownames(h$coefs) <- h$basis$names
  FICA <- list(h, kurt, zi)
  names(FICA) <- c("eigenbasis", "kurtosis", "scores")
  return(FICA)
}
