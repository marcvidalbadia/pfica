ffobi <- function(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
                  pr = c("fdx", "fdx.st"),
                  shrinkage = FALSE, center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "fdx.st"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project")
  }
  if (center) {
    fdx <- center.fd(fdx)}

  Lfdobj <- eigenfPar$Lfd
  lambda <- eigenfPar$lambda
  phi <- fdx$basis
  rphi <- eigenfPar$fd$basis
  L <- eval.penalty(rphi, 0)
  if (lambda > 0) {
    R <- eval.penalty(rphi, Lfdobj)
    L <- L + lambda * R
  }; L <- (L + t(L))/2
  W <- solve(chol(L))
  J <- inprod(rphi, phi)
  rGram <- crossprod(W, J)

  a <- fdx$coefs
  if (shrinkage == TRUE) {
    svda <- La.svd(corpcor::cov.shrink(t(a), verbose = F))
  } else {
    svda <- La.svd(tcrossprod(a)/ncol(a)) }
  wa <-  (W %*% svda$u) %*% diag(c(1/sqrt(svda$d))) %*% t(W%*%svda$u)
  wa <- sqrt(ncol(a)-1) * wa
  asta <- wa %*% a
  nr <- numeric()
  for (i in 1:ncol(asta)) {
    nr[i] <- (t(asta[,i])%*%rGram %*% (rGram%*%asta[,i]))}
  ast <- nr * asta
  kurt <- tcrossprod(ast)/(ncol(ast) * (ncomp + 2))
  C4 <- rGram %*% kurt %*% t(rGram); C4  <- (C4 + t(C4))/2

  svdk <- La.svd(C4)
  u <- as.matrix(svdk$u[, 1:ncomp])
  b <- W %*% u
  h <- fd(b, rphi)
  svdk <- La.svd(C4)
  u <- as.matrix(svdk$u[, 1:ncomp])
  b <- W %*% u
  h <- fd(b, rphi)
  xst <- fd(asta, phi)
  project <- list(fdx, xst)
  names(project) <- c("fdx", "fdx.st")
  zi <- inprod(project[[paste(pr)]], h)
  kurt <- vector()
  for (i in 1:ncomp) {
    K <- moments::kurtosis(zi[,i])
    kurt[i]=K}

  if (plotfd) {
    oldpar <- par(mfrow=c(2,3), no.readonly = TRUE)
    on.exit(par(oldpar))
    for (j in 1:ncomp){
      plot(h[j])
      title(paste('IC', j, "Kurt:", round(kurt[j],3)))
      par(ask=T)}; par(ask=F)}
  colnames(h$coefs) <- paste("eigenf.", c(1:ncomp), sep = "-")
  rownames(h$coefs) <- h$basis$names
  FICA <- list(h, kurt, zi)
  names(FICA) <- c("eigenbasis", "kurtosis", "scores")
  return(FICA)
}



