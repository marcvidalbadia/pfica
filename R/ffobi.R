ffobi <- function(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
                  pr = c("fdx", "fdx.st"),
                  shrinkage = FALSE, center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (length(pr) != 1 & is.character(pr))
    pr <- "fdx.st"
  else if (!is.character(pr)) 
    stop("Select a functional data object to project")
  if (center) 
    fdx <- center.fd(fdx);
  
  a <- fdx$coefs 
  nrep <- dim(a)[2]
  phi <- fdx$basis
  J <- inprod(phi, phi)
  rGram1 <- chol(J)
  W1 <- solve(rGram1)
  
  if (shrinkage == TRUE) 
    covc <- corpcor::cov.shrink(t(a), verbose = F)
  else covc <- tcrossprod(a);
  
  C2 <-  rGram1 %*% covc %*% t(rGram1)
  C2 <- (C2 + t(C2))/2
  dC2 <- La.svd(C2) 
  wa <- dC2$u %*% tcrossprod(diag((1/dC2$d)^0.5), dC2$u)
  asta <- W1 %*% wa %*% rGram1 %*% a
  
  Lfdobj <- eigenfPar$Lfd
  lambda <- eigenfPar$lambda
  rphi <- eigenfPar$fd$basis
  L <- eval.penalty(rphi, 0)
  if (lambda > 0) {
    R <- eval.penalty(rphi, Lfdobj)
    L <- L + lambda * R
  }; L <- (L + t(L))/2
  J <- inprod(rphi, phi)
  W <- solve(chol(L))
  rGram <- crossprod(W, J)
  
  nr <- numeric()
  for (i in 1:ncol(asta)) nr[i] <- (t(asta[,i]) %*% J %*% asta[,i]);
  ast <- asta %*% diag(nr)
  kurt <- tcrossprod(ast)/nrep
  C4 <- rGram %*% kurt %*% t(rGram)
  C4  <- (C4 + t(C4))/2
  
  #Alternative discretized version of the Kurtosis operator
  #C4 <- rGram%*%asta%*%t(asta)%*%t(rGram)%*%rGram%*%asta%*%t(asta)%*%t(rGram)
  #C4  <- (C4 + t(C4))/2
  
  svdk <- La.svd(C4)
  u <- as.matrix(svdk$u[, 1:ncomp])
  eigv <- svdk$d[1:ncomp]
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
  FICA <- list(h, eigv, kurt, zi)
  names(FICA) <- c("eigenbasis", "eigenvalues", "kurtosis", "scores")
  return(FICA)
}
