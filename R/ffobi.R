ffobi <- function(fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
                  pr = c("fdx", "fdx.st"),
                  shrinkage = FALSE, center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (length(pr) != 1 & is.character(pr))
    pr <- "fdx.st"
  else if (!is.character(pr)) 
    stop("Select a functional data object to project")
  
  if (center) fdx <- center.fd(fdx)
  a <- fdx$coefs
  nrep <- ncol(a)
  
  if (nrep < 2)
    stop("ICA not possible without replications.")
  else if (!is.character(pr))
    stop("Select a functional data object to project")

  phi <- fdx$basis
  G <- inprod(phi, phi)
  L1 <- chol(G)
  W1 <- solve(L1)

  if (shrinkage == TRUE) covc <- corpcor::cov.shrink(t(a), verbose = F)/nrep
  else covc <- tcrossprod(a)/nrep

  C2 <-  L1 %*% tcrossprod(covc, L1)
  C2 <- (C2 + t(C2))/2
  dC2 <- La.svd(C2)
  wa <- dC2$u %*% tcrossprod(diag((1/dC2$d)^0.5), dC2$u)
  asta <- W1 %*% wa %*% L1 %*% a

  Lfdobj <- eigenfPar$Lfd
  lambda <- eigenfPar$lambda
  rphi <- eigenfPar$fd$basis
  R <- eval.penalty(rphi, Lfdobj)
  L <- G + lambda * R
  L <- (L + t(L))/2
  J <- inprod(rphi, phi)
  W <- solve(chol(L))
  rGram <- crossprod(W, J)

  nr <- c()
  for (i in 1:ncol(asta)) nr[i] <- (t(asta[,i]) %*% J %*% asta[,i]);
  ast <- asta %*% diag(nr)
  kurt <- tcrossprod(ast)/nrep
  C4 <- rGram %*% kurt %*% t(rGram)
  C4  <- (C4 + t(C4))/2

  #C4 <- rGram%*%asta%*%t(asta)%*%t(rGram)%*%rGram%*%asta%*%t(asta)%*%t(rGram)
  #C4  <- (C4 + t(C4))/2

  svdk <- La.svd(C4)
  u <- as.matrix(svdk$u[, 1:ncomp])
  eigenk <- svdk$d[1:ncomp]/sum(svdk$d)
  b <- W %*% u
  psi <- fd(b, rphi)
  xst <- fd(asta, phi)
  project <- list(fdx, xst)
  names(project) <- c("fdx", "fdx.st")
  zi <- inprod(project[[paste(pr)]], psi)

  if (plotfd) {
    oldpar <- par(mfrow=c(1,2), no.readonly = TRUE)
    on.exit(par(oldpar))
    for (j in 1:ncomp){
      plot(psi[j])
      title(paste('IC', j, "Kurt:", round(eigenk[j],3)))
      par(ask=T)}; par(ask=F)}
  colnames(psi$coefs) <- paste("eigenf.", c(1:ncomp), sep = "-")
  rownames(psi$coefs) <- psi$basis$names
  FICA <- list(psi, eigenk, zi)
  names(FICA) <- c("eigenbasis", "kurtosis", "scores")
  return(FICA)
}
