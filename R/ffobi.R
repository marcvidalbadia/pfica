ffobi  <- function (fdx, ncomp = fdx$basis$nbasis, eigenfPar = fdPar(fdx),
                    w = c("PCA", "PCA-cor","ZCA", "ZCA-cor",
                          "Varimax",  "Varimax-cor", "Cholesky"),
                    pr = c("fdx", "wfdx"), center = TRUE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (length(pr) != 1 & is.character(pr))
    pr <- "wfdx"
  else if (!is.character(pr))
    stop("Select a functional data object to project")
  if (length(w) != 1 & is.character(w))
    w <- "ZCA"
  else if (!is.character(w))
    stop("Select a whitening method")
  nrep <- length(fdx$fdnames$reps)
  if (nrep < 2)
    stop("ICA not possible without replications.")

  if (center) fdx <- center.fd(fdx)
  phi <- fdx$basis
  G <- inprod(phi, phi)

  wa <- whiten.fd(fdx, w = w)$coefs

  Lfdobj <- eigenfPar$Lfd
  pp <- eigenfPar$lambda
  rphi <- eigenfPar$fd$basis
  Gs <- eval.penalty(rphi, 0)
  if (pp > 0) {
    P <- eval.penalty(rphi, Lfdobj)
    Gs <- Gs + pp * P
  }
  Gs <- (Gs + t(Gs))/2
  J <- inprod(rphi, phi)
  Li <- solve(chol(Gs))
  LiJ <- crossprod(Li, J)

  nr <- diag(t(wa) %*% G %*% wa)
  kurt <- wa %*% diag(nr) %*% t(wa)/nrep
  C4 <- LiJ  %*% kurt %*% t(LiJ)
  C4 <- C4 + t(C4)/2

  svdk <- eigen(C4, symmetric = TRUE)
  eigenvalues <- svdk$values
  u <- as.matrix(svdk$vectors[, 1:ncomp])
  Q <- solve(Gs) %*% J
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(t(Li) %*% J) %*% u
  psi <- fd(b, rphi)

  if (pr == "fdx") {
    zi <- inprod(fdx, psi)
  } else {
    zi <- t(wa)%*%G%*%b
  }

  colnames(psi$coefs) <- paste("eigenf.", c(1:ncomp), sep = " ")
  rownames(psi$coefs) <- psi$basis$names
  FICA <- list(eigenvalues, psi, zi)
  names(FICA) <- c("ICA.eigv", "ICA.basis", "ICA.scores")
  return(FICA)
}

