pspline.kffobi <- function(fdx, ncomp = fdx$basis$nbasis, pp = 0, r = 2,
                           w = c("PCA", "PCA-cor", "ZCA",
                                 "ZCA-cor", "Cholesky"),
                           pr = c("fdx", "wfdx", "KL", "wKL"),
                           center = TRUE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (ncomp <= 1) ncomp <- 2
  if (length(w) != 1 & is.character(w))
    w <- "ZCA"
  else if (!is.character(w))
    stop("Select a whitening method.")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "wKL"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project.")
  }
  nrep <- length(fdx$fdnames$reps)
  if (nrep < 2)
    stop("ICA not possible without replications.")

  if (center) fdx <- center.fd(fdx)

  a <- fdx$coefs
  phi <- fdx$basis
  cov <- tcrossprod(a)/nrep
  delta <- diff(diag(nrow(a)), differences = r)
  P <- crossprod(delta)
  G <- inprod(phi, phi)
  Gl <- G + pp*P
  L <- chol(Gl)
  Li <- solve(L)
  LiG <- crossprod(Li,G)
  C2 <- LiG %*% cov %*% t(LiG)
  C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2) #(!)
  u <- as.matrix(svdc$u[, 1:ncomp])

  Q <- solve(Gl) %*% G
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(t(Li) %*% G) %*% u
  beta <- fd(b,phi)
  z <-  t(t(u)%*%LiG%*%a)

  W <- whitening::whiteningMatrix(crossprod(z)/nrep, method = w)
  wz <- z %*% W

  nr <- sqrt(rowSums(wz^2))
  w.z <- nr * wz
  C4 <- crossprod(w.z)/(nrep * (ncomp + 2))
  svdk <- La.svd(C4)
  v <- svdk$u
  c <- b %*% v
  psi <- fd(c, phi)

  wfdx <- whiten.fd(fdx, w = w)
  KL <- fd(b %*% t(z), phi)
  wKL <- fd(b %*% t(wz), phi)
  project <- list(fdx, KL)
  names(project) <- c("fdx", "KL")
  if (pr == "wKL") {
    zi <- wz %*% v
  } else if (pr ==  "wfdx") {
    zi <- t(v%*%t(b)%*%G%*%wfdx$coefs)
  }  else {
    zi <- inprod(project[[paste(pr)]], psi)
  }

  colnames(psi$coefs) <- paste("eigenf.", c(1:ncomp), sep = " ")
  rownames(psi$coefs) <- psi$basis$names

  FICA <- list(svdc$d,beta,z,svdk$d,psi,zi,wKL)
  names(FICA) <- c("PCA.eigv",
                   "PCA.basis",
                   "PCA.scores",
                   "ICA.eigv",
                   "ICA.basis",
                   "ICA.scores",
                   "wKL")
  return(FICA)
}
