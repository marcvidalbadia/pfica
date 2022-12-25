whiten.fd  <- function (fdx,
                        w = c("PCA",
                              "PCA-cor",
                              "ZCA",
                              "ZCA-cor",
                              "Varimax",
                              "Varimax-cor",
                              "Cholesky"),
                        pre.center = TRUE,
                        post.center = FALSE) {

  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object. See fda package")
  if (length(w) != 1 & is.character(w))
    w <- "ZCA"
  else if (!is.character(w))
    stop("Select a whitening method.")
  nrep <- length(fdx$fdnames$reps)
  if (nrep < 2)
    stop("ICA not possible without replications.")

  PD = function(U) return( sweep(U, 2, sign(diag(U)), "*") )

  if (pre.center) fdx <- center.fd(fdx)

  a <- fdx$coefs
  phi <- fdx$basis
  nbasis <- fdx$basis$nbasis
  G <- inprod(phi, phi)
  GChol <- chol(G)
  iGChol <- solve(GChol)

  covc <- tcrossprod(a)/nrep
  C2 <- GChol %*% covc %*% t(GChol)
  C2 <- (C2 + t(C2))/2
  dC2 <- eigen(C2, symmetric = TRUE)

  if (w == "PCA") {

    U <- PD(dC2$vectors)
    W <- diag(1/sqrt(dC2$values)) %*% t(U)

  } else if (w == "PCA-cor") {

    v <- diag(C2)
    R <- cov2cor(C2)
    eR <- eigen(R, symmetric = TRUE)
    varphi <- PD(eR$vectors)
    W <- diag(1/sqrt(eR$values)) %*% t(varphi) %*% diag(1/sqrt(v))

  } else if (w == "ZCA") {

    U <- dC2$vectors
    W <- U %*% diag(1/sqrt(dC2$values)) %*% t(U)

  } else if (w == "ZCA-cor") {

    v <- diag(C2)
    R <- cov2cor(C2)
    eR <- eigen(R, symmetric = TRUE)
    varphi <- eR$vectors
    W <- varphi %*% diag(1/sqrt(eR$values)) %*% t(varphi) %*% diag(1/sqrt(v))

  } else if (w == "Varimax") {

    rotmat <- varmx.pca.fd(pca.fd(fdx,nbasis), nharm=nbasis, nx=length(fdx$fdnames$time))$rotmat
    W <- PD(rotmat) %*% diag(1/sqrt(dC2$values)) %*% t(dC2$vectors)

  } else if (w == "Varimax-cor") {

    rotmat <- varmx.pca.fd(pca.fd(fdx,nbasis), nharm=nbasis, nx=length(fdx$fdnames$time))$rotmat
    v <- diag(C2)
    R <- cov2cor(C2)
    eR <- eigen(R, symmetric = TRUE)
    varphi <- eR$vectors
    W <- PD(rotmat) %*% diag(1/sqrt(eR$values)) %*% t(varphi) %*% diag(1/sqrt(v))

  } else if (w == "Cholesky") {

    W <- chol(solve(C2))
}

  wa <- iGChol %*% W %*% GChol %*% a

  wfdx <- fd(wa,phi)
  if (post.center) wfdx <- center.fd(wfdx);
  return(wfdx)
}
