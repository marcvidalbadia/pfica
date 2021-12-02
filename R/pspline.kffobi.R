pspline.kffobi <- function(fdx, ncomp = fdx$basis$nbasis, pp = 0, r = 2,
                           pr = c("fdx", "fdx.st", "KL", "KL.st"),
                           shrinkage = FALSE, center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "KL.st"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project.")
  }
  if (center) {
    fdx <- center.fd(fdx)
  }
  a <- fdx$coefs
  nrep <- ncol(a)
  if (nrep < 2)
    stop("ICA not possible without replications.")
  phi <- fdx$basis
  if (shrinkage == TRUE) {
    cov <- corpcor::cov.shrink(t(a), verbose = F)/nrep
  } else { cov <- tcrossprod(a)/nrep }
  delta <- diff(diag(nrow(a)), differences = r)
  Pr <- crossprod(delta)
  G <- inprod(phi, phi)
  Gl <- G + (pp*Pr)
  L <- chol(Gl)
  Li <- solve(L)
  rGram <- crossprod(Li,G)
  C2 <- rGram %*% cov %*% t(rGram); C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2)
  u <- as.matrix(svdc$u[, 1:ncomp])
  Q <- solve(Gl) %*% G
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(t(Li) %*% G) %*% u
  #diag(t(b)%*%G%*%b) #check norms
  beta <- fd(b,phi)
  z <- inprod(fdx, beta)
  svdz <- La.svd(crossprod(z)/nrep)
  wz <- svdz$u %*% diag(c(1/sqrt(svdz$d)))%*% t(svdz$u)
  zst <- z %*% wz
  nr <- sqrt(rowSums(zst^2))
  z.st <- nr * zst
  C4 <- crossprod(z.st)/(nrep * (ncomp + 2))
  svdk <- La.svd(C4)
  eigenk  <- svdk$d/sum(svdk$d)
  v <- svdk$u
  c <- b %*% v
  #diag(t(c)%*%G%*%c) #check norms
  psi <- fd(c, phi)
  Ls <- chol(G)
  V2 <- Ls %*% cov %*% t(Ls)
  V2  <- (V2 + t(V2))/2
  V <- La.svd(V2)
  wa <- V$u %*% diag(c(1/sqrt(V$d))) %*% t(V$u)
  ast <- solve(Ls) %*% wa %*% Ls %*% a
  xst <- fd(ast, phi)
  KL <- fd(b %*% t(z), phi)
  KLst <- fd(b %*% t(zst), phi)
  project <- list(fdx, xst, KL, KLst)
  names(project) <- c("fdx", "fdx.st", "KL", "KL.st")
  zi <- inprod(project[[paste(pr)]], psi)

  kurt <- vector()
  for (i in 1:ncomp) kurt[i] <- moments::kurtosis(zi[,i])

  if (plotfd) {
    oldpar <- par(mfrow=c(1,2), no.readonly = TRUE)
    on.exit(par(oldpar))
    for (j in 1:ncomp){
      plot(psi[j])
      title(paste('IC', j, "Kurt:", round(eigenk[j],3)))
      par(ask=T)}; par(ask=F)}
  colnames(psi$coefs) <- paste("eigenf.", c(1:ncomp), sep = "-")
  rownames(psi$coefs) <- psi$basis$names

  FICA <- list(svdc$d,beta,z,svdk$d,psi,zi,kurt)
  names(FICA) <- c("PCA.eigv",
                   "PCA.basis",
                   "PCA.scores",
                   "ICA.eigv",
                   "ICA.basis",
                   "ICA.scores",
                   "ICA.kurtosis")
  return(FICA)
}
