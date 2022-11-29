pspline.kffobi <- function(fdx, ncomp = fdx$basis$nbasis, pp = 0, r = 2,
                           pr = c("fdx", "wfdx", "KL", "wKL"),
                           center = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "wKL"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project.")
  }
  
  if (center) fdx <- center.fd(fdx)

  a <- fdx$coefs
  nrep <- ncol(a)
  if (nrep < 2)
    stop("ICA not possible without replications.")

  phi <- fdx$basis
  cov <- tcrossprod(a)/nrep
  delta <- diff(diag(nrow(a)), differences = r)
  Pr <- crossprod(delta)
  G <- inprod(phi, phi)
  Gl <- G + (pp*Pr)
  L <- chol(Gl)
  Li <- solve(L)
  rGram <- crossprod(Li,G)
  C2 <- rGram %*% cov %*% t(rGram)
  C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2)
  u <- as.matrix(svdc$u[, 1:ncomp])

  Q <- solve(Gl) %*% G
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(t(Li) %*% G) %*% u
  #print(diag(t(b)%*%G%*%b)) #check norms

  b3 <-  Li %*% u
  f3 <- fd(b3,phi);
  G2 <- inprod(f3,f3)
  Gram2 <- chol(G2)
  Li3  <- solve(Gram2)
  b2 <- b3%*%Li3
  beta <- fd(b2,phi)
  #print(inprod(beta,beta)) #check orthonormality

  z <-  t(t(u)%*%rGram%*%a)
  #print(crossprod(z)/nrep) #check uncorrelatedness

  svdz <- La.svd(crossprod(z)/nrep)
  wz <- svdz$u %*% diag(c(1/sqrt(svdz$d)))%*% t(svdz$u)
  zst <- z %*% wz
  #crossprod(zst)/nrep
  nr <- sqrt(rowSums(zst^2))
  z.st <- nr * zst
  C4 <- crossprod(z.st)/(nrep * (ncomp + 2))
  svdk <- La.svd(C4)
  eigenk  <- svdk$d/sum(svdk$d)
  v <- svdk$u
  c <- b %*% v
  #diag(t(c)%*%G%*%c) #check norms
  psi <- fd(c, phi)
  #print(inprod(psi,psi)) #check orthonormality
  Ls <- chol(G)
  V2 <- Ls %*% cov %*% t(Ls)
  V2  <- (V2 + t(V2))/2
  V <- La.svd(V2)
  wa <- V$u %*% diag(c(1/sqrt(V$d))) %*% t(V$u)
  ast <- solve(Ls) %*% wa %*% Ls %*% a
  print(diag(Ls%*%crossprod(t(ast))%*%t(Ls))/nrep)
  xst <- fd(ast, phi)

  KL <- fd(b %*% t(z), phi)
  KLst <- fd(b %*% t(zst), phi)
  project <- list(fdx, KL)
  names(project) <- c("fdx", "KL")

  if (pr ==  "wKL") {
    zi <- zst %*% v
  } else if (pr ==  "wfdx") {
    zi <- t(v%*%t(b)%*%G%*%ast)
  }  else {
    zi <- inprod(project[[paste(pr)]], psi)
  }
  #print(crossprod(zi)/nrep) #check orthonormality

  FICA <- list(svdc$d,beta,z,svdk$d,psi,zi)
  names(FICA) <- c("PCA.eigv",
                   "PCA.basis",
                   "PCA.scores",
                   "ICA.eigv",
                   "ICA.basis",
                   "ICA.scores")
  return(FICA)
}
