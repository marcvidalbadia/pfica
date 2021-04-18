bcv <- function(a, phi, G, r = 2, ncomp, baseline, value = 0.1) {

  if (baseline >= 10^6) value <- 100
  if (baseline >= 10^9) value <- 1000
  if (baseline < 1) {vlambda <- ceiling(baseline)
  } else {vlambda <- baseline+value}

  fdx <- fd(a,phi)
  cov <- corpcor::cov.shrink(t(a), verbose = F)/ncol(a)
  delta <- diff(diag(nrow(a)), differences = r)
  Pr <- crossprod(delta)
  Gl <- G + (baseline*Pr)
  L <- chol(Gl)
  Li <- solve(L)
  rGram <- crossprod(Li,G)
  C2 <- rGram %*% cov %*% t(rGram); C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2)
  u <- as.matrix(svdc$u[, 1:ncomp])
  Q <- solve(Gl) %*% G
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(crossprod(Li,G)) %*% u
  gamma <- fd(b,phi)
  z <- inprod(fdx, gamma)
  KL1 <- tcrossprod(b,z)

  lambda <- vlambda
  Gl <- G + (lambda*Pr)
  L <- chol(Gl)
  Li <- solve(L)
  rGram <- crossprod(Li,G)
  C2 <- rGram %*% cov %*% t(rGram); C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2)
  u <- as.matrix(svdc$u[, 1:ncomp])
  Q <- solve(Gl) %*% G
  Qs <- expm::sqrtm(Q)
  b <-  Qs %*% solve(crossprod(Li,G)) %*% u
  gamma <- fd(b,phi)
  z <- inprod(fdx, gamma)
  KL2 <- tcrossprod(b,z)
  res.coef <- KL1-KL2
  cres <- corpcor::cov.shrink(t(res.coef),  verbose = F)/ncol(a)
  cholF <-  chol(cres)
  res.scoef <- solve(cholF) %*% res.coef
  cv <- mean(log(diag(t(res.scoef) %*% G %*% res.scoef)));
  return(cv)
}
