pspline.kffobi <- function(fdx, ncomp = fdx$basis$nbasis, rho = 0, r = 2,
                           pr = c("fdx", "fdx.st", "KL", "KL.st"),
                           center = FALSE, plotfd = FALSE) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "KL.st"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project")
  }
  if (center) {
    fdx <- center.fd(fdx)
  }
  a <- fdx$coefs
  phi <- fdx$basis
  cov <- tcrossprod(a)/ncol(a)
  delta <- diff(diag(nrow(a)), differences = r)
  Pr <- t(delta) %*% delta
  J <- inprod(phi, phi)
  L <- J + (rho*Pr)
  Li <- solve(chol(L))
  rGram <- crossprod(Li,J)
  C2 <- rGram %*% cov %*% t(rGram); C2  <- (C2 + t(C2))/2
  svdc <- La.svd(C2)
  u <- as.matrix(svdc$u[, 1:ncomp])
  b <- Li %*% u
  f <- fd(b, phi)
  z <- inprod(fdx, f)
  svdz <- La.svd(crossprod(z))
  ds <- diag(c(1/sqrt(svdz$d)))
  wz <- ds %*% t(svdz$u)
  zst <- z %*% wz
  nr <- sqrt(rowSums(zst^2))
  zst <- nr * zst
  C4 <- crossprod(zst)/(ncol(a) * (ncomp + 2))
  svdk <- La.svd(C4)
  eigenk  <- svdk$d
  v <- svdk$u
  c <- b %*% v
  h <- fd(c, phi)
  V <- La.svd(cov)
  wa <- V$u %*% diag(c(1/sqrt(V$d))) %*% t(V$u)
  wa <- sqrt(ncol(a)-1)*wa
  ast <- wa %*% a
  xst <- fd(ast, phi)
  KL <- fd(b %*% t(z), phi)
  KLst <- fd(b %*% t(zst), phi)
  project <- list(fdx, xst, KL, KLst)
  names(project) <- c("fdx", "fdx.st", "KL", "KL.st")
  zi <- inprod(project[[paste(pr)]], h)
  kurt <- vector()
  for (i in 1:ncomp) {
    K <- kurtosis(zi[,i])
    kurt[i]=K}
  if (plotfd) {
    oldpar <- par(mfrow=c(2,3), no.readonly = TRUE)
    on.exit(par(oldpar))
    for (j in 1:ncomp){
      txt <- paste('EF', j, "IC Kurt:", round(kurt[j],3))
      plot(h[j], col="#6495ED", bty="n", main=txt)
      par(ask=T)}; par(ask=F)}
  colnames(h$coefs) <- paste("eigenf.", c(1:ncomp), sep = "-")
  rownames(h$coefs) <- h$basis$names
  FICA <- list(h, kurt, zi)
  names(FICA) <- c("eigenbasis", "kurtosis", "scores")
  return(FICA)
}
