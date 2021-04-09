kd <- function(fdx, hm = fdPar(fdx), pp = NULL, r = 2,
               pr = c("fdx", "fdx.st", "KL", "KL.st"),
               centerfd = FALSE, qmin = 2, qmax = 5) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (qmin <= 1)
    stop("init must be greater or equal than 2")
  if (qmax > fdx$basis$nbasis)
    stop("qmax is too large")
  if (length(pr) != 1 & is.character(pr)) {
    pr <- "KL.st"
  } else if (!is.character(pr)) {
    stop("Select a functional data object to project")
  }

  K <- list()
  k <- numeric()
  kurt.d <- numeric()
  nrep <- ncol(fdx$coefs)

  if (!is.null(pp)) {
    for (i in qmin:qmax) {
      sc <- pspline.kffobi(fdx, i, pp, r, pr = pr, centerfd)$scores
      for (j in 1:i) k[j] <- moments::kurtosis(sc[,j])
      K[[i]] <- k }
    } else {
    for (i in qmin:qmax) {
      sc <- kffobi(fdx, i, eigenfPar=hm, pr = pr, center = centerfd)$scores
      for (j in 1:i) k[j] <- moments::kurtosis(sc[,j])
      K[[i]] <- k }
  }

  for (i in qmin:qmax) {
    max <- which.max(unlist(K[i]))
    min <- which.min(unlist(K[i]))
    res <-  sqrt(unlist(K[i])[max]^2-unlist(K[i])[min]^2)
    kurt.d[[i]] <- res
  }

  kurt.d <- kurt.d[qmin:qmax]
  kurt.d <- as.data.frame(kurt.d)
  rownames(kurt.d) <- qmin:qmax
  return(kurt.d)
}
