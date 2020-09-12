kd <- function(fdx, hm = fdPar(fdx), rho = NULL, r = 2, centerfd = FALSE, qmin = 2, qmax = 5) {
  if (!(inherits(fdx, "fd")))
    stop("Argument FD  not a functional data object; see fda package.")
  if (qmin <= 1)
    stop("init must be greater or equal than 2")
  if (qmax > fdx$basis$nbasis)
    stop("qmax is too large")
  K <- list()
  kurt.d <- numeric()
  if (!is.null(rho)) {
    for (i in qmin:qmax){
    K[[i]] <- pspline.kffobi(fdx, i, rho, r, pr = "fdx.st", centerfd)$kurtosis}
  } else {
    for (i in qmin:qmax){
    K[[i]] <- kffobi(fdx, i, hm, pr = "fdx.st", center = centerfd)$kurtosis
  }}
  for (i in qmin:qmax) {
    max <- which.max(unlist(K[i]))
    min <- which.min(unlist(K[i]))
    res <-  unlist(K[i])[max]-unlist(K[i])[min]
    kurt.d[[i]] <- res
  }
  kurt.d <- kurt.d[qmin:qmax]
  kurt.d <- as.data.frame(kurt.d)
  rownames(kurt.d) <- qmin:qmax
  return(kurt.d)
}
