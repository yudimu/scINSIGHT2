#' Estimate process for specific candidate p and seed
#'
#' @description
#' This function input the concatenated count matrix, covariate information and log library size factor,
#' and outputs a list of
#'
#' @param Y Concatenated count matrix.
#' @param X Covariate information.
#' @param logs Log library size factor using whole gene list instead of only highly variable genes.
#' @param p A specific hyper parameter p.
#' @param seed A specific initialization seed.
#' @param maxIter The number of maximum iteration. (default 5000)
#' @param alpha Step size. (default 0.01)
#' @param tol1 The stopping rule for RMSE. (default 1e-5)
#' @param tol2 The stopping rule for mean RMSE. (default 1.5e-5)
#' @param out.dir Output directory.
#'
#' @importFrom stats rnorm
#'
#' @return A list of estimates.



estimation = function(Y,
                      X,
                      logs,
                      p,
                      seed,
                      maxIter,
                      alpha,
                      tol1,
                      tol2,
                      out.dir){
  #start time
  start = Sys.time()
  # Derive dimensions
  n = dim(Y)[1]
  m = dim(Y)[2]

  # Add the intercept as a column of ones in X
  X_new = cbind(matrix(1,n,1),X)


  # X still may be NULL if no intercept and empty input X. Then d = 0
  d = 0
  if (!is.null(X_new))
    d = dim(X_new)[2]

  # Initialization
  set.seed(seed)
  sd = 1e-1
  udim = c(n,p)
  vdim = c(m,p)
  betadim = c(d,m)
  U = array(rnorm(prod(udim))/prod(dim(udim))*sd, udim)
  V = array(rnorm(prod(vdim))/prod(dim(vdim))*sd, vdim)
  beta = array(rnorm(prod(betadim))/prod(betadim)*sd, betadim)
  sigma = mean(apply(U, 1, FUN = sd))

  # Remember the last deviance
  llast = Inf

  # remember where is missing data so that we can keep replacing it
  # with more and more accurate models
  isna = is.na(Y)
  Y[isna] = mean(Y,na.rm=TRUE)

  # Since regression and latent predictions are used multiple times
  # throughout the computation, we will store them
  regpart = 0
  if (d)
    regpart = X_new %*% beta
  latentpart = multMat(U, t(V))

  #Run estimation using Rcpp
  results = iteration(U, V,  Y, beta, X_new, alpha, sigma,
                      latentpart, regpart, logs, maxIter, d,
                      n, m, llast, tol1, tol2, p)

  U = results$U
  V = results$V
  squaredU = results$squaredU
  meansquaredU = results$meansquaredU
  likelihood = results$likelihood
  beta = results$beta
  sigma = results$sigma
  latentpart = results$latentpart
  regpart = results$regpart
  logs = results$logs
  l = results$l
  i = results$i
  deviance_vec = results$dev

  #correction on U and V
  decomp = correct.uv(U, V, sigma)
  U = decomp$U
  V = decomp$V

  print(paste("Stopped after", i,"iterations"))

  #end time
  end = Sys.time()

  output = list(beta = beta,
                U = U,
                V = V,
                likelihood = likelihood,
                iter = i,
                sigma = sigma,
                time = end - start)


  #save ouput to out.dir
  saveRDS(output, file = paste0(out.dir, "Results_p=", p, "_seed=", seed, ".rds"))

  #save the three parameters for initialization in next step
  return(list(U = output$U, V = output$V, beta = output$beta, likelihood = output$likelihood[1:output$iter]))
}
