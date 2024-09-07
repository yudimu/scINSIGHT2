#' Integration index
#'
#' @description
#' Calculate integration index for each covariate
#'
#' @param matrix The normalized latent variables.
#' @param s The number of cells that are closest neighbors of cell $i$.
#' @param l The number of cells that are in a larger neighborhood of cell $i$.
#' @param covariates A list of covariate vectors.

#' @return A list of integration index scores.
#' @importFrom RANN nn2
#'
#' @export




integration_score = function(matrix, s = 50, l = 500, covariates){

  U = matrix
  L = nrow(U)
  nns_small = nn2(data = U, k = s)
  nns_large = nn2(data = U, k = l)

  covariates_num  = length(covariates)

  scores = list()
  for (j in 1:covariates_num){
      index = covariates[[j]]
      score = c()
      for (i in 1:L){
        knn_small = nns_small$nn.idx[i,]
        knn_large = nns_large$nn.idx[i,]
        smallgroup = index[knn_small]
        largegroup = index[knn_large]

        a = 1/s*sum(smallgroup[2:s] == smallgroup[1]) - 1/l*sum(largegroup[2:l] == largegroup[1])
        score[i] = 1 - abs(a)
      }
      scores[[j]] = 1/L*sum(score)
  }
  return(scores)
}





