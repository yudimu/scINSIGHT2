#' U and V correction

#' @description
#' Rotates U and V such that the covariance of U is diagonal, and V is upper-triangular with a positive diagonal
#'
#' @param U Estimated latent factor matrix.
#' @param V Estimated factor loadings.
#' @param sigma Estimated variance of latent factor.
#' @return Corrected U and V
#'
#' @importFrom whitening whiteningMatrix
#' @importFrom stats cov

correct.uv = function(U, V, sigma){
  S = cov(U)

  if (ncol(U)==1){
    return(list(u=U/sqrt(c(S)),v=V*sqrt(c(S))))
  }

  # Make cov of U identity
  W = whiteningMatrix(S)
  U = U %*% W*sigma
  V = V %*% t(solve(W))

  # Make V lower triangular
  V.qr = qr(t(V))
  U = U %*% qr.Q(V.qr)
  V = t(qr.R(V.qr))

  # Positive diagonal of V
  d = diag(V)
  V = t(sign(d)*t(V))
  U = t(sign(d)*t(U))

  list(U=U,V=V)
}

