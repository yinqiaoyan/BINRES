#' @useDynLib BINRES, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
NULL



#' Draw samples from truncated beta distribution
#' @param num Number of samples
#' @param s1 Shape1
#' @param s2 Shape2
#' @param w1 Lower bound
#' @param w2 Upper bound
#' @param by Distance between two nearby points
#' @return A vector of samples
#' @export
trunc_rbeta <- function(num, s1, s2, w1, w2, by=10^(-7)){
  x <- seq(w1, w2, by=by)
  x <- x[-1]
  lf <- (s1-1)*log(x) + (s2-1)*log(1-x)
  max_lf <- max(lf)
  ss <- sum(exp(lf-max_lf))
  prob <- exp(lf-max_lf) / ss
  elements <- sample(x, size = num, prob = prob, replace = TRUE)
  return(elements)
}



#' Convert beta random variables to weights pi_k in stick-breaking process
#' @param sticks A beta random variable vector
#' @return A vector of weights pi_k
#' @export
SticksToPi = function(sticks) {
  edges = 1.0 - cumprod(1.0 - sticks)
  return(diff(c(0, edges)))
}



#' Determine the finite set {k: P(L_i^{(2)} = k | C_i) >= v_i} for each spot i
#' @param v Augmented variable for spot i.
#' @param c Consensus clustering indicator for spot i.
#' @param K2_max Maximum in the gene-expression-specific cluster indicators L2.
#' @param G Number of genes (feature dimension of the gene expression matrix).
#' @param eta_k Mean value matrix for gene expression data. Rows: gene-expression cluster number. Columns: genes.
#' @param a_eta Mean of the normal prior for \eqn{\eta}.
#' @param b_eta Standard deviation of the normal prior for \eqn{\eta}.
#' @param psi2 Gene-expression contribution weight.
#' @param geoq Hyperparameter of the geometric distribution ranging in (0,1).
#' @return FindFiniteSet_L2 returns an R list including the following information.
#' \item{finiteSet}{vector, the finite set this function aims to find.}
#' \item{eta_k}{matrix, updated mean value matrix (maybe add some samples drawn from prior).}
#' \item{K2_max}{integer, updated value of K2_max (after possibly adding some samples).}
#' @export
FindFiniteSet_L2 = function(v, c, K2_max, G, eta_k, a_eta, b_eta, psi2, geoq){
  mass <- 0
  k <- 1
  # K_max = max(L2_Ids)   ### 这里有问题，不应该是 max(L2_Ids) 而应该是 nrow(mu2_k)！！！！！！
  finiteSet <- c()
  while(TRUE){
    if (k > K2_max) {
      eta_k <- rbind(eta_k, rnorm(G, a_eta, b_eta))
      K2_max = K2_max + 1
    }

    p_tmp = if (k == c) {psi2} else {
      (1 - psi2) * (1 - geoq)^(k - 1) * geoq / (1 - (1 - geoq)^(c - 1) * geoq)
    }

    if (p_tmp >= v) finiteSet = c(finiteSet, k)
    mass = mass + p_tmp
    rest = 1 - mass
    if (rest < v) break else k <- k + 1
  }

  return(list(finiteSet = finiteSet,
              eta_k = eta_k,
              K2_max = K2_max))
}
