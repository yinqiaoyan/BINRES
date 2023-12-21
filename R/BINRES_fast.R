#' Bayesian Integrative Region Segmentation in Spatially Resolved Transcriptomic Studies
#'
#' The function BINRES-fast is a scalable extension of BINRES in which PCA
#'     is used to reduce the dimensionality for the normalized gene expressions
#'     in logarithmic scale. Normal priors are assigned to top principal
#'     components.
#'
#' @param image_data p1*n image data matrix. There are p1 features (rows) and n spots (columns).
#' @param gene_data_pc n.PCs*n preprocessed gene expression matrix. Obtained by normalizing ST raw count matrix, taking logarithm, and conducting PCA.
#' @param coord Coordinates dataframe (4 columns). 1st column: first spot coordinate. 2nd column: second spot coordinate. 3rd column: first pixel coordinate. 4th column: second pixel coordinate.
#' @param platform Spatial sequencing platform. Used to determine neighborhood structure (ST = square, Visium = hex).
#' @param num_init Initial region number. Default is 5.
#' @param a_mu_elem Element of vector \eqn{a_\mu}. Default is 0.
#' @param B_mu_elem Diagonal element of matrix \eqn{B_\mu}. Default is 1.
#' @param d1 Degree of freedom for the inverse Wishart prior of \eqn{\Lambda_k}. Default is 3.
#' @param R1_elem Diagonal element of matrix R1. Default is 0.5.
#' @param a_eta Mean of the normal prior for \eqn{\eta}. Default is 0.
#' @param b_eta Standard deviation of the normal prior for \eqn{\eta}. Default is 1.5.
#' @param IGkappa Shape parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 2.
#' @param IGtau Scale parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 10.
#' @param dpAlpha Hyperparameter of the GEM distribution for the stick-breaking prior of \eqn{\pi_k}. That is, \eqn{\xi_i} are drawn from Be(1, dpAlpha). Default is 1.
#' @param a_beta Mean of the normal distribution before truncation for the spatial interaction parameter \eqn{\beta}. Default is 0.7.
#' @param tau_beta Standard deviation of the normal distribution before truncation for \eqn{\beta}. Default is 1.
#' @param tau0 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\xi_k^*} when k < M0. Default is 0.01.
#' @param tau1 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\beta}. Default is 0.05.
#' @param M0 A relatively large fixed positive integer. Used to determine proposal distribution form of \eqn{\xi_k^*}. Default is 50.
#' @param geoq Hyperparameter of the geometric distribution ranging in (0,1). Default is 0.5.
#' @param minPsi1 lower bound of the uniform prior for the image contribution weight \eqn{\psi_1}. Default is 0.
#' @param maxPsi1 upper bound of the uniform prior for the image contribution weight \eqn{\psi_1}. Default is 0.9.
#' @param minPsi2 lower bound of the uniform prior for the gene-expression contribution weight \eqn{\psi_2}. Default is 0.
#' @param maxPsi2 upper bound of the uniform prior for the gene-expression contribution weight \eqn{\psi_2}. Default is 0.9.
#' @param numOfMCMC Number of MCMC iterations. Default is 6000.
#' @param burnIn Number of iterations in burn-in. After burnIn the posterior samples are used and saved to estimate the unknown parameters. Default is 3000.
#' @param trunc_rbeta_by Argument "by" in the function "trunc_rbeta". Default is 10^(-5).
#' @param Is_beta_zero Logical; if TRUE, \eqn{\beta} is fixed at zero. Default is FALSE.
#' @param Is_warm_start Logical; if TRUE, warm start steps by KMeans are used to initialize L1 and L2. Default is FALSE.
#' @param Is_print Logical; if TRUE, iteration time information of each update step are printed. Default is TRUE.
#' @param print_gap Length of iteration interval to print the number of iterations. Default is 10.
#' @param Is_random_seed Logical; if TRUE, a random seed is used for reproducibility. Default is TRUE.
#' @param random_seed Random seed. Required if Is_random_seed is TRUE. Default is 30.
#'
#' @return BINRES returns an R list including the following information.
#' \item{clIds_mcmc}{matrix, the posterior samples of integrative region indicators for each spot. Rows: MCMC samples. Columns: n spots.}
#' \item{L1_Ids_mcmc}{matrix, the posterior samples of image-specific cluster indicators for each spot. Rows: MCMC samples. Columns: n spots.}
#' \item{L2_Ids_mcmc}{matrix, the posterior samples of gene-expression-specific cluster indicators for each spot. Rows: MCMC samples. Columns: n spots.}
#' \item{mu_k_mcmc}{list, each element contains the posterior sample of \eqn{\mu_k} for all image-specific clusters in each MCMC iteration.}
#' \item{Lambda_k_mcmc}{list, each element contains the posterior sample of \eqn{\Lambda_k} for all image-specific clusters in each MCMC iteration.}
#' \item{eta_k_mcmc}{list, each element contains the posterior sample of \eqn{\eta_k} for all gene-expression-specific clusters in each MCMC iteration.}
#' \item{sigma_g_mcmc}{matrix, the posterior samples of \eqn{\sigma_g} for each gene. Rows: MCMC samples. Columns: PCs.}
#' \item{pottsBeta_mcmc}{vector, the posterior samples of spatial interaction parameter \eqn{\beta}.}
#' \item{psi1_mcmc}{vector, the posterior samples of image contribution weight \eqn{\psi_1}.}
#' \item{psi2_mcmc}{vector, the posterior samples of gene-expression contribution weight \eqn{\psi_2}.}
#' \item{dpXi_mcmc}{list, each element contains the posterior sample of \eqn{\xi_k} in each MCMC iteration.}
#' \item{exeTime}{Total execution time of running the code.}
#'
#' @references
#' @export
#' @importFrom stats kmeans rnorm runif dbeta dnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish rinvgamma
#' @importFrom truncnorm rtruncnorm dtruncnorm
#' @importFrom rBeta2009 rbeta
BINRES_fast <- function(image_data, gene_data_pc, coord, platform=c("ST", "Visium"),
                        num_init=5, a_mu_elem=0, B_mu_elem=1, d1=3, R1_elem=0.5,
                        a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10,
                        dpAlpha=1, a_beta=0.7, tau_beta=1, tau0=0.01, tau1=0.05,
                        M0=50, geoq=0.5,
                        minPsi1=0, maxPsi1=0.9,
                        minPsi2=0, maxPsi2=0.9,
                        numOfMCMC=6000, burnIn=3000, trunc_rbeta_by=10^(-5),
                        Is_beta_zero=FALSE, Is_warm_start=FALSE,
                        Is_print=TRUE, print_gap=10,
                        Is_random_seed=TRUE, random_seed=30) {

  # two types of data
  image_data = t(as.matrix(image_data))
  gene_data_pc = t(as.matrix(gene_data_pc))

  # location coordinates
  X_loc = as.integer(coord[, 1])
  Y_loc = as.integer(coord[, 2])

  # data shapes
  numOfData = nrow(image_data)
  dimImage = ncol(image_data)
  G = ncol(gene_data_pc)


  ##########################
  ###  Training process  ###
  ##########################
  cat(paste0("=== Initialization ===\n"))

  if (Is_random_seed) set.seed(random_seed)
  if (Is_warm_start) {
    # warm start using KMeans
    kmres_L2 <- kmeans(x = gene_data_pc, centers = num_init)
    kmres_L1 <- kmeans(x = image_data, centers = num_init)

    # label switching
    init_label1 = kmres_L1[["cluster"]]
    uniq = unique(init_label1)
    for (i in 1:length(uniq)){
      init_label1[kmres_L1[["cluster"]] == uniq[i]] = i
    }

    init_label2 = kmres_L2[["cluster"]]
    uniq = unique(init_label2)
    for (i in 1:length(uniq)){
      init_label2[kmres_L2[["cluster"]] == uniq[i]] = i
    }

    L2_Ids = init_label2
    L1_Ids = init_label1
    clIds  = L2_Ids

    K1_max = max(L1_Ids)
    K2_max = max(L2_Ids)
  } else {
    K1_max = num_init
    K2_max = num_init

    clIds = sample(1:num_init, numOfData, replace = T)
    L1_Ids = clIds
    L2_Ids = clIds
  }


  # data source 1
  mu_k = vector(mode = "list", length = K1_max)
  Lambda_k = vector(mode = "list", length = K1_max)
  a_mu = matrix(a_mu_elem, nrow = dimImage, ncol = 1)
  B_mu = diag(B_mu_elem, dimImage, dimImage)
  # d1 = dimImage
  R1 = diag(R1_elem, dimImage, dimImage)
  for (k in 1:K1_max){
    mu_k[[k]] = rmvnorm(1, a_mu, B_mu)
    Lambda_k[[k]] = riwish(d1, R1)
  }

  # data source 2
  eta_k = matrix(rnorm(K2_max * G, a_eta, b_eta), nrow = K2_max)
  sigma_g = sqrt(rinvgamma(G, shape = IGkappa, scale = IGtau))

  if (!Is_beta_zero) {
    pottsBeta = rtruncnorm(1, a = 0, mean = a_beta, sd = tau_beta)
  } else {
    pottsBeta = 0
  }



  psi1 = runif(1, minPsi1, maxPsi1)
  psi2 = runif(1, minPsi2, maxPsi2)

  # initialize dpXi
  dpXi <- rbeta(M0, 1, dpAlpha)



  ### store MCMC samples
  clIds_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = numOfData)
  L1_Ids_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = numOfData)
  L2_Ids_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = numOfData)
  mu_k_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)
  Lambda_k_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)
  eta_k_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)
  sigma_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  pottsBeta_mcmc = array(0, dim = numOfMCMC - burnIn)
  psi1_mcmc = array(0, dim = numOfMCMC - burnIn)
  psi2_mcmc = array(0, dim = numOfMCMC - burnIn)
  dpXi_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)



  cat(paste0("=== MCMC Iterations ===\n"))

  ssTime = Sys.time()
  for (mcmc in 0:numOfMCMC) {

    if (!Is_beta_zero) {
      # update dpXi and pottsBeta (double MH)
      # [-] proposal of \xi_k: tnorm (k <= M0); Be(1, alpha) (k > M0)
      # [-] proposal of \beta: tnorm

      tmplen = length(dpXi)
      dpXi_proposal = numeric(tmplen)
      dpXi_proposal[1:M0] = rtruncnorm(M0, a = 0, b = 1, mean = dpXi[1:M0], sd = tau0)
      beta_proposal = rtruncnorm(1, a = 0, mean = pottsBeta, sd = tau1)
      if (tmplen > M0) dpXi_proposal[(M0 + 1):tmplen] = rbeta(tmplen - M0, 1, dpAlpha)

      if ( !(beta_proposal < 0 | sum(dpXi_proposal < 0 | dpXi_proposal > 1) > 0) ) {
        # Step 2: Propose new C* from P(C | beta*, {pi_k*})
        aux_c = clIds
        for (i in 1:numOfData){
          tmpNei = FindNeighbors_Rcpp(i - 1, X_loc, Y_loc, platform) + 1
          uniq_tmpNei_c = unique(aux_c[tmpNei])
          aux_pi = SticksToPi(dpXi_proposal)

          # compute normalizing constant NC_i
          NCi = 1
          for (m in uniq_tmpNei_c) {
            NCi = NCi + aux_pi[m] * ( exp(beta_proposal * sum(aux_c[tmpNei] == m)) - 1 )
          }

          aux_p = aux_pi / NCi
          for (m in uniq_tmpNei_c) {
            aux_p[m] = aux_pi[m] * exp(beta_proposal * sum(aux_c[tmpNei] == m)) / NCi
          }

          # inverse-cdf
          ui = runif(1)
          sumP = sum(aux_p)

          while (ui > sumP) {
            newTau = rbeta(2, 1, dpAlpha)  # two new sampled tau for dpXi and dpXi_proposal
            dpXi = c(dpXi, newTau[1])
            newPi = newTau[2] * prod(1.0 - dpXi_proposal)
            dpXi_proposal = c(dpXi_proposal, newTau[2])
            aux_p = c(aux_p, newPi / NCi)
            sumP = sumP + newPi / NCi
          }

          csum_aux_p = cumsum(aux_p)
          aux_c[i] = sum(ui > csum_aux_p) + 1

        }

        pi_old = SticksToPi(dpXi)
        pi_new = SticksToPi(dpXi_proposal)


        # Step 3: accept proposal with probability min(1,r)
        # log ratio
        # prior
        logfrac1 = sum(dbeta(dpXi_proposal[1:M0], 1, dpAlpha, log = TRUE)) +
          log(dtruncnorm(beta_proposal, a=0, mean = a_beta, sd = tau_beta))
        logfrac2 = sum(dbeta(dpXi[1:M0], 1, dpAlpha, log = TRUE)) +
          log(dtruncnorm(pottsBeta, a=0, mean = a_beta, sd = tau_beta))

        # C
        logfrac1 = logfrac1 +
          sum(log(pi_new[clIds])) + sum(log(pi_old[aux_c])) +
          ComputePottsDist_Rcpp(beta_proposal, clIds, X_loc, Y_loc, platform) + ComputePottsDist_Rcpp(pottsBeta, aux_c, X_loc, Y_loc, platform)
        logfrac2 = logfrac2 +
          sum(log(pi_old[clIds])) + sum(log(pi_new[aux_c])) +
          ComputePottsDist_Rcpp(pottsBeta, clIds, X_loc, Y_loc, platform) + ComputePottsDist_Rcpp(beta_proposal, aux_c, X_loc, Y_loc, platform)

        # proposal
        logfrac1 = logfrac1 +
          sum(log(dtruncnorm(dpXi[1:M0], a=0, b=1, mean = dpXi_proposal[1:M0], sd = tau0))) +
          log(dtruncnorm(pottsBeta, a=0, mean = beta_proposal, sd = tau1))
        logfrac2 = logfrac2 +
          sum(log(dtruncnorm(dpXi_proposal[1:M0], a=0, b=1, mean = dpXi[1:M0], sd = tau0))) +
          log(dtruncnorm(beta_proposal, a=0, mean = pottsBeta, sd = tau1))

        ratio = logfrac1 - logfrac2
        ratio = exp(ratio)
        prob = min(1, ratio)

        tmpU = runif(1)
        if (tmpU < prob) {
          pottsBeta = beta_proposal
          dpXi = dpXi_proposal
          # print("accept")
        } else {
          # print("reject")
        }

      }
    } else {
      # pottsBeta is fixed at zero
      # update dpXi (same procedure in Dirichlet process)

      dpXi = dpXi[1:max(clIds)]
      for (k in 1:max(clIds)) {
        Nk_dpXi = sum(clIds == k)
        N_larger_k_dpXi = sum(clIds > k)
        dpXi[k] = rbeta(1, 1 + Nk_dpXi, dpAlpha + N_larger_k_dpXi)
      }
    }





    # update psi1 and psi2
    num_eq1 = sum(L1_Ids == clIds)
    num_eq2 = sum(L2_Ids == clIds)
    psi1 = trunc_rbeta(1, 1 + num_eq1, 1 + numOfData - num_eq1, minPsi1, maxPsi1, trunc_rbeta_by)
    psi2 = trunc_rbeta(1, 1 + num_eq2, 1 + numOfData - num_eq2, minPsi2, maxPsi2, trunc_rbeta_by)



    # update c_lw
    Pi = SticksToPi(dpXi)
    resC = UpdateClIds(dpAlpha, pottsBeta, psi1, psi2, geoq, Pi, dpXi, clIds,
                       L1_Ids, L2_Ids, X_loc, Y_loc, platform)

    clIds = resC$clIds
    dpXi = resC$dpXi



    # update L1_lw
    resL1 = UpdateL1(K1_max, image_data, mu_k, Lambda_k, a_mu, B_mu, d1, R1,
                     psi1, geoq, clIds, L1_Ids)
    L1_Ids = resL1$L1_Ids
    mu_k = resL1$mu_k
    Lambda_k = resL1$Lambda_k

    K1_max = max(L1_Ids)
    mu_k <- mu_k[1:K1_max]
    Lambda_k <- Lambda_k[1:K1_max]



    # update L2_lw
    for (i in 1:numOfData){
      # update v2_lw
      p_tmp = if (L2_Ids[i] == clIds[i]) {psi2} else {
        (1 - psi2) * (1 - geoq)^(L2_Ids[i] - 1) * geoq / (1 - (1 - geoq)^(clIds[i] - 1) * geoq)
      }
      v_i = runif(1, min=0, max=p_tmp)

      res         = FindFiniteSet_L2(v_i, clIds[i], K2_max, G, eta_k, a_eta, b_eta, psi2, geoq)
      finiteSet_L = res$finiteSet
      eta_k       = res$eta_k
      K2_max      = res$K2_max

      q_L = numeric(length(finiteSet_L))
      for (kk in 1:length(finiteSet_L)) {
        q_L[kk] = sum(dnorm(gene_data_pc[i, ], mean = eta_k[finiteSet_L[kk], ],
                            sd = sigma_g, log = TRUE))
      }
      q_L = q_L - max(q_L)
      q_L = exp(q_L) / sum(exp(q_L))
      L2_Ids[i] = if (length(finiteSet_L) > 1) sample(finiteSet_L, 1, prob=q_L) else finiteSet_L
    }
    K2_max = max(L2_Ids)
    eta_k <- eta_k[1:K2_max, ]



    # update mu_k and Lambda_k
    for (k in 1:K1_max){
      tmpInds = which(L1_Ids == k)
      Nk = length(tmpInds)

      if (Nk != 0) {
        tmpSum = if (Nk > 1) as.matrix(colSums(image_data[tmpInds, ])) else as.matrix(image_data[tmpInds, ])
        tmpSigma = solve(solve(B_mu) + Nk * solve(Lambda_k[[k]]))
        tmpMu = tmpSigma %*% (solve(B_mu) %*% a_mu + solve(Lambda_k[[k]]) %*% tmpSum)
        mu_k[[k]] = rmvnorm(1, mean = tmpMu, sigma = tmpSigma)

        if (Nk > 1) {
          tmpSumSumT = Reduce('+', Map(function(i){
            t(image_data[i, ] - mu_k[[k]]) %*% (image_data[i, ] - mu_k[[k]])
          }, tmpInds))  # use '+' not sum
        } else {
          tmpSumSumT = t(image_data[tmpInds, ] - mu_k[[k]]) %*% (image_data[tmpInds, ] - mu_k[[k]])
        }

        Sigman = R1 + tmpSumSumT
        Lambda_k[[k]] = riwish(d1 + Nk, Sigman)
      } else {
        mu_k[[k]] = rmvnorm(1, mean = a_mu, sigma = B_mu)
        Lambda_k[[k]] = riwish(d1, R1)
      }

    }



    # update eta_k
    for (k in 1:K2_max){
      tmpInds = which(L2_Ids == k)
      Nk = length(tmpInds)

      if (Nk != 0) {
        tmpSum = if (Nk > 1) colSums(gene_data_pc[tmpInds, ]) else gene_data_pc[tmpInds, ]
        tmpSigma2 = 1 / (Nk / sigma_g^2 + 1 / b_eta^2)
        tmpMu = tmpSigma2 * (tmpSum / sigma_g^2 + a_eta / b_eta^2)
        eta_k[k, ] = rnorm(G, mean = tmpMu, sd = sqrt(tmpSigma2))
      } else {
        eta_k[k, ] = rnorm(G, mean = a_eta, sd = b_eta)
      }
    }



    # update sigma_g
    for (g in 1:G){
      tmpSumSqu = sum((gene_data_pc[, g] - eta_k[L2_Ids, g])^2) / 2
      sigma_g[g] = sqrt(rinvgamma(1, shape = numOfData / 2 + IGkappa, scale = IGtau + tmpSumSqu)) # sigma_g is std
    }



    # Results --------------------------------------------------------------------
    if (mcmc > burnIn) {
      clIds_mcmc[mcmc - burnIn, ] = clIds
      L1_Ids_mcmc[mcmc - burnIn, ] = L1_Ids
      L2_Ids_mcmc[mcmc - burnIn, ] = L2_Ids
      mu_k_mcmc[[mcmc - burnIn]] <- mu_k
      Lambda_k_mcmc[[mcmc - burnIn]] <- Lambda_k
      eta_k_mcmc[[mcmc - burnIn]] <- eta_k
      sigma_g_mcmc[mcmc - burnIn, ] = sigma_g
      pottsBeta_mcmc[mcmc - burnIn] = pottsBeta
      psi1_mcmc[mcmc - burnIn] = psi1
      psi2_mcmc[mcmc - burnIn] = psi2
      dpXi_mcmc[[mcmc - burnIn]] = dpXi
    }

    ## Output
    if (Is_print) {
      if (mcmc <= burnIn) {
        if (mcmc==0) cat(" Burn-in:") else if (mcmc/print_gap == floor(mcmc/print_gap)) cat(paste0("+++", mcmc))
      } else {
        if (mcmc==burnIn+1) cat("\n MCMC sampling:") else if (mcmc > burnIn & mcmc/print_gap == floor(mcmc/print_gap)) cat(paste0("...", mcmc))
      }
    }

  }
  eeTime = Sys.time()
  exeTime = eeTime - ssTime

  cat(paste0("\n=== End Train ===\n"))

  # return results list
  res_list = vector("list")
  res_list[["clIds_mcmc"]]      = clIds_mcmc
  res_list[["L1_Ids_mcmc"]]     = L1_Ids_mcmc
  res_list[["L2_Ids_mcmc"]]     = L2_Ids_mcmc
  res_list[["mu_k_mcmc"]]       = mu_k_mcmc
  res_list[["Lambda_k_mcmc"]]   = Lambda_k_mcmc
  res_list[["eta_k_mcmc"]]      = eta_k_mcmc
  res_list[["sigma_g_mcmc"]]    = sigma_g_mcmc
  res_list[["pottsBeta_mcmc"]]  = pottsBeta_mcmc
  res_list[["psi1_mcmc"]]       = psi1_mcmc
  res_list[["psi2_mcmc"]]       = psi2_mcmc
  res_list[["dpXi_mcmc"]]       = dpXi_mcmc
  res_list[["exeTime"]]         = exeTime

  return(res_list)
}






