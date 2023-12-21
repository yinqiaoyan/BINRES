#' Bayesian Integrative Region Segmentation in Spatially Resolved Transcriptomic Studies
#'
#' The function BINRES implements a nonparametric Bayesian integrative region
#'     segmentation method. The ST raw count matrix is normalized to account
#'     for library sizes WITHOUT taking logarithm. Dropouts and informative
#'     genes in gene expression data are taken into consideration.
#'
#' @param image_data p1*n image data matrix. There are p1 features (rows) and n spots (columns).
#' @param gene_data p2*n gene expression data matrix. Please DO NOT take logarithm in data preprocessing step..
#' @param coord Coordinates dataframe (4 columns). 1st column: first spot coordinate. 2nd column: second spot coordinate. 3rd column: first pixel coordinate. 4th column: second pixel coordinate.
#' @param platform Spatial sequencing platform. Used to determine neighborhood structure (ST = square, Visium = hex).
#' @param num_init Initial region number. Default is 5.
#' @param a_mu_elem Element of vector \eqn{a_\mu}. Default is 0.
#' @param B_mu_elem Diagonal element of matrix \eqn{B_\mu}. Default is 1.
#' @param d1 Degree of freedom for the inverse Wishart prior of \eqn{\Lambda_k}. Default is 3.
#' @param R1_elem Diagonal element of matrix R1. Default is 0.5.
#' @param a_eta Mean of the normal prior for \eqn{\eta}. Default is 1.
#' @param b_eta Standard deviation of the normal prior for \eqn{\eta}. Default is 1.
#' @param a0 Mean of the normal prior for \eqn{\lambda_0}. Default is 0.
#' @param b0 Standard deviation of the normal prior for \eqn{\lambda_0}. Default is 1.
#' @param a1 Mean of the normal distribution before truncation for \eqn{\lambda_1}. Default is -10.
#' @param b1 Standard deviation of the normal distribution before truncation for \eqn{\lambda_1}. Default is 5.
#' @param IGkappa Shape parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 2.
#' @param IGtau Scale parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 10.
#' @param dpAlpha Hyperparameter of the GEM distribution for the stick-breaking prior of \eqn{\pi_k}. That is, \eqn{\xi_i} are drawn from Be(1, dpAlpha). Default is 1.
#' @param a_beta Mean of the normal distribution before truncation for the spatial interaction parameter \eqn{\beta}. Default is 0.7.
#' @param tau_beta Standard deviation of the normal distribution before truncation for \eqn{\beta}. Default is 1.
#' @param tau0 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\xi_k^*} when k < M0. Default is 0.01.
#' @param tau1 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\beta}. Default is 0.05.
#' @param M0 A relatively large fixed positive integer. Used to determine proposal distribution form of \eqn{\xi_k^*}. Default is 50.
#' @param geoq Hyperparameter of the geometric distribution ranging in (0,1). Default is 0.5.
#' @param bernRho Hyperparameter of the Bernoulli prior for the marker gene indicator \eqn{\gamma_g}. Default is 0.1.
#' @param minPsi1 lower bound of the uniform prior for the image contribution weight \eqn{\psi_1}. Default is 0.
#' @param maxPsi1 upper bound of the uniform prior for the image contribution weight \eqn{\psi_1}. Default is 0.9.
#' @param minPsi2 lower bound of the uniform prior for the gene-expression contribution weight \eqn{\psi_2}. Default is 0.
#' @param maxPsi2 upper bound of the uniform prior for the gene-expression contribution weight \eqn{\psi_2}. Default is 0.9.
#' @param numOfMCMC Number of MCMC iterations. Default is 6000.
#' @param burnIn Number of iterations in burn-in. After burnIn the posterior samples are used and saved to estimate the unknown parameters. Default is 3000.
#' @param trunc_rbeta_by Argument "by" in the function "trunc_rbeta". Default is 10^(-5).
#' @param Is_beta_zero Logical; if TRUE, \eqn{\beta} is fixed at zero. Default is FALSE.
#' @param Is_warm_start Logical; if TRUE, warm start steps by KMeans are used to initialize L1 and L2. Default is FALSE.
#' @param Is_logNormPcaGene Logical; if TRUE, a preprocessed gene expression matrix (gene_data_pc) is used for warm start to initialize L2. Default is FALSE.
#' @param gene_data_pc n.PCs*n preprocessed gene expression matrix. Obtained by normalizing ST raw count matrix, taking logarithm, and conducting PCA. Required if Is_logNormPcaGene is TRUE. Default is NULL.
#' @param Is_print Logical; if TRUE, iteration time information of each update step are printed. Default is TRUE.
#' @param print_gap Length of iteration interval to print the current number of iterations. Default is 10.
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
#' \item{eta_0_mcmc}{matrix, the posterior samples of \eqn{\eta_0} for each gene. Rows: MCMC samples. Columns: genes.}
#' \item{sigma_g_mcmc}{matrix, the posterior samples of \eqn{\sigma_g} for each gene. Rows: MCMC samples. Columns: genes.}
#' \item{gamma_g_mcmc}{matrix, the posterior samples of \eqn{\gamma_g} for each gene. Rows: MCMC samples. Columns: genes.}
#' \item{lam1_g_mcmc}{matrix, the posterior samples of \eqn{\lambda_1} for each gene. Rows: MCMC samples. Columns: genes.}
#' \item{lam0_g_mcmc}{matrix, the posterior samples of \eqn{\lambda_0} for each gene. Rows: MCMC samples. Columns: genes.}
#' \item{pottsBeta_mcmc}{vector, the posterior samples of spatial interaction parameter \eqn{\beta}.}
#' \item{psi1_mcmc}{vector, the posterior samples of image contribution weight \eqn{\psi_1}.}
#' \item{psi2_mcmc}{vector, the posterior samples of gene-expression contribution weight \eqn{\psi_2}.}
#' \item{dpXi_mcmc}{list, each element contains the posterior sample of \eqn{\xi_k} in each MCMC iteration.}
#' \item{exeTime}{Total execution time of running the code.}
#' \item{zero_prop_spots}{vector, the spot-wise zero proportions.}
#'
#' @examples
#' library(BINRES)
#' library(aricode)
#' library(ggplot2)
#' # Import example data
#' # (1) coord: Spatial coordinates
#' # (2) image_data: Simulated image data (normally distributed)
#' # (3) gene_data: Simulated gene expression data (excess zeros and log-normally distributed)
#' # (4) true_label: True region labels of all spots
#' # (5) marker_id: ID number of marker genes
#' data(example_data_BINRES)
#' # Dimension of spatial coordinates
#' dim(coord)
#' # Dimension of image data
#' dim(image_data)
#' # Dimension of gene expression data
#' dim(gene_data)
#' # Auxiliary functions
#' getmode <- function(v) {
#'   uniqv <- unique(v)
#'   res <- uniqv[which.max(tabulate(match(v, uniqv)))]
#'   return(res)
#' }
#' # --- run BINRES ---
#' res_list = BINRES(image_data = image_data, gene_data = gene_data, coord = coord,
#'                   platform="ST", num_init=5,
#'                   minPsi1=0, maxPsi1=0.9, minPsi2=0, maxPsi2=0.9,
#'                   a_beta=0.7, tau_beta=1, tau0=0.01, tau1=0.05,
#'                   a0=0, b0=1, a1=-10, b1=5,
#'                   a_mu_elem=0, B_mu_elem=1, d1=3, R1_elem=0.5,
#'                   a_eta=1, b_eta=1, IGkappa=2, IGtau=10,
#'                   bernRho=0.1, dpAlpha=1, geoq=0.5,
#'                   Is_beta_zero=FALSE, Is_logNormPcaGene=FALSE,
#'                   numOfMCMC=1000, burnIn=500, print_gap=10,
#'                   Is_warm_start=FALSE, Is_random_seed=TRUE, random_seed=60)
#' # Posterior mode of consensus clustering C and marker-gene indicator \gamma
#' clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
#' gamma_mode = apply(res_list$gamma_g_mcmc, 2, getmode)
#' # Compared with true labels
#' table(clIds_mode, true_label)
#' cat("ARI value:", ARI(clIds_mode, true_label))
#' # Visualization
#' tmpc = clIds_mode
#' tmpc2 = tmpc
#' tmpc2[tmpc == 1] = 2
#' tmpc2[tmpc == 2] = 4
#' tmpc2[tmpc == 3] = 3
#' tmpc2[tmpc == 4] = 1
#' tmpc = tmpc2
#' plot_color=c("#CC6677", "#53B8EA", "#E2C845", "#03AF3C")
#' plot_dat <- data.frame(x = coord[,1], y = coord[,2], c = tmpc)
#' p <- ggplot(data = plot_dat, aes(x=x, y=y)) +
#'   geom_point(aes(color=factor(c)), size = 3) +
#'   theme(panel.background = element_blank(),
#'         axis.title.x=element_blank(),
#'         axis.text.x=element_blank(),
#'         axis.ticks.x=element_blank(),
#'         axis.title.y=element_blank(),
#'         axis.text.y=element_blank(),
#'         axis.ticks.y=element_blank(),
#'         legend.title = element_blank()) +
#'   scale_color_manual(values=plot_color)
#' p
#' # 95% credible interval for spatial interaction parameter \beta
#' quantile(res_list$pottsBeta_mcmc, c(0.025, 0.975))
#' # True marker gene ID
#' cat("True marker gene ID:", marker_id)
#' # Estimated marker gene ID
#' cat("Estimated marker gene ID:", which(gamma_mode == 1))
#' # Computational time
#' res_list$zero_prop_spots
#'
#' @references
#' @export
#' @importFrom stats kmeans rnorm runif dbeta dnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish rinvgamma
#' @importFrom truncnorm rtruncnorm dtruncnorm
#' @importFrom rBeta2009 rbeta
BINRES <- function(image_data, gene_data, coord, platform=c("ST", "Visium"),
                   num_init=5, a_mu_elem=0, B_mu_elem=1, d1=3, R1_elem=0.5,
                   a_eta=1, b_eta=1, a0=0, b0=1, a1=-10, b1=5, IGkappa=2, IGtau=10,
                   dpAlpha=1, a_beta=0.7, tau_beta=1, tau0=0.01, tau1=0.05,
                   M0=50, geoq=0.5, bernRho=0.1,
                   minPsi1=0, maxPsi1=0.9,
                   minPsi2=0, maxPsi2=0.9,
                   numOfMCMC=6000, burnIn=3000, trunc_rbeta_by=10^(-5),
                   Is_beta_zero=FALSE,
                   Is_warm_start=FALSE, Is_logNormPcaGene=FALSE, gene_data_pc=NULL,
                   Is_print=TRUE, print_gap=10, Is_random_seed=TRUE, random_seed=30) {

  # two types of data
  image_data = t(as.matrix(image_data))
  gene_data = t(as.matrix(gene_data))
  if (Is_logNormPcaGene) gene_data_pc = t(as.matrix(gene_data_pc))

  # location coordinates
  X_loc = as.integer(coord[, 1])
  Y_loc = as.integer(coord[, 2])

  # spot-wise zero proportions
  zero_prop_spots = unlist(Map(function(i){
    sum(gene_data[i, ] ==  0) / ncol(gene_data)
  }, 1:nrow(gene_data)))

  # data shapes
  numOfData = nrow(image_data)
  dimImage = ncol(image_data)
  G = ncol(gene_data)


  ##########################
  ###  Training process  ###
  ##########################
  cat(paste0("=== Initialization ===\n"))

  if (Is_random_seed) set.seed(random_seed)
  if (Is_warm_start) {
    # warm start using KMeans
    if (Is_logNormPcaGene) {
      kmres_L2 <- kmeans(x = gene_data_pc, centers = num_init)
    } else {
      kmres_L2 <- kmeans(x = gene_data, centers = num_init)
    }
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
  eta_0 = rnorm(G, a_eta, b_eta)

  sigma_g = sqrt(rinvgamma(G, shape = IGkappa, scale = IGtau))
  gamma_g = sample(c(1,0), G, replace = TRUE)
  lam0_g = rnorm(G, a0, b0)
  lam1_g = rtruncnorm(G, b=0, a1, b1)

  theta_latent = matrix(0, nrow=numOfData, ncol=G)
  for (i in 1:numOfData) {
    for (g in 1:G) {
      if (gene_data[i, g] != 0) {
        theta_latent[i, g] = log(gene_data[i, g])
      } else if (gamma_g[g] == 0) {
        theta_latent[i, g] = rnorm(1, mean = eta_0[g], sd = sigma_g[g])
      } else {
        theta_latent[i, g] = rnorm(1, mean = eta_k[L2_Ids[i], g], sd = sigma_g[g])
      }
    }
  }

  Z_latent = matrix(0, nrow=numOfData, ncol=G)
  for (i in 1:numOfData) {
    for (g in 1:G){
      if (gene_data[i, g] == 0) {
        Z_latent[i, g] = rtruncnorm(1, a = 0, mean = lam0_g[g] + lam1_g[g] * theta_latent[i, g], sd = 1)
      } else {
        Z_latent[i, g] = rtruncnorm(1, b = 0, mean = lam0_g[g] + lam1_g[g] * theta_latent[i, g], sd = 1)
      }

    }
  }

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
  eta_0_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  sigma_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  gamma_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  lam1_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  lam0_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
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
    act_Inds = which(gamma_g == 1)
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
        q_L[kk] = sum(dnorm(theta_latent[i, act_Inds], mean = eta_k[finiteSet_L[kk], act_Inds],
                           sd = sigma_g[act_Inds], log = TRUE))
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



    # update gamma_g (integrate out mu_{g_})
    for (g in 1:G) {
      # compute tmp1
      sumOfLogSk = Reduce(sum, Map(function(k){
        tmpInds = which(L2_Ids == k)
        Nk = length(tmpInds)
        tmpSum = sum(theta_latent[tmpInds, g])
        tmpSum2 = sum(theta_latent[tmpInds, g]^2)
        taun2 = 1 / (Nk / sigma_g[g]^2 + 1 / b_eta^2)
        mu_taun2 = tmpSum / sigma_g[g]^2 + a_eta / b_eta^2

        logSk = 0.5 * log(taun2) -
          log(b_eta) -
          (tmpSum2 / (2*sigma_g[g]^2) + (a_eta / b_eta)^2 / 2) +
          (mu_taun2^2 * taun2 / 2)
        return(logSk)
      }, unique(L2_Ids)), 0)
      tmp1 = log(bernRho) + sumOfLogSk

      # compute tmp0
      tmpSum = sum(theta_latent[, g])
      tmpSum2 = sum(theta_latent[, g]^2)
      taun2 = 1 / (numOfData / sigma_g[g]^2 + 1 / b_eta^2)
      mu_taun2 = tmpSum / sigma_g[g]^2 + a_eta / b_eta^2

      logS0 = 0.5 * log(taun2) -
        log(b_eta) -
        (tmpSum2 / (2*sigma_g[g]^2) + (a_eta / b_eta)^2 / 2) +
        (mu_taun2^2 * taun2 / 2)
      tmp0 = log(1 - bernRho) + logS0

      qq = c(tmp1, tmp0)
      qq = qq - max(qq)
      qq = exp(qq) / sum(exp(qq))
      gamma_g[g] = sample(c(1,0), 1, prob = qq)
    }


    # update eta_k
    tmpGamma_1 = which(gamma_g == 1)
    for (k in 1:K2_max){
      tmpInds = which(L2_Ids == k)
      Nk = length(tmpInds)

      if (Nk != 0) {
        tmpSum = if (Nk > 1) colSums(as.matrix(theta_latent[tmpInds, tmpGamma_1])) else theta_latent[tmpInds, tmpGamma_1]
        tmpSigma2 = 1 / (Nk / sigma_g[tmpGamma_1]^2 + 1 / b_eta^2)
        tmpMu = tmpSigma2 * (tmpSum / sigma_g[tmpGamma_1]^2 + a_eta / b_eta^2)
        eta_k[k, tmpGamma_1] = rnorm(length(tmpGamma_1), mean = tmpMu, sd = sqrt(tmpSigma2))
      } else {
        eta_k[k, tmpGamma_1] = rnorm(length(tmpGamma_1), mean = a_eta, sd = b_eta)
      }
    }



    # update eta_0
    tmpGamma_0 = which(gamma_g == 0)
    tmpSum = colSums(theta_latent[, tmpGamma_0])
    tmpSigma2 = 1 / (numOfData / sigma_g[tmpGamma_0]^2 + 1 / b_eta^2)
    tmpMu = tmpSigma2 * (tmpSum / sigma_g[tmpGamma_0]^2 + a_eta / b_eta^2)
    eta_0[tmpGamma_0] = rnorm(length(tmpGamma_0), mean = tmpMu, sd = sqrt(tmpSigma2))



    # update sigma_g
    for (g in 1:G){
      tmpMu = if (gamma_g[g] == 1) eta_k[L2_Ids, g] else eta_0[g]
      tmpSumSqu = sum((theta_latent[, g] - tmpMu)^2) / 2
      sigma_g[g] = sqrt(rinvgamma(1, shape = numOfData / 2 + IGkappa, scale = IGtau + tmpSumSqu)) # sigma_g is std!!!
    }



    # update theta_latent
    for (g in 1:G) {
      tmpTest = which(gene_data[, g] == 0)
      if (gamma_g[g] == 1) {
        ttmp = eta_k[L2_Ids[tmpTest], g]
      } else {
        ttmp = eta_0[g]
      }
      tmpSigmaSqu = 1 / (lam1_g[g]^2 + 1 / sigma_g[g]^2)
      tmpMu = ((Z_latent[tmpTest, g] - lam0_g[g]) * lam1_g[g] + ttmp / sigma_g[g]^2) * tmpSigmaSqu
      theta_latent[tmpTest, g] = rnorm(length(tmpTest), tmpMu, sqrt(tmpSigmaSqu))
    }



    # update lambda_g1 and lambda_g0
    for (g in 1:G) {
      taun0Squ = 1 / (numOfData + 1 / b0^2)

      # update lambda_g1 (integrate out lambda_g0)
      taun1Squ = 1 / (sum(theta_latent[, g]^2) - (sum(theta_latent[, g]))^2 * taun0Squ + 1 / b1^2)
      tmpInter = (sum(Z_latent[, g]) + a0 / b0^2) * sum(theta_latent[, g]) * taun0Squ
      an1 = ( sum(Z_latent[, g] * theta_latent[, g]) - tmpInter + a1 / b1^2 ) * taun1Squ
      newLam1 = rtruncnorm(1, b=0, mean = an1, sd = sqrt(taun1Squ))
      lam1_g[g] = newLam1

      # update lambda_g0
      an0 = ( sum(Z_latent[, g] - newLam1 * theta_latent[, g]) + a0 / b0^2 ) * taun0Squ
      lam0_g[g] = rnorm(1, mean = an0, sd = sqrt(taun0Squ))
    }



    # update Z_latent
    for (i in 1:numOfData) {
      tmpTest = (gene_data[i, ] == 0)
      a = b = numeric(G)
      a[!tmpTest] = -Inf
      b[tmpTest] = Inf
      Z_latent[i, ] = rtruncnorm(G, a=a, b=b, mean = lam0_g + lam1_g * theta_latent[i, ], sd = 1)
    }



    # Results --------------------------------------------------------------------
    if (mcmc > burnIn) {
      clIds_mcmc[mcmc - burnIn, ] = clIds
      L1_Ids_mcmc[mcmc - burnIn, ] = L1_Ids
      L2_Ids_mcmc[mcmc - burnIn, ] = L2_Ids
      mu_k_mcmc[[mcmc - burnIn]] <- mu_k
      Lambda_k_mcmc[[mcmc - burnIn]] <- Lambda_k
      eta_k_mcmc[[mcmc - burnIn]] <- eta_k
      eta_0_mcmc[mcmc - burnIn, ] = eta_0
      sigma_g_mcmc[mcmc - burnIn, ] = sigma_g
      gamma_g_mcmc[mcmc - burnIn, ] = gamma_g
      lam1_g_mcmc[mcmc - burnIn, ] = lam1_g
      lam0_g_mcmc[mcmc - burnIn, ] = lam0_g
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
  res_list[["eta_0_mcmc"]]      = eta_0_mcmc
  res_list[["sigma_g_mcmc"]]    = sigma_g_mcmc
  res_list[["gamma_g_mcmc"]]    = gamma_g_mcmc
  res_list[["lam1_g_mcmc"]]     = lam1_g_mcmc
  res_list[["lam0_g_mcmc"]]     = lam0_g_mcmc
  res_list[["pottsBeta_mcmc"]]  = pottsBeta_mcmc
  res_list[["psi1_mcmc"]]       = psi1_mcmc
  res_list[["psi2_mcmc"]]       = psi2_mcmc
  res_list[["dpXi_mcmc"]]       = dpXi_mcmc
  res_list[["exeTime"]]         = exeTime
  res_list[["zero_prop_spots"]] = zero_prop_spots

  return(res_list)
}






