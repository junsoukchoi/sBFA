#' Implementation of the data-augmneted Gibbs sampler for the sparse Bayesian factor models with the Dirichlet-Laplace priors
#'
#' @param x a matrix containing data
#' @param q the number of factors to be fitted.
#' @param starting a list with each tag corresponding to a parameter name.
#' Valid tags are 'Lambda', 'U', 'Sigma', 'Phi', 'tau', 'Psi'.
#' The value portion of each tag is the parameters' starting values for MCMC.
#' @param priors a list with each tag corresponding to a parameter name.
#' Valid tags are 'Sigma' and 'Phi'.
#' The value portion of each tag defines the hyperparameters.
#' @param nmcmc the number of MCMC iterations.
#' @param nburnin the number of burn-in samples.
#' @param verbose if TRUE, progress of the sampler is printed to the screen.
#' Otherwise, nothing is printed to the screen.
#' @param nreport the interval to report the MCMC progress.
#'
#' @return a list that contains MCMC samples from the posterior distribution of the sparse Bayesian factor models with the Dirichlet-Laplace priors .
#' @export
#'
#' @examples
#' set.seed(7)
#'
#' # set the sample size n, dimension p, and the number of factors q
#' n = 100
#' p = 100
#' q = floor(log(p))
#'
#' # generate true Lambda
#' # s = log(p) nonzero elements are drawn uniformly between 1 and 2 per columns
#' Lambda = matrix(0, p, q)
#' for (k in 1 : q)
#' {
#'    id_nonzero = sample(1 : p, floor(log(p)))
#'    Lambda[id_nonzero, k] = runif(floor(log(p)), 1, 2)
#' }
#'
#' # set true Sigma to be identity
#' Sigma = diag(1, p)
#'
#' # generate data from a sparse factor model with the given Lambda and Sigma
#' Y = matrix(NA, n, p)
#' U = matrix(NA, n, q)
#' for (i in 1 : n)
#' {
#'    U[i, ] = rnorm(q)
#'    Y[i, ] = Lambda %*% U[i, ] + sqrt(Sigma) %*% rnorm(p)
#' }
#'
#' # choose the values of hyperparameters
#' priors = list()
#' priors$Sigma = c(0.1, 0.1)
#' priors$Phi   = 0.5
#'
#' # obtain starting values for MCMC from the prior distribution
#' starting = list()
#' starting$Phi    = matrix(1, p, q)
#' starting$Phi    = starting$Phi / sum(starting$Phi)
#' starting$tau    = rgamma(1, shape = p * q * priors$Phi, rate = 0.5)
#' starting$Psi    = matrix(rexp(p * q, rate = 0.5), p, q)
#' starting$Lambda = matrix(rnorm(p * q, sd = c(sqrt(starting$Psi * starting$tau^2 * starting$Phi^2))), p, q)
#' starting$U      = matrix(rnorm(n * q), n, q)
#' starting$Sigma  = diag(1, p)
#'
# run MCMC for the sparse Bayesian factor models with the Dirichlet-Laplace priors
#out = sBFA_DL(Y, q, starting, priors)
#
## calculate the posterior mean of the loading matrix
#Lambda_est = apply(out$Lambda, c(1, 2), mean)
sBFA_DL = function(x, q, starting, priors, nmcmc = 10000, nburnin = 5000, verbose = TRUE, nreport = 500)
{
   # sample size and dimension
   n = nrow(x)
   p = ncol(x)

   # initialize parameters with provided starting values
   Lambda = starting$Lambda
   U      = starting$U
   Sigma  = starting$Sigma
   Phi    = starting$Phi
   tau    = starting$tau
   Psi    = starting$Psi

   # initialize MCMC samples
   Lambda_mcmc = array(NA, dim = c(p, q, nmcmc))
   U_mcmc      = array(NA, dim = c(n, q, nmcmc))
   Sigma_mcmc  = array(NA, dim = c(p, p, nmcmc))
   Phi_mcmc    = array(NA, dim = c(p, q, nmcmc))
   tau_mcmc    = rep(NA, nmcmc)
   Psi_mcmc    = array(NA, dim = c(p, q, nmcmc))

   # initialize H's that are necessary for sampling \phi's
   H = matrix(1, p, q)

   # iterate the data-augmented Gibbs sampler for the sparse Bayesian factor models with the DL priors
   for (iter in 1 : nmcmc)
   {
      # sample \Lambda conditional on U, \Sigma, \Phi, \tau, \Psi, x
      # draw \lambda_j independently from a N(\eta_j, \Xi_j), where
      # \eta_j = \Xi_j (\sigma_j^{-2} \sum_i x_{ij} u_i),
      # \Xi_j = (\sigma_j^{-2} \sum_i u_i u_i^\top + P_j^{-1})^{-1},
      # P_j = diag(\psi_{j1} \tau^2 \phi_{j1}^1, \ldots, \psi_{jk} \tau^2 \phi_{jk}^1)
      for (j in 1 : p)
      {
         invXi = crossprod(U) / Sigma[j, j] + diag(1 / (Psi[j, ] * tau^2 * Phi[j, ]^2), q)
         invXi = (invXi + t(invXi)) / 2   # guarantee symmetry
         invXi = as(invXi, "sparseMatrix")
         chol_invXi  = Matrix::Cholesky(invXi, LDL = FALSE)
         e_Xi        = as.vector(Matrix::solve(chol_invXi, rnorm(q), system = "Lt"))
         invXi_eta   = crossprod(U, Y[ , j]) / Sigma[j, j]
         eta         = as.vector(Matrix::solve(chol_invXi, Matrix::solve(chol_invXi, invXi_eta, system = "L"), system = "Lt"))
         Lambda[j, ] = eta + e_Xi
      }

      # sample U conditional on \Lambda, \Sigma, x
      # draw U from a MN(M, I_n, V), where
      # M = (m_1^\top, \ldots, m_n^\top)^\top
      # m_i= V (\Lambda^\top \Sigma^{-1} x_i)
      # V = (\Lambda^\top \Sigma^{-1} \Lambda + I_k)^{-1}.
      invV_i    = crossprod(Lambda, Lambda / diag(Sigma)) + diag(1, q)
      invV_i    = (invV_i + t(invV_i)) / 2   # guarantee symmetry
      invV      = do.call(Matrix::bdiag, rep(list(invV_i), n))
      chol_invV = Matrix::Cholesky(invV, LDL = FALSE)
      e_V       = as.vector(Matrix::solve(chol_invV, rnorm(n * q), system = "Lt"))
      invV_m    = c(crossprod(Lambda, t(Y) / diag(Sigma)))
      nu        = as.vector(Matrix::solve(chol_invV, Matrix::solve(chol_invV, invV_m, system = "L"), system = "Lt"))
      U         = matrix(nu + e_V, n, q, byrow = TRUE)

      # sample \Sigma conditional on \Lambda, U, x
      # draw \sigma_j^{-2} independently from a Inv-Gamma(a_{\sigma^2} + n / 2, b_{\sigma^2} + \sum_{i} (y_{ij} - \lambda_j^\top u_i)^2/2)
      Sigma[diag(TRUE, p)] = 1 / rgamma(p, shape = priors$Sigma[1] + rep(n / 2, p),
                                        rate = priors$Sigma[2] + 0.5 * colSums((Y - tcrossprod(U, Lambda))^2))

      # sample \Phi conditional on \Lambda
      # draw H_{jh} independently from a giG(\alpha - 1, 1 , 2|\lambda_{jh}|)
      # set \phi_{jh} = H_{jh} / sum_{j,h} H_{jh}
      for (j in 1 : p)
      {
         for (k in 1 : q)
         {
            H[j, k] = rgig(1, lambda = priors$Phi - 1, chi = 2 * abs(Lambda[j, k]), psi = 1)
         }
      }
      Phi = H / sum(H)

      # sample \tau conditional on \Phi, \Lambda
      # draw \tau from a giG(pk(\alpha - 1), 1, 2 \sum_{j, h} |\lambda_{jh}| / \phi_{jh})
      tau = rgig(1, lambda = p * q * (priors$Phi - 1), chi = 2 * sum(abs(Lambda) / Phi), psi = 1)

      # sample \Psi conditional on \Phi, \tau, \Lambda
      # draw \psi_{jh} independently from a giG(1/2, 1, \lambda_{jl}^2 / (\tau^2 \phi_{jh}^2))
      Psi_tilde = rinvgauss(p * q, mean = c(Phi * tau / abs(Lambda)))
      Psi = matrix(1 / Psi_tilde, p, q)

      # save MCMC samples at iteration iter
      Lambda_mcmc[ , , iter] = Lambda
      U_mcmc[ , , iter]      = U
      Sigma_mcmc[ , , iter]  = Sigma
      Phi_mcmc[ , , iter]    = Phi
      tau_mcmc[iter]         = tau
      Psi_mcmc[ , , iter]    = Psi

      # print progress of our sampler
      if (verbose)
      {
         if (iter %% nreport == 0) cat("iter =", iter, "\n")
      }
   }

   # return a list of MCMC samples (afte burn-in period) for the model paramters
   out = list()
   out$Lambda = Lambda_mcmc[ , , (nburnin + 1) : nmcmc]
   out$U      = U_mcmc[ , , (nburnin + 1) : nmcmc]
   out$Sigma  = Sigma_mcmc[ , , (nburnin + 1) : nmcmc]
   out$Phi    = Phi_mcmc[ , , (nburnin + 1) : nmcmc]
   out$tau    = tau_mcmc[(nburnin + 1) : nmcmc]
   out$Psi    = Psi_mcmc[ , , (nburnin + 1) : nmcmc]
   return(out)
}
