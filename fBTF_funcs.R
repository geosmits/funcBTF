## To install dsp and dspCP, run the below -----
# library(devtools)
# devtools::install_github("drkowal/dsp")
# devtools::install_github("haoxuanwu/dspCP")

# Helper functions ####
preprocess_ties <- function(y_k,   # original observations
                            t_k) { # original observation times
  if (length(y_k) != length(t_k)) {
    stop('observations and observation times must have same length')
  }
  # aggregate by mean over shared observation times
  ord <- order(t_k)
  t_k <- t_k[ord]
  y_k <- y_k[ord]
  
  t_i <- unique(t_k)
  n <- length(t_i)
  
  y_i <- rep(0, n)
  m_i <- rep(0, n)
  for (i in seq_along(t_i)) {
    index <- which(t_k == t_i[i])
    m_i[i] <- length(index)
    y_i[i] <- mean(y_k[index])
  }
  
  return(list(y = y_i, t = t_i, m = m_i))
}

get_deltas <- function(t_i, # non-tied observation times
                       D) { # degree of differencing
  dt <- diff(t_i)
  if (D == 0) {
    return(NULL)
  } else if (D == 1) {
    return(dt)
  } else if (D == 2) {
    dt_prod <- dt[-1] * dt[-length(dt)]
    stopifnot(length(dt_prod) == (length(t_i) - D)) # ensure lengths are n-2
    return(dt_prod)
  } else {stop('D must be 0, 1, or 2')}
}

rescale_evol_error <- function(mu,  # beta_i
                               t_i, # non-tied observation times
                               D) { # degree of differencing
  dt <- diff(t_i)
  if (D == 0) {
    return(mu)
  } else if (D == 1) {
    return(diff(mu)/dt)
  } else if (D == 2) {
    dt_prod <- dt[-1] * dt[-length(dt)]
    return(diff(mu, differences = 2)/dt_prod)
  } else {stop('D must be 0, 1, or 2')}
}

# Sampling functions ####
t_sampleEvolParams <- function(omega, evolParams) { # takes already re-scaled evolution error as input
  nd <- length(omega)
  tau <- 0.001/sqrt(nd) # maybe can fix tau universally if we pre-scale the data? - David
  omega_norm <- omega/tau
  lambda_2 <- evolParams$lambda_2
  x_lambda_t <- extraDistr::rinvgamma(nd, 1, 1 + 1 / lambda_2)
  lambda_2 <- extraDistr::rinvgamma(nd, 1, 1 / x_lambda_t + omega_norm^2 / 2)
  sigma_wt <- pmax(sqrt(lambda_2) * tau, 1e-8)
  return(list(sigma_wt = sigma_wt, tau = tau, lambda_2 = lambda_2))
}

t_initEvolParams <- function(omega) { # takes already re-scaled evolution error as input
  nd <- length(omega)
  tau <- 0.001/sqrt(nd) # maybe can fix tau universally if we pre-scale the data? - David
  omega_norm <- omega/tau
  x_lambda_t <- rep(100, nd)
  lambda_2 <- extraDistr::rinvgamma(nd, 1, 1 / x_lambda_t + omega_norm^2 / 2)
  sigma_wt <- pmax(sqrt(lambda_2) * tau, 1e-8)
  return(list(sigma_wt = sigma_wt, tau = tau, lambda_2 = lambda_2))
}

t_initEvolParams_no <- function(y, D, omega) { # takes already re-scaled evolution error as input
  n <- length(y)
  nd <- n-D
  tau = 0.001/sqrt(nd) # maybe can fix tau universally if we pre-scale the data? - David
  omega_norm <- omega/tau
  x_lambda_t <- rep(100, nd)
  lambda_2 <- extraDistr::rinvgamma(nd, 1, 1 / x_lambda_t + omega_norm^2 / 2)
  sigma_wt <- pmax(sqrt(lambda_2) * tau, 1e-8)
  return(list(sigma_wt = sigma_wt, tau = tau, lambda_2 = lambda_2))
}

# Functional BTF function ####
fbtf <- function (y_k,        # original observations
                  t_k = NULL, # original observation times
                  D = 2,      # degree of differencing
                  nsave = 1000, nburn = 1000, nskip = 4,
                  mcmc_params = list("mu", "yhat",
                                     "evol_sigma_t2", "obs_sigma_t2",
                                     "dhs_phi", "dhs_mean"),
                  computeDIC = TRUE, verbose = TRUE) {
  
  # preprocessing and initialization
  if (D == 0) {stop('D must be 1 or 2')} # will modify btf0 later, not really needed at the moment
  
  if (is.null(t_k)) {t_k <- seq_along(y_k)} # default equally spaced times
  
  # remove ties
  tie_adj_dat <- preprocess_ties(y_k, t_k)
  y <- tie_adj_dat$y
  t <- tie_adj_dat$t
  m <- tie_adj_dat$m
  n <- length(y)
  
  deltas <- get_deltas(t, D)
  delta_sq <- if (is.null(deltas)) {NULL} else {deltas^2}
  
  loc <- t_create_loc(n - D, 1)
  
  is.missing <- which(is.na(y))
  any.missing <- (length(is.missing) > 0)
  y <- approxfun(t, y, rule = 2)(t)
  
  # shrink sigma by square root of number of ties
  sigma_e <- sd(y, na.rm = TRUE) / sqrt(m)
  sigma_et <- sigma_e
  
  chol0 <- NULL
  
  evol_sigma_t2_init <- 0.01 * sigma_et^2
  if (D > 0){
    evol_sigma_t2_init[(D+1):n] <- evol_sigma_t2_init[(D+1):n] * delta_sq
  }
  
  mu <- dsp::sampleBTF(y, obs_sigma_t2 = sigma_et^2,
                       evol_sigma_t2 = evol_sigma_t2_init,
                       D = D, chol0 = chol0)
  
  # re-scale evolution error
  omega <- rescale_evol_error(mu, t, D)
  mu0 <- as.matrix(mu[1:D, ])
  evolParams <- t_initEvolParams_no(y, D, omega)
  evolParams0 <- initEvol0(mu0)
  
  mcmc_output <- vector("list", length(mcmc_params))
  names(mcmc_output) <- mcmc_params
  if (!is.na(match("mu", mcmc_params)) || computeDIC) 
    post_mu <- array(NA, c(nsave, n))
  if (!is.na(match("yhat", mcmc_params))) 
    post_yhat <- array(NA, c(nsave, n))
  if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) 
    post_obs_sigma_t2 <- array(NA, c(nsave, n))
  if (!is.na(match("evol_sigma_t2", mcmc_params))) 
    post_evol_sigma_t2 <- array(NA, c(nsave, n))
  post_loglike <- numeric(nsave)
  
  nstot <- nburn + (nskip + 1) * (nsave)
  skipcount <- 0
  isave <- 0
  if (verbose) {timer0 = proc.time()[3]}
  
  # sampling
  for (nsi in seq_len(nstot)) {
    if (any.missing) 
      y[is.missing] <- mu[is.missing] + sigma_et[is.missing] * rnorm(length(is.missing))
    
    mu <- dsp::sampleBTF(y, obs_sigma_t2 = sigma_et^2,
                         evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                           evolParams$sigma_wt^2 * if (D > 0) delta_sq else NULL),
                         D = D, chol0 = chol0)
    
    omega <- rescale_evol_error(mu, t, D)
    mu0 <- as.matrix(mu[1:D, ])
    evolParams0 <- sampleEvol0(mu0, evolParams0, A = 1)
    evolParams <- t_sampleEvolParams(omega, evolParams)
    sigma_e <- sqrt(extraDistr::rinvgamma(n = 1, alpha = n/2, 
                                          beta = sum(m*(y-mu)^2)/2))
    sigma_et <- rep(sigma_e, n)
    
    if (nsi > nburn) {
      skipcount = skipcount + 1
      if (skipcount > nskip) {
        isave = isave + 1
        if (!is.na(match("mu", mcmc_params)) || computeDIC) 
          post_mu[isave, ] <- mu
        if (!is.na(match("yhat", mcmc_params))) 
          post_yhat[isave, ] <- mu + sigma_et * rnorm(n)
        if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) 
          post_obs_sigma_t2[isave, ] <- sigma_et^2
        if (!is.na(match("evol_sigma_t2", mcmc_params))) 
          post_evol_sigma_t2[isave, ] <- c(evolParams0$sigma_w0^2, 
                                           evolParams$sigma_wt^2 * if (D > 0) delta_sq else NULL)
        post_loglike[isave] = sum(dnorm(y, mu, sigma_et, log = TRUE))
        skipcount = 0
      }
    }
    if (verbose) {computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)}
  }
  
  # results
  if (!is.na(match("mu", mcmc_params))) 
    mcmc_output$mu <- post_mu
  if (!is.na(match("yhat", mcmc_params))) 
    mcmc_output$yhat <- post_yhat
  if (!is.na(match("obs_sigma_t2", mcmc_params))) 
    mcmc_output$obs_sigma_t2 <- post_obs_sigma_t2
  if (!is.na(match("evol_sigma_t2", mcmc_params))) 
    mcmc_output$evol_sigma_t2 <- post_evol_sigma_t2
  mcmc_output$loglike <- post_loglike
  if (computeDIC) {
    loglike_hat <- sum(dnorm(y, mean = colMeans(post_mu), 
                             sd = colMeans(sqrt(post_obs_sigma_t2)), log = TRUE))
    p_d <- c(2 * (loglike_hat - mean(post_loglike)), 2 * 
               var(post_loglike))
    DIC <- -2 * loglike_hat + 2 * p_d
    mcmc_output$DIC <- DIC
    mcmc_output$p_d <- p_d
  }
  mcmc_output$t <- t
  if (verbose) {print(paste("Total time: ", round((proc.time()[3] - timer0)), "seconds"))}
  return(mcmc_output)
}

# Unedited functions ####
# From dspCP
t_create_loc <- function (len, D) {
  if (D == 0 || D == 1) {
    row_ind = c()
    col_ind = c()
    for (i in 0:(len - 1)) {
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len - 2)) {
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i + 1)
    }
    for (i in 0:(len - 2)) {
      row_ind = c(row_ind, i + 1)
      col_ind = c(col_ind, i)
    }
  }
  if (D == 2) {
    row_ind = c()
    col_ind = c()
    for (i in 0:(len - 1)) {
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len - 2)) {
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i + 1)
    }
    for (i in 0:(len - 2)) {
      row_ind = c(row_ind, i + 1)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len - 3)) {
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i + 2)
    }
    for (i in 0:(len - 3)) {
      row_ind = c(row_ind, i + 2)
      col_ind = c(col_ind, i)
    }
  }
  list(r = row_ind, c = col_ind)
}

computeTimeRemaining <- function (nsi, timer0, nsims, nrep = 1000) 
{
  if (nsi%%nrep == 0 || nsi == 20) {
    timer = proc.time()[3]
    simsPerSec = nsi/(timer - timer0)
    secRemaining = (nsims - nsi - 1)/simsPerSec
    if (secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    }
    else {
      if (secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      }
      else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}

