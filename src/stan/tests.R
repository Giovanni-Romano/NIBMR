# Run with Rscript --vanilla src/stan/tests.R from the R directory.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
source(file.path(here, "prepare-data.R"))
stopifnot(requireNamespace("rstan", quietly = TRUE))

# Load only function definitions, avoiding library calls or execution of samplers.
original <- new.env(parent = globalenv())
for (file in c("01a-utils.R", "01b-mcmc.R")) {
  for (expr in parse(file.path(here, "..", file))) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=") &&
        length(expr) == 3L && is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("function"))) eval(expr, original)
  }
}
close <- function(x, y, tolerance = 1e-7) {
  if (!isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))
    stop("Numerical mismatch: ", paste(x, collapse = ","), " versus ", paste(y, collapse = ","))
}
rejects <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))

set.seed(410)
n <- 18L
x <- matrix(rnorm(n), ncol = 1)
z <- matrix(rnorm(n * 2), ncol = 2)
spatial <- matrix(rnorm(n * 2), ncol = 2)
y <- rnorm(n)
freq <- rep(c(1, 2, 0), length.out = n)
group <- rep(c("b", "a", "c"), length.out = n)
base_args <- list(y = y, k = c(1.3, 2.1), c = 0.03,
                  X_cont = x, Z_cont = list(z), d = 2L,
                  a_tau = c(2, 3), b_tau = c(1, 2),
                  intercept_mean = 0.2, intercept_sd = 2, freq = freq)
prepared <- do.call(prepare_mcmc_iwls_stan_data, c(base_args, list(w = 2)))
stopifnot(identical(as.integer(prepared$data$tau_idx), c(1L, 2L, 2L)))
intercept <- prepare_mcmc_iwls_stan_data(y, c(1, 1), 0.01, w = 1)
linear <- prepare_mcmc_iwls_stan_data(y, c(1, 1), 0.01, X_cont = x,
                                    d = 0L, a_tau = 2, b_tau = 1, w = 1)
stopifnot(intercept$data$P == 1, intercept$data$n_tau == 0,
          linear$data$P == 2, linear$data$n_tau == 1)
rejects(prepare_mcmc_iwls_stan_data(y, c(1, 1), 0.01, freq = -1, w = 1))
rejects(prepare_mcmc_iwls_stan_data(y, c(1, 1), 0.01, X_cont = x, d = 2,
                                  a_tau = 2, b_tau = 1, w = 1))

# Compare CV with the actual sampler, run with no sampling iterations.
# Z_cont is a matrix for the original sampler's cbind call.
for (link_name in c("identity", "logit")) {
  mapping <- if (link_name == "identity") function(x) x else plogis
  args <- base_args
  args$link <- link_name
  args$beta0 <- rep(0, 4)
  args$K_fold <- 3L
  args$seed <- 72
  cv <- do.call(prepare_mcmc_iwls_stan_data, args)
  r_args <- args
  r_args$link <- mapping
  r_args$Z_cont <- z
  r_args$niter <- 0L
  r_args$k_deriv <- args$k
  r_args$c_deriv <- args$c
  invisible(capture.output(r_fit <- do.call(original$MCMC_IWLS, r_args)))
  close(cv$calibration$w, r_fit$calibration$w)
  close(cv$calibration$mean_cv, r_fit$calibration$mean_cv)
}
cat("Data layouts and original-sampler CV agree.\n")

# Keep compilation products outside the project; disable rstan's .rds cache.
build_dir <- tempfile("mcmc-iwls-stan-")
dir.create(build_dir)
old_dir <- getwd()
setwd(build_dir)
sm <- rstan::stan_model(file = file.path(here, "MCMC_IWLS.stan"),
                        auto_write = FALSE, save_dso = FALSE)
setwd(old_dir)

ig_log <- function(v, a, b) sum(dgamma(1 / v, shape = a, rate = b, log = TRUE) - 2 * log(v))
r_target <- function(par, dat) {
  eta <- drop(dat$X %*% par$beta)
  if (dat$use_re) eta <- eta + par$u[dat$id]
  mu <- if (dat$link_id == 1L) eta else plogis(eta)
  out <- -dat$w * sum(dat$freq * original$loss_asymm(dat$y - mu, dat$k, dat$c))
  out <- out + dnorm(par$beta[1], dat$intercept_mean, dat$intercept_sd, log = TRUE)
  ni <- dat$P - 1L - dat$s_gmrf
  if (ni > 0L) out <- out + sum(dnorm(par$beta[1L + seq_len(ni)], 0,
                                     sqrt(par$tau2[dat$tau_idx[seq_len(ni)]]), log = TRUE))
  ng <- dat$n_tau - dat$use_gmrf
  if (ng > 0L) out <- out + ig_log(par$tau2[seq_len(ng)], dat$a_tau[seq_len(ng)], dat$b_tau[seq_len(ng)])
  if (dat$use_gmrf) {
    b <- tail(par$beta, dat$s_gmrf)
    v <- tail(par$tau2, 1)
    out <- out - dat$s_gmrf * log(v)/2 - drop(crossprod(b, dat$K_GMRF %*% b))/(2*v)
    out <- out + dweibull(v, shape = 0.5, scale = dat$theta_GMRF, log = TRUE)
  }
  if (dat$use_re) {
    out <- out + sum(dnorm(par$u, 0, sqrt(par$sigma2u), log = TRUE))
    out <- out + dweibull(par$sigma2u, 0.5, dat$theta_re, log = TRUE)
  }
  out
}

cases <- list(intercept, linear, prepared)
for (link_name in c("identity", "logit")) {
  args <- base_args
  args$link <- link_name
  args$id <- group
  args$theta_re <- 1.7
  args$Z_spat <- spatial
  args$a_tau <- c(2, 3, 4)
  args$b_tau <- c(1, 2, 3)
  args$w <- 2
  cases[[length(cases) + 1L]] <- do.call(prepare_mcmc_iwls_stan_data, args)
  args$Z_spat <- NULL
  args$X_GMRF <- spatial
  args$K_GMRF <- matrix(c(2, -0.3, -0.3, 1), 2)
  args$theta_GMRF <- 0.9
  cases[[length(cases) + 1L]] <- do.call(prepare_mcmc_iwls_stan_data, args)
}

for (i in seq_along(cases)) {
  dat <- cases[[i]]$data
  par <- list(beta = seq_len(dat$P) / 10, tau2 = rep(1.2, dat$n_tau),
              u = seq_len(dat$n_id) / 20, sigma2u = rep(0.8, dat$use_re))
  par <- lapply(par, function(x) array(x, dim = length(x)))
  fit <- rstan::sampling(sm, data = dat, chains = 1, iter = 1, warmup = 0,
                         algorithm = "Fixed_param", init = list(par), seed = 19,
                         refresh = 0)
  stopifnot(length(rstan::extract(fit)$weighted_loss) == 1L)
  lp <- function(p) rstan::log_prob(fit, rstan::unconstrain_pars(fit, p), adjust_transform = FALSE)
  base_lp <- lp(par)
  for (nm in names(par)) {
    for (j in seq_along(par[[nm]])) {
      candidate <- par
      candidate[[nm]][j] <- candidate[[nm]][j] + 0.17
      close(lp(candidate) - base_lp, r_target(candidate, dat) - r_target(par, dat))
    }
  }
  # Compare the original SDP conditional on log variance, including its Jacobian.
  if (dat$use_re) {
    candidate <- par
    candidate$sigma2u <- par$sigma2u * 1.4
    old_sdp <- original$eval_logfullcond_SDP(log(par$sigma2u), diag(dat$n_id), dat$n_id, par$u, dat$theta_re)$logpost
    new_sdp <- original$eval_logfullcond_SDP(log(candidate$sigma2u), diag(dat$n_id), dat$n_id, par$u, dat$theta_re)$logpost
    close(lp(candidate) - base_lp + log(1.4), new_sdp - old_sdp)
  }
  cat("Stan density case", i, "passed.\n")
}
cat("All Stan target checks passed. No posterior convergence claim is made.\n")
