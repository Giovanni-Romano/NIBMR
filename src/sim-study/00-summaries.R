# Compact archives: empirical posterior summaries, not posterior draws.
# Ten fixed probabilities include the median, 95% interval and upper tail.
POSTERIOR_PROBS <- c(.025, .05, .10, .25, .50, .75, .90, .95, .975, .99)

parameter_summary <- function(draws, chain = rep(1L, nrow(draws)),
                              probs = POSTERIOR_PROBS,
                              mean_exists = rep(NA, ncol(draws)),
                              sd_exists = rep(NA, ncol(draws))) {
  draws <- as.matrix(draws)
  ids <- unique(chain)
  columns <- lapply(seq_len(ncol(draws)), function(j) {
    x <- draws[, j]
    stopifnot(all(is.finite(x)))
    chains <- vapply(ids, function(id) x[chain == id], numeric(sum(chain == ids[1L])))
    diagnostics <- if (nrow(chains) >= 4L) c(
      ess_bulk = posterior::ess_bulk(chains), ess_tail = posterior::ess_tail(chains),
      rhat = posterior::rhat(chains),
      mcse_mean = if (isFALSE(sd_exists[j])) NA_real_ else posterior::mcse_mean(chains)) else
      c(ess_bulk = NA_real_, ess_tail = NA_real_, rhat = NA_real_, mcse_mean = NA_real_)
    # mean/sd describe these finite draws. With bounded NIL loss, some prior
    # moments remain infinite in the posterior; flags prevent misinterpreting
    # finite empirical summaries as estimates of finite theoretical moments.
    c(mean = mean(x), sd = sd(x), setNames(quantile(x, probs, names = FALSE),
      paste0("q", format(100 * probs, trim = TRUE))), diagnostics)
  })
  out <- data.frame(parameter = colnames(draws), do.call(rbind, columns), row.names = NULL,
                    check.names = FALSE)
  out$mean_exists <- mean_exists; out$sd_exists <- sd_exists
  out
}

chain_description <- function(chain) list(
  n_chains = length(unique(chain)), n_draws = length(chain),
  rhat_type = if (length(unique(chain)) == 1L) "single-chain split; no between-chain check" else
    "rank-normalized split/folded R-hat")

stan_diagnostics <- function(samples, max_treedepth = 10L, thin = 1L) {
  if (is.null(samples)) return(NULL)
  if (is.null(max_treedepth)) max_treedepth <- 10L
  do.call(rbind, lapply(seq_along(samples), function(j) {
    x <- samples[[j]]
    energy_var <- var(x[, "energy__"])
    data.frame(chain = j, iterations = nrow(x), stored_thin = thin,
      mean_accept_stat = mean(x[, "accept_stat__"]),
      divergences = sum(x[, "divergent__"]),
      treedepth_hits = sum(x[, "treedepth__"] >= max_treedepth),
      max_treedepth = max_treedepth, mean_leapfrog = mean(x[, "n_leapfrog__"]),
      ebfmi = if (is.finite(energy_var) && energy_var > 0)
        mean(diff(x[, "energy__"])^2) / energy_var else NA_real_)
  }))
}

# Called while each M4 candidate's Stan object is still available. Nuisance
# parameters and diagnostics are retained for losing candidates as well.
compact_gud_candidate <- function(draws, distribution_draws, stanfit, X,
                                  probs = POSTERIOR_PROBS) {
  B <- draws
  colnames(B) <- sub("^\\(Intercept\\)$", "Intercept", colnames(X))
  B <- cbind(B, distribution_draws)
  chain <- rep(1L, nrow(B)) # fit_gud currently runs one chain per candidate.
  depth <- stanfit@stan_args[[1L]]$control$max_treedepth
  if (is.null(depth)) depth <- 10L
  list(type = "posterior", parameters = parameter_summary(B, chain, probs),
       chains = chain_description(chain),
       sampler = stan_diagnostics(rstan::get_sampler_params(stanfit, inc_warmup = FALSE),
                                  depth, stanfit@sim$thin))
}

# Moment conditions follow the M0/M2/M5 proper joint prior. Since the NIL
# kernel is bounded above and below by positive constants, it cannot change
# whether a nonnegative prior moment is finite. IG(a,b): E[v^r] exists iff r<a.
nils_moment_flags <- function(fit, design) {
  alpha <- unname(fit$tau_prior$a_tau)
  p <- ncol(design$X_cont); q <- design$q
  if (is.null(q)) q <- 0L
  means <- sds <- rep(TRUE, ncol(fit$draws)) # Gaussian intercept/spatial prior.
  linear <- 1L + seq_len(q + p)
  means[linear] <- alpha[seq_len(q + p)] > .5
  sds[linear] <- alpha[seq_len(q + p)] > 1
  if (sum(design$d) > 0L) {
    nonlinear <- 1L + q + p + seq_len(sum(design$d))
    a <- rep(alpha[q + p + seq_len(p)], design$d)
    means[nonlinear] <- a > .5; sds[nonlinear] <- a > 1
  }
  if (design$s > 0L) alpha <- c(alpha, Inf) # SDP spatial variance has all moments.
  list(mean = c(means, alpha > .5, alpha > 1),
       sd = c(sds, alpha > 1, alpha > 2))
}

compact_fit <- function(fit, model, design, probs = POSTERIOR_PROBS) {
  if (model == "m3") return(list(type = "frequentist", estimates = fit$beta,
    optimization = fit[intersect(c("method", "value", "convergence", "message", "all",
                                   "optim_settings", "loss"), names(fit))]))
  if (model == "m4" && !is.null(fit$candidates)) return(list(type = "posterior-selection",
    selected = fit$selected, elpd = fit$elpd, candidates = fit$candidates))
  if (is.null(fit$draws)) return(list(type = "unavailable", selected = fit$selected))
  B <- fit$draws; chain <- rep(1L, nrow(B))
  moments <- NULL; sampler <- NULL; native_family <- NULL
  if (model %in% c("m0", "m2", "m5")) {
    if (!is.null(fit$tau)) {
      tau <- fit$tau
      B <- cbind(B, setNames(as.data.frame(tau), paste0("tau.", colnames(tau))),
                 setNames(as.data.frame(tau^2), paste0("tau2.", colnames(tau))))
      B <- as.matrix(B)
      if (!is.null(fit$tau_prior)) moments <- nils_moment_flags(fit, design)
    }
    if (!is.null(fit$chain)) chain <- fit$chain
    if (model == "m5") {
      sampler <- stan_diagnostics(fit$sampler_diagnostics, fit$stan_settings$control$max_treedepth)
    } else {
      acceptance <- fit$acc
      if (!is.null(acceptance$beta)) names(acceptance$beta) <- c("Intercept",
        colnames(design$X_cont), if (sum(design$d) > 0L) paste0(colnames(design$X_cont), ".nonlinear"),
        if (design$s > 0L) "GMRF")
      sampler <- list(acceptance = acceptance, acceptance_period = "all iterations, including burn-in")
    }
  } else if (model == "m1" && !is.null(fit$fit)) {
    raw <- as.matrix(fit$fit$samples)
    if (inherits(fit$fit$samples, "mcmc.list"))
      chain <- rep(seq_along(fit$fit$samples), vapply(fit$fit$samples, nrow, integer(1)))
    diagnostic <- grepl("\\.(alpha|accepted|edf)$|^logLik$", colnames(raw))
    B <- raw[, !diagnostic, drop = FALSE]
    # Native sigma coefficients use the log link. In this bamlss ALD family,
    # sigma^2 is the density's scale, not the error variance.
    if ("sigma.p.(Intercept)" %in% colnames(B)) B <- cbind(B,
      sigma = exp(B[, "sigma.p.(Intercept)"]),
      ald_scale = exp(2 * B[, "sigma.p.(Intercept)"]))
    accepted <- grep("\\.accepted$", colnames(raw))
    alpha <- grep("\\.alpha$", colnames(raw))
    sampler <- list(acceptance = colMeans(raw[, accepted, drop = FALSE]),
      native_alpha_mean = colMeans(raw[, alpha, drop = FALSE]),
      acceptance_period = "retained post-burn-in iterations")
    native_family <- list(family = "ALD", tau = .5,
      eps = get("eps", environment(fit$fit$family$d)), links = fit$fit$family$links)
  } else if (model == "m4") B <- cbind(B, fit$distribution_draws)
  args <- list(draws = B, chain = chain, probs = probs)
  if (!is.null(moments)) {
    args$mean_exists <- moments$mean; args$sd_exists <- moments$sd
  }
  metadata <- fit[intersect(c("calibration", "loss", "intercept_prior", "tau_prior",
                             "stan_settings", "schedule", "selected", "elpd"), names(fit))]
  if (!is.null(native_family)) metadata$family <- native_family
  list(type = "posterior", parameters = do.call(parameter_summary, args),
    chains = chain_description(chain), sampler = sampler,
    metadata = metadata)
}

compact_models <- function(s, probs = POSTERIOR_PROBS) {
  setNames(lapply(names(s$fits), function(model)
    compact_fit(s$fits[[model]], model, s$designs[[model]], probs)), names(s$fits))
}
