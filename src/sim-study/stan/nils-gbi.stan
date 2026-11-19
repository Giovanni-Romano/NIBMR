// M0 / M2 in ../02-worker.R: identity link, DR nonlinear blocks, optional
// MMR spatial block. Its coefficient penalty is I after reparameterization.
functions {
  real nil_loss(real residual, vector k, real c) {
    real kr = k[residual < 0 ? 1 : 2] * residual;
    real u = square(kr);
    real a = sqrt(u + c) - sqrt(c);
    // Same loss as R: rationalize the small difference near zero, then use
    // expm1 for 1-exp(-a). This changes numerical accuracy, not the target.
    if (u <= c) a = u / (sqrt(u + c) + sqrt(c));
    return -expm1(-a);
  }
}

data {
  int<lower=1> n;
  int<lower=3> P;
  matrix[n, P] X; // (1, scaled x1/x2, DR blocks for M0, MMR spatial block for X3)
  vector[n] y;
  vector<lower=0>[2] k;
  real<lower=0> c;
  real<lower=0> w;
  int<lower=2> n_tau;
  int<lower=2, upper=n_tau> n_ig;
  array[P - 1] int<lower=1, upper=n_tau> tau_idx;
  vector<lower=0>[n_ig] a_tau;
  vector<lower=0>[n_ig] b_tau;
  int<lower=0, upper=P - 3> s; // Number of penalized MMR coefficients, NOT regions.
  real<lower=0> theta_GMRF;
  real intercept_mean;
  real<lower=0> intercept_sd;
}

transformed data {
  if (c <= 0 || w <= 0 || intercept_sd <= 0)
    reject("c, w and intercept_sd must be strictly positive.");
  if (min(a_tau) <= 0 || min(b_tau) <= 0)
    reject("Inverse-gamma shapes and scales must be positive.");
  if (n_tau - n_ig != (s > 0 ? 1 : 0))
    reject("Only the optional final MMR block has an SDP variance prior.");
  if (s > 0 && theta_GMRF <= 0) reject("theta_GMRF must be positive.");
  for (i in 1:n) {
    if (X[i, 1] != 1) reject("The first X column must be ones.");
  }
  for (j in 1:(P - 1)) {
    if (s > 0 && j > P - 1 - s) {
      if (tau_idx[j] != n_tau) reject("MMR coefficients must share the final tau.");
    } else if (tau_idx[j] > n_ig) {
      reject("Only MMR coefficients may use the SDP spatial scale.");
    }
  }
}

parameters {
  vector[P] beta;
  vector<lower=0>[n_tau] tau2;
}

transformed parameters {
  vector<lower=0>[n_tau] tau = sqrt(tau2);
}

model {
  vector[n] residual = y - X * beta;
  beta[1] ~ normal(intercept_mean, intercept_sd);
  // The MMR decomposition gives K_GMRF = I, including in the R sampler.
  tail(beta, P - 1) ~ normal(0, tau[tau_idx]);
  head(tau2, n_ig) ~ inv_gamma(a_tau, b_tau);
  if (s > 0) tau2[n_tau] ~ weibull(0.5, theta_GMRF);
  for (i in 1:n) target += -w * nil_loss(residual[i], k, c);
}

generated quantities {
  real total_loss = 0;
  real log_loss_kernel;
  {
    vector[n] residual = y - X * beta;
    for (i in 1:n) total_loss += nil_loss(residual[i], k, c);
  }
  log_loss_kernel = -w * total_loss;
}
