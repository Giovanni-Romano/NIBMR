functions {
  // Safe sqrt(k^2 r^2 + c)
  real s_val(real r, real k, real c) {
    return sqrt(square(k) * square(r) + c);
  }
}

data {
  int<lower=1> n; // number of obs.
  int<lower=1> J; // number of covariates
  array[J] int<lower=1> q;  // q[j] = number of basis for covariate j
  int<lower=1> p; // number of spline basis (\sum_{j = 1}^J q_j)
  matrix[n, p] B; // basis design matrix
  // matrix[n, p] Z; // basis design matrix after orthogonalization
  vector[n] y; // response
  
  real<lower=0> k; // const. for Nest.Int.Loss
  real<lower=0> c; // const. for Nest.Int.Loss
}

transformed data {
  array[J] int start;
  array[J] int end;
  int pos;
  
  pos = 1;
  for (j in 1:J) {
    start[j] = pos;
    end[j]   = pos + q[j] - 1;
    pos      = end[j] + 1;
  }
  
  if (pos - 1 != p) reject("p must equal sum(q). Found p = ", p, ", sum(q) = ", pos - 1);
  
}

parameters {
  vector[p] beta;
  // vector<lower=0.0>[J] lambda;
}

model {
  // ---- compute intermediate quantities ----
  vector[n] eta = B * beta;
  vector[n] r   = y - eta;
  
  // Objects needed to compute m-bar
  vector[n] s;
  vector[n] t; // weights t_i
  
  // M has rows m_i' = -t_i * B_[i, .]
  matrix[n, p] M;
  
  vector[p] mbar;
  matrix[p, p] W;
  
  for (i in 1:n) {
    s[i] = s_val(r[i], k, c);
    // t_i = exp(-s_i) * k^2 * r_i / s_i
    t[i] = exp(-s[i]) * square(k) * r[i] / s[i];
    M[i] = -t[i] * B[i];
  }
  
  # // Row i: -t[i] * B[i,]
  #M = diag_pre_multiply(-t, B);
  
  // mbar = mean of rows of M (as a column vector)
  mbar = (M' * rep_vector(1.0, n)) / n;
  
  // Center rows of M
  {
    matrix[n, p] Mc;
    Mc = M - rep_matrix(mbar', n);
    W  = (Mc' * Mc) / n; // VarCov matrix of mbar is VarCov of M[i] divided by n
  }
  
  // ---- prior on the coefficients ----
  for (j in 1:J) {
    vector[q[j]] beta_j = beta[start[j]:end[j]];
    // target += normal_lpdf(beta_j | 0, inv_sqrt(lambda[j]));
    target += normal_lpdf(beta_j | 0, 1);
  }
  // ---- prior on the smoothing parameters ----
  // target += gamma_lpdf(lambda | 1, 5*1e-3);
  
  // ---- calibrated GBI pseudo-likelihood ----
  {
    // log|W| = 2 * sum(log(diag(L)))
    real log_det_W = log_determinant(W);
    
    // Q = 0.5 mbar' W^{-1} mbar
    real Q = 0.5 * quad_form(inverse(W), mbar);
    
    // Add pseudo-likelihood to prior
    target += -0.5 * log_det_W - n * Q;
  }
}

generated quantities {
  vector[n] r;
  
  {
    vector[n] eta = B * beta;
    r = y - eta;
  }
}
