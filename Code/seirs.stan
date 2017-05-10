functions {
  vector derivs(int t,
               vector y,
               real rv,
               real rov,
               real alpha,
               real ro,
               real gamma,
               real dv,
               real delta,
               real lambda,
               int pop) {
    
    /**
    * documentation block
    */
    
    vector[8] dydt;
    
    real b;
    real Sv;
    real R;
    
    real foi_vh;
    real foi_hv;
    real infectious;
    real infect_mosq;
    
    // Assigning data
    b = 7.0 / (76 * 365);

    // Computing mosquito population size
    Sv = y[6] - y[4] - y[5];
    R = 1 - y[1] - y[2] - y[3];

    // Compute transition rates
    // Mosquito to human foi
    foi_vh = alpha * y[5] * y[1];

    // Human to mosquito foi
    foi_hv = alpha * y[3] * Sv;

    // Infectious humans
    infectious = ro * y[2];
    
    // Infectious mosq
    infect_mosq = rov * y[4];
    
    // Compute derivatives
    /*S*/  dydt[1] = b - b * y[1] - foi_vh + delta * R;
    
    /*E*/ dydt[2] = foi_vh - infectious - b * y[2];
    /*I*/ dydt[3] = infectious - (gamma + b) * y[3];

    /*VE*/ dydt[4] = foi_hv - infect_mosq - dv * y[4];
    /*VI*/ dydt[5] = infect_mosq - dv * y[5];
    /*VN*/ dydt[6] = rv * y[6] - lambda * y[6];
    
    /*VC*/ dydt[7] = lambda * y[6];
    /*cases*/ dydt[8] = infectious;
    
    return dydt;
  }
  // Softmax function
  vector softmax_id(vector alpha) {
    vector[num_elements(alpha) + 1] alphac;
      for (k in 1:num_elements(alpha))
        alphac[k] = alpha[k];
    alphac[num_elements(alphac)] = 0;
    return softmax(alphac);
  }
}
data {
  int<lower=1> T;
  int y[T];
  int q[T];
  vector[T] tau;
  real rov[T];
  matrix[T, T] sinmat;
  int pop;
}
transformed data {
  vector[T] mu = rep_vector(0.0, T);
  real alpha = 4.87;
  real phi_y = 1.0 / 12.0;
}
parameters {
  vector[3] p0_raw;                    // untransformed initial conditions
  real<upper=0> log_phi_q;             // log per-trap capture rate
  real<lower=0> eta_inv_y;             // overdispersion of case reports
  real<lower=0> eta_inv_q;             // overdispersion of mosquito capture
  real<lower=0> ro_c;                  // human latenet period
  real<lower=0> gamma_c;               // human infectious period
  real<lower=0> delta_c;               // cross-immune period
  real logNv0;
  real<lower=0> hsq_r;
  real<lower=0> omegasq_r;
  real<lower=0> hsq_d;
  real<lower=0> omegasq_d;
  vector[T] z_r;
  vector[T] z_d;
}
transformed parameters {
  vector[4] p0;
  vector[8] y0;
  vector[T] dv;
  vector[T] rv;
  vector[T] lambda;
  matrix[T, T] Sigma_r;
  matrix[T, T] L_r;
  matrix[T, T] Sigma_d;
  matrix[T, T] L_d;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;

  // initial conditions
  
  p0 = softmax_id(p0_raw);

  for(i in 1:3){
    y0[i] = p0[i];
  }
  y0[4] = 0.0;
  y0[5] = 0.0;
  y0[6] = exp(logNv0);
  y0[7] = 0.0;
  y0[8] = 0.0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // mosquito capture rate
  lambda = exp(log_phi_q) * tau;
  
  // GP covariance matrix
  for (i in 1:(T - 1)){
    for (j in (i + 1):T){
      Sigma_r[i, j] = hsq_r * exp(-1 / (2 * omegasq_r) * sinmat[i, j]);
      Sigma_r[j, i] = Sigma_r[i, j];
      
      Sigma_d[i, j] = hsq_d * exp(-1 / (2 * omegasq_d) * sinmat[i, j]);
      Sigma_d[j, i] = Sigma_d[i, j];
    }
  }
  
  for (i in 1:T){
      Sigma_r[i, i] = hsq_r + 0.001;
      Sigma_d[i, i] = hsq_d + 0.001;
  }
  
  L_r = cholesky_decompose(Sigma_r);
  L_d = cholesky_decompose(Sigma_d);

  rv = mu + L_r * z_r;
  dv = 1.47 * exp(L_d * z_d);
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  vector[8] state[T * 7 + 1];
  int idx = 1;
  
  // Priors
  
  // Initial conditions
  p0_raw[1] ~ normal(-0.4, 0.2);
  p0_raw[2] ~ normal(-9, 0.6);
  p0_raw[3] ~ normal(-9, 0.6);

  // Measurement models
  log_phi_q ~ normal(-13, 0.5);
  eta_inv_y ~ normal(0, 5);
  eta_inv_q ~ normal(0, 5);
  
  // Epidemiological parameters
  ro_c ~ gamma(8.3, 8.3); 
  gamma_c ~ gamma(100, 100);
  delta_c ~ gamma(10, 10);
  
  // Initial mosquito pop size
  logNv0 ~ normal(0.7, 0.5);
  
  // GP stuff
  hsq_r ~ normal(0, 5);
  omegasq_r ~ normal(0, 5);
  hsq_d ~ normal(0, 5);
  omegasq_d ~ normal(0, 5);

  z_r ~ normal(0, 1);
  z_d ~ normal(0, 1);

  // Process model
  
  state[1] = y0;

  for (t in 1:T){
    for(j in 1:7){
      
      state[idx + 1] = state[idx] + 1.0 / 7.0 *
        derivs(t, state[idx], rv[t], rov[t], alpha, ro, gamma, dv[t], delta, lambda[t], pop);
        
      idx = idx + 1;
    }
    
    q_hat[t] = state[idx, 7] - state[idx - 7, 7];
    y_hat[t] = state[idx, 8] - state[idx - 7, 8];
  }
  
  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {

  vector[T] y_hat;
  vector[T] q_hat;
  vector[8] state;
  
  state = y0;

  // Estimated trajectories
  for (t in 1:T){
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * 
        derivs(t, state, rv[t], rov[t], alpha, ro, gamma, dv[t], delta, lambda[t], pop);
    }

    q_hat[t] = neg_binomial_2_rng(state[7] * pop, eta_q);
    y_hat[t] = neg_binomial_2_rng(phi_y * pop * state[8], eta_y);

    state[7] = 0;
    state[8] = 0;
    
  }
}
