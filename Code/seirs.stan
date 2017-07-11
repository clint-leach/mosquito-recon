functions {
  vector derivs(int t,
               vector y,
               real rv,
               real rov,
               real lambda,
               real ro,
               real gamma,
               real dv,
               real delta,
               real cap) {
    
    /**
    * documentation block
    */
    
    vector[8] dydt;
    
    real b;
    real loss;
    real Sv;
    real R;
    
    real foi_vh;
    real foi_hv;
    real infectious;
    real infect_mosq;
    
    // Assigning data
    b = 7.0 / (76 * 365);

    // Computing mosquito population size
    Sv = y[6] - y[5] - y[4];
    R = 1 - y[1] - y[2] - y[3];

    loss = dv + cap;
    
    // Compute transition rates
    // Mosquito to human foi
    foi_vh = lambda * y[5] * y[1];

    // Human to mosquito foi
    foi_hv = lambda * y[3] * Sv;

    // Infectious humans
    infectious = ro * y[2];
    
    // Infectious mosq
    infect_mosq = rov * y[4];
    
    // Compute derivatives
    /*S*/  dydt[1] = b - b * y[1] - foi_vh + delta * R;
    
    /*E*/ dydt[2] = foi_vh - infectious - b * y[2];
    /*I*/ dydt[3] = infectious - (gamma + b) * y[3];
    
    /*VE*/ dydt[4] = foi_hv - infect_mosq - loss * y[4];
    /*VI*/ dydt[5] = infect_mosq - loss * y[5];
    /*VN*/ dydt[6] = rv * y[6];
    
    /*VC*/ dydt[7] = cap * y[6];
    /*cases*/ dydt[8] = infectious;
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=0> T_pred;
  int y[T];
  int q[T];
  vector[T] tau;
  vector[T + T_pred] rov;
  int pop;
}
transformed data {
  real lambda = 4.87;
  real phi_y = 1.0 / 12.0;
}
parameters {
  real<lower=0,upper=1> S0;            // untransformed initial conditions
  real<lower=0> E0;
  real<lower=0> I0;
  real<upper=0> log_phi_q;             // log per-trap capture rate
  real<lower=0> eta_inv_y;             // overdispersion of case reports
  real<lower=0> eta_inv_q;             // overdispersion of mosquito capture
  real<lower=0> ro_c;                  // human latenet period
  real<lower=0> gamma_c;               // human infectious period
  real<lower=0> delta_c;               // cross-immune period
  real logNv;                         // initial mosquito population size
  real psi0;
  real rv0;
  real<lower=-1,upper=1> alpha_rv;
  real<lower=-1,upper=1> alpha_psi;
  real<lower=-1,upper=1> theta_rv;
  real<lower=-1,upper=1> theta_psi;
  real<lower=0> sigma0_psi;
  real<lower=0> sigma0_rv;
  real<lower=0> sigmapsi;
  real<lower=0> sigmarv;
  vector[T - 1] eps_psi;
  vector[T - 1] eps_rv;
}
transformed parameters {
  vector[8] y0;
  vector[T] dv;
  vector[T] rv;
  vector[T] psi_raw;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;

  // initial conditions
  y0[1] = S0 * (pop - E0 - I0) / pop;
  y0[2] = E0 / pop;
  y0[3] = I0 / pop;
  y0[4] = 0.0;
  y0[5] = 0.0;
  y0[6] = exp(logNv);
  y0[7] = 0.0;
  y0[8] = 0.0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // mosquito demographic parameters
  {
    vector[T - 1] z_psi;
    vector[T - 1] z_rv;
    
    z_psi[1:52] = sigma0_psi * eps_psi[1:52];
    z_rv[1:52] = sigma0_rv * eps_rv[1:52];

    for(i in 53:(T - 1)){
      z_psi[i] = theta_psi * z_psi[i - 52] + sigmapsi * eps_psi[i];      
      z_rv[i] = theta_rv * z_rv[i - 52] + sigmarv * eps_rv[i];
    }
    
    psi_raw[1] = psi0;
    rv[1] = rv0;
    
    for(i in 2:T){
      psi_raw[i] = alpha_psi * psi_raw[i - 1] + z_psi[i - 1];
      rv[i] = alpha_rv * rv[i - 1] + z_rv[i - 1];
    }
  }
  
  dv = 3.56 * rov[1:T] .* exp(psi_raw);
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  vector[8] state[T * 7 + 1];
  int idx = 1;
  
  // Priors
  
  // Initial conditions
  S0 ~ beta(4, 6);
  E0 ~ gamma(10, 0.1);
  I0 ~ gamma(6, 0.1);

  // Measurement models
  log_phi_q ~ normal(-13, 0.5);
  eta_inv_y ~ normal(0, 5);
  eta_inv_q ~ normal(0, 5);
  
  // Epidemiological parameters
  ro_c ~ gamma(8.3, 8.3); 
  gamma_c ~ gamma(100, 100);
  delta_c ~ gamma(10, 10);
  
  // Initial mosquito pop size
  logNv ~ normal(0.7, 0.3);
  
  // Mosquito demographic series
  
  // Initial values
  psi0 ~ normal(0, 0.5);
  rv0 ~ normal(0, 0.5);
  
  // SAR parameters
  alpha_rv ~ normal(0, 1);
  alpha_psi ~ normal(0, 1);
  
  theta_rv ~ normal(0, 1);
  theta_psi ~ normal(0, 1);
  
  // Error component
  eps_rv ~ normal(0, 1);
  eps_psi ~ normal(0, 1);
  
  // Variance parameters
  sigma0_psi ~ normal(0, 0.1);
  sigma0_rv ~ normal(0, 0.1);
  
  sigmapsi ~ normal(0, 0.1);
  sigmarv ~ normal(0, 0.1);
  
  // Process model
  
  state[1] = y0;

  for (t in 1:T){
    for(j in 1:7){
      
      state[idx + 1] = state[idx] + 1.0 / 7.0 * derivs(t, 
                                                       state[idx], 
                                                       rv[t], 
                                                       rov[t], 
                                                       lambda, 
                                                       ro, 
                                                       gamma, 
                                                       dv[t], 
                                                       delta, 
                                                       phi_q * tau[t]);
        
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

  vector[T + T_pred] y_hat;
  vector[T + T_pred] q_hat;
  vector[T + T_pred] psi;
  vector[T + T_pred] risk;
  vector[T + T_pred] psi_raw_full;
  vector[T + T_pred] rv_full;
  // vector[8] system[T + T_pred];
  vector[8] state;
  
  state = y0;
  
  // Estimated trajectories
  for (t in 1:T){
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(t, 
                                         state, 
                                         rv[t],
                                         rov[t], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         dv[t], 
                                         delta, 
                                         phi_q * tau[t]);
    }
    
    // system[t] = state;

    q_hat[t] = neg_binomial_2_rng(state[7] * pop, eta_q);
    y_hat[t] = neg_binomial_2_rng(phi_y * pop * state[8], eta_y);

    state[7] = 0;
    state[8] = 0;
    
    psi_raw_full[t] = psi_raw[t];
    rv_full[t] = rv[t];
    
    psi[t] = inv_logit(-1.27 + psi_raw[t]);
    risk[t] = (state[6] - state[5] - state[4]) * psi[t];
  }
  
  // Predicted trajectory
  for (k in 1:T_pred){
    
    psi_raw_full[T + k] = psi_raw_full[T + k - 1] + 
                          normal_rng(0, sigmapsi); 
                          
    rv_full[T + k] = rv_full[T + k - 1] + 
                     normal_rng(0, sigmarv); 

    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(T + k, 
                                         state, 
                                         rv_full[T + k],
                                         rov[T + k], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         rov[T + k] * exp(1.27 - psi_raw_full[T + k]), 
                                         delta, 
                                         phi_q * tau[T]);
    }

    // system[T + k] = state;
    
    q_hat[T + k] = neg_binomial_2_rng(state[7] * pop, eta_q);
    y_hat[T + k] = neg_binomial_2_rng(phi_y * pop * state[8], eta_y);
    
    state[7] = 0;
    state[8] = 0;
    
    psi[T + k] = inv_logit(-1.27 + psi_raw_full[T + k]);
    risk[T + k] = (state[6] - state[5] - state[4]) * psi[T + k];

  }
}
