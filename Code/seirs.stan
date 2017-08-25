functions {
  vector derivs(int t,
               vector y,
               real eps_rv,
               real rov,
               real lambda,
               real ro,
               real gamma,
               real eps_dv,
               real delta,
               real cap,
               vector control) {
    
    /**
    * documentation block
    */
    
    vector[12] dydt;
    
    real b;
    real dv;
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

    dv = exp(y[9]);
    
    // Compute transition rates
    // Mosquito to human foi
    foi_vh = control[3] * lambda * y[5] * y[1];

    // Human to mosquito foi
    foi_hv = control[3] * lambda * y[3] * Sv;

    // Infectious humans
    infectious = ro * y[2];
    
    // Infectious mosq
    infect_mosq = rov * y[4];
    
    // Compute derivatives
    /*S*/  dydt[1] = b - b * y[1] - foi_vh + delta * R;
    
    /*E*/ dydt[2] = foi_vh - infectious - b * y[2];
    /*I*/ dydt[3] = infectious - (gamma + b) * y[3];
    
    /*VE*/ dydt[4] = foi_hv - infect_mosq - (control[1] * dv + cap) * y[4];
    /*VI*/ dydt[5] = infect_mosq - (control[1] * dv + cap) * y[5];
    /*VN*/ dydt[6] = control[2] * (y[11] + dv) * y[6] - (control[1] * dv + cap) * y[6];
    
    /*VC*/ dydt[7] = cap * y[6];
    /*cases*/ dydt[8] = infectious;
    
    // Parameter processes
    /*logdv*/ dydt[9] = y[10];
    /*dlogdv*/ dydt[10] = - (pi() / 26) ^ 2 * y[9] + eps_dv;
    /*rv*/ dydt[11] = y[12];
    /*drv*/ dydt[12] = - (pi() / 26) ^ 2 * y[11] + eps_rv;
    
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
  vector[3] control[T + T_pred];
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
  real dv0;
  real rv0;
  real<lower=0> sigmadv;
  real<lower=0> sigmarv;
  vector[T] eps_dv;
  vector[T] eps_rv;
}
transformed parameters {
  vector[12] y0;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;
  vector[T + T_pred] mu_rv;
  vector[T + T_pred] mu_dv;

  // initial conditions
  y0[1] = S0 * (pop - E0 - I0) / pop;
  y0[2] = E0 / pop;
  y0[3] = I0 / pop;
  y0[4] = 0.0;
  y0[5] = 0.0;
  y0[6] = exp(logNv);
  y0[7] = 0.0;
  y0[8] = 0.0;
  y0[9] = dv0;
  y0[10] = 0;
  y0[11] = rv0;
  y0[12] = 0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // Mosquito process
  mu_rv = sigmarv * eps_rv;
  mu_dv = sigmadv * eps_dv;
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  vector[12] state[T * 7 + 1];
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
  dv0 ~ normal(0.39, 0.5);
  rv0 ~ normal(0, 0.5);
  
  // Error component
  eps_rv ~ normal(0, 1);
  eps_dv ~ normal(0, 1);
  
  // Variance parameters
  sigmadv ~ normal(0, 0.2);
  sigmarv ~ normal(0, 0.2);

  // Process model
  
  state[1] = y0;

  for (t in 1:T){
    for(j in 1:7){
      
      state[idx + 1] = state[idx] + 1.0 / 7.0 * derivs(t, 
                                                       state[idx], 
                                                       mu_rv[t], 
                                                       rov[t], 
                                                       lambda, 
                                                       ro, 
                                                       gamma, 
                                                       mu_dv[t], 
                                                       delta, 
                                                       phi_q * tau[t],
                                                       control[t]);
        
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
  vector[12] system[T + T_pred];
  vector[12] state;
  
  state = y0;
  
  // Estimated trajectories
  for (t in 1:T){
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(t, 
                                         state, 
                                         mu_rv[t],
                                         rov[t], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         mu_dv[t], 
                                         delta, 
                                         phi_q * tau[t],
                                         control[t]);
    }
    
    system[t] = state;

    q_hat[t] = neg_binomial_2_rng(state[7] * pop, eta_q);
    y_hat[t] = neg_binomial_2_rng(phi_y * pop * state[8], eta_y);

    state[7] = 0;
    state[8] = 0;
  }
  
  // Predicted trajectory
  for (k in 1:T_pred){

    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(T + k, 
                                         state, 
                                         normal_rng(0, sigmarv),
                                         rov[T + k], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         normal_rng(0, sigmadv), 
                                         delta, 
                                         phi_q * tau[T],
                                         control[T + k]);
    }

    system[T + k] = state;
    
    q_hat[T + k] = neg_binomial_2_rng(state[7] * pop, eta_q);
    y_hat[T + k] = neg_binomial_2_rng(phi_y * pop * state[8], eta_y);
    
    state[7] = 0;
    state[8] = 0;
  }
}
