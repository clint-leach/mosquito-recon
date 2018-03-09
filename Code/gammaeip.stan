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
               real omega,
               real damp, 
               vector control) {
    
    /**
    * documentation block
    */
    
    vector[15] dydt;
    
    real b;
    real dv;
    real Sv;
    real R;
    
    real foi_vh;
    real foi_hv;
    real infectious;

    // Assigning data
    b = 7.0 / (76 * 365);
    
    // Computing mosquito population size
    Sv = y[9]- y[8] - y[7] - y[6] - y[5] - y[4];
    R = 1 - y[1] - y[2] - y[3];

    dv = exp(y[12] + 0.39);
    
    // Compute transition rates
    // Mosquito to human foi
    foi_vh = lambda * y[5] * y[1];

    // Human to mosquito foi
    foi_hv = lambda * y[3] * Sv;

    // Infectious humans
    infectious = ro * y[2];
    
    // Compute derivatives
    /*S*/  dydt[1] = b - b * y[1] - foi_vh + delta * R;
    
    /*E*/ dydt[2] = foi_vh - infectious - b * y[2];
    /*I*/ dydt[3] = infectious - (gamma + b) * y[3];
    
    /*VE1*/ dydt[4] = foi_hv - (control[1] * dv + cap + 4 * rov) * y[4];
    /*VE2*/ dydt[5] = 4 * rov * y[4] - (control[1] * dv + cap + 4 * rov) * y[5];
    /*VE3*/ dydt[6] = 4 * rov * y[5] - (control[1] * dv + cap + 4 * rov) * y[6];
    /*VE4*/ dydt[7] = 4 * rov * y[6] - (control[1] * dv + cap + 4 * rov) * y[7];
    /*VI*/ dydt[8] = 4 * rov * y[7] - (control[1] * dv + cap) * y[8];
    /*VN*/ dydt[9] = control[2] * (y[14] + dv) * y[9] - (control[1] * dv + cap) * y[9];
    
    /*VC*/ dydt[10] = cap * y[9];
    /*cases*/ dydt[11] = infectious;
    
    // Parameter processes
    /*logdv*/ dydt[12] = y[13];
    /*dlogdv*/ dydt[13] = -2 * damp * omega * y[13] - omega ^ 2 * y[12] + eps_dv;
    /*rv*/ dydt[14] = y[15];
    /*drv*/ dydt[15] = -2 * damp * omega * y[15] - omega ^ 2 * y[14] + eps_rv;
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> steps;
  int y[T];
  int q[T];
  vector[T] tau;
  vector[T] rov;
  vector[3] control[T];
  int pop;
}
transformed data {
  real lambda = 4.87;
  real phi_y = 1.0 / 12.0;
  real damp = 0.0;
  real omega = (pi() / 26) / sqrt(1 - damp ^ 2);
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
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;
  vector[T] y_hat;
  vector[T] q_hat;
  vector[15] state[T + 1];
  
  // initial conditions
  state[1, 1] = S0 * (pop - E0 - I0) / pop;
  state[1, 2] = E0 / pop;
  state[1, 3] = I0 / pop;
  state[1, 4] = 0.0;
  state[1, 5] = 0.0;
  state[1, 6] = 0.0;
  state[1, 7] = 0.0;
  state[1, 8] = 0.0;
  state[1, 9] = exp(logNv);
  state[1, 10] = 0.0;
  state[1, 11] = 0.0;
  state[1, 12] = dv0;
  state[1, 13] = 0;
  state[1, 14] = rv0;
  state[1, 15] = 0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // Process model
  for (t in 2:(T + 1)){
    
    state[t] = state[t - 1];
    
    state[t, 4] = state[t, 4] * control[t - 1, 3];
    state[t, 5] = state[t, 5] * control[t - 1, 3];
    state[t, 6] = state[t, 6] * control[t - 1, 3];
    state[t, 7] = state[t, 7] * control[t - 1, 3];
    state[t, 8] = state[t, 8] * control[t - 1, 3];
    state[t, 9] = state[t, 9] * control[t - 1, 3];
    
    for(j in 1:steps){
      
      state[t] = state[t] + 1.0 / steps * derivs(t, 
                                                 state[t], 
                                                 sigmarv * eps_rv[t - 1], 
                                                 rov[t - 1], 
                                                 lambda, 
                                                 ro, 
                                                 gamma, 
                                                 sigmadv * eps_dv[t - 1], 
                                                 delta, 
                                                 phi_q * tau[t - 1], 
                                                 omega, 
                                                 damp, 
                                                 control[t - 1]);
    }
                                                      
    q_hat[t - 1] = state[t, 10] - state[t - 1, 10];
    y_hat[t - 1] = state[t, 11] - state[t - 1, 11];

  }
}
model {

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
  
  // Mean and initial values
  dv0 ~ normal(0, 0.5);
  rv0 ~ normal(0, 0.5);
  
  // Error component
  eps_rv ~ normal(0, 1);
  eps_dv ~ normal(0, 1);
  
  // Variance parameters
  sigmadv ~ normal(0, 0.2);
  sigmarv ~ normal(0, 0.2);

  
  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {
  vector[T] y_meas;
  vector[T] q_meas;

  // Estimated trajectories
  for (t in 1:T){
    
    q_meas[t] = neg_binomial_2_rng(q_hat[t] * pop, eta_q);
    y_meas[t] = neg_binomial_2_rng(phi_y * y_hat[t] * pop, eta_y);
    
  }
}
