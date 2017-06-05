functions {
  vector derivs(int t,
               vector y,
               real alpha0,
               real alpha1,
               real alpha2,
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
    
    vector[21] dydt;
    
    real b;
    real bv;
    real Sv;
    real R;
    
    real foi_vh_1;
    real foi_vh_2;
    
    real foi_hv_1;
    real foi_hv_2;
    
    // Assigning data
    b = 7.0 / (76 * 365);
    bv = exp(alpha0 + alpha1 * sin(2 * pi() * t / 52) + alpha2 * cos(2 * pi() * t / 52));

    // Computing mosquito population size
    Sv = y[15] - y[16] - y[17] - y[18] - y[19];
    
    // Compute forces of infection
    // Mosquito to human
    foi_vh_1 = lambda * y[17];
    foi_vh_2 = lambda * y[19];
    
    // Human to mosquito
    foi_hv_1 = lambda * (y[3] + y[13] + 1e-6);
    foi_hv_2 = lambda * (y[7] + y[9] + 1e-6);
    
    // Compute derivatives
    /*S0*/  dydt[1] = b - b * y[1] - (foi_vh_1 + foi_vh_2) * y[1];
    
    /*E1a*/ dydt[2] = foi_vh_1 * y[1] - (ro + b) * y[2];
    /*I1a*/ dydt[3] = ro * y[2] - (gamma + b) * y[3];
    /*C1*/  dydt[4] = gamma * y[3] - (delta + b) * y[4];
    /*S2*/  dydt[5] = delta * y[4] - b * y[5] - foi_vh_2 * y[5];
    /*E2b*/ dydt[6] = foi_vh_2 * y[5] - (ro + b) * y[6];
    /*I2b*/ dydt[7] = ro * y[6] - (gamma + b) * y[7];
    
    /*E2a*/ dydt[8] = foi_vh_2 * y[1] - (ro + b) * y[8];
    /*I2a*/ dydt[9] = ro * y[8]  - (gamma + b) * y[9];
    /*C2*/  dydt[10] = gamma * y[9] - (delta + b) * y[10];
    /*S1*/  dydt[11] = delta * y[10] - b * y[11] - foi_vh_1 * y[11];
    /*E1b*/ dydt[12] = foi_vh_1 * y[11] - (ro + b) * y[12];
    /*I1b*/ dydt[13] = ro * y[12] - (gamma + b) * y[13];
    
    /*R*/   dydt[14] = gamma * (y[7] + y[13]) - b * y[14];
    
    /*VN*/  dydt[15] = bv - dv * y[15] - cap * y[15];
    /*VE1*/ dydt[16] = foi_hv_1 * Sv - (rov + dv + cap) * y[16];
    /*VI1*/ dydt[17] = rov * y[16] - (dv + cap) * y[17];
    /*VE2*/ dydt[18] = foi_hv_2 * Sv - (rov + dv + cap) * y[18];
    /*VI2*/ dydt[19] = rov * y[18] - (dv + cap) * y[19];
    
    /*cases*/ dydt[20] = ro * (y[2] + y[8] + y[6] + y[12]);
    /*VC*/ dydt[21] = cap * y[15];
    
    return dydt;
  }
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
  int<lower=0> T_pred;
  int<lower=1> D;
  int y[T];
  int q[T];
  vector[T] tau;
  real rov[T + T_pred];
  matrix[T + T_pred, D] covars;
  int pop;
}
transformed data {
  real lambda = 4.87;
  real phi_y = 1.0 / 12.0;
}
parameters {
  vector[5] p0_raw;                   // untransformed initial conditions
  real<upper=0> log_phi_q;             // log per-trap capture rate
  real<lower=0> eta_inv_y;             // overdispersion of case reports
  real<lower=0> eta_inv_q;             // overdispersion of mosquito capture
  real<lower=0> ro_c;                  // human latenet period
  real<lower=0> gamma_c;               // human infectious period
  real<lower=0> delta_c;               // cross-immune period
  real logNv0;                         // initial mosquito population size
  real alpha0;
  real alpha1;
  real alpha2;
  real<lower=0> beta0;
  vector[D] beta;
  real<lower=0> sigmad;
  vector[T] z_d;                       // mosquito death rate series
}
transformed parameters {
  vector[6] p0;
  vector[21] y0;
  vector[T] dv;
  vector[T + T_pred] mu_log_dv;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;

  // initial conditions
  
  p0 = softmax_id(p0_raw);

  y0[1] = p0[1];
  y0[2] = 0;
  y0[3] = 0;
  y0[4] = p0[2];
  y0[5] = p0[3];
  y0[6] = 0;
  y0[7] = 0;
  y0[8] = p0[4];
  y0[9] = 0;
  y0[10] = 0;
  y0[11] = p0[5];
  y0[12] = 0;
  y0[13] = 0;
  y0[14] = p0[6];
  y0[15] = exp(logNv0);
  y0[16] = 0.0;
  y0[17] = 0.0;
  y0[18] = 0.0;
  y0[19] = 0.0;
  y0[20] = 0.0;
  y0[21] = 0.0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // mosquito demographic parameters
  mu_log_dv = covars * beta;
  dv = 1.47 * beta0 * exp(head(mu_log_dv, T) + sigmad * z_d);
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  vector[21] state[T * 7 + 1];
  int idx = 1;
  
  // Priors
  
  // Initial conditions
  p0_raw[1] ~ normal(0.8, 0.2); // fully susceptible
  p0_raw[2] ~ normal(-11, 0.6); // cross-immune to 2
  p0_raw[3] ~ normal(0.3, 0.2); // susceptible to 2, immune to 1
  p0_raw[4] ~ normal(-11, 0.6); // exposed to 2
  p0_raw[5] ~ normal(0.3, 0.2); // susceptible to 1, immune to 2

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
  
  // Mosquito demographic series
  sigmad ~ normal(0, 0.1);
  
  alpha0 ~ normal(-1, 1);
  alpha1 ~ normal(0, 2);
  alpha2 ~ normal(0, 2);
  
  beta0 ~ gamma(100, 100);
  beta ~ normal(0, 2);

  z_d ~ normal(0, 1);
  
  // Process model
  
  state[1] = y0;
  
  for (t in 1:T){
    for(j in 1:7){
      
      state[idx + 1] = state[idx] + 1.0 / 7.0 * derivs(t, state[idx], 
                                                       alpha0,
                                                       alpha1, 
                                                       alpha2, 
                                                       rov[t], 
                                                       lambda, 
                                                       ro, 
                                                       gamma, 
                                                       dv[t], 
                                                       delta, 
                                                       phi_q * tau[t]);
        
      idx = idx + 1;
    }
    
    q_hat[t] = state[idx, 21] - state[idx - 7, 21];
    y_hat[t] = state[idx, 20] - state[idx - 7, 20];
  }
  
  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {

  vector[T + T_pred] y_hat;
  vector[T + T_pred] q_hat;
  vector[21] state;
  #vector[21] system[T];
  vector[T_pred] d_pred;
  
  state = y0;

  for (t in 1:T){
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(t, state, 
                                         alpha0,
                                         alpha1,
                                         alpha2,
                                         rov[t],
                                         lambda,
                                         ro,
                                         gamma,
                                         dv[t],
                                         delta,
                                         phi_q * tau[t]);
    }

    q_hat[t] = neg_binomial_2_rng(state[21] * pop, eta_q);
    y_hat[t] = neg_binomial_2_rng(phi_y * pop * state[20], eta_y);
    
    #system[t] = state;

    state[20] = 0;
    state[21] = 0;
    
  }
  for (k in 1:T_pred){
    
    d_pred[k] = 1.47 * beta0 * exp(mu_log_dv[T + k] + sigmad * normal_rng(0, 1));
    
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(T + k, state, 
                                         alpha0,
                                         alpha1,
                                         alpha2,
                                         rov[T + k],
                                         lambda,
                                         ro,
                                         gamma,
                                         d_pred[k],
                                         delta,
                                         phi_q * tau[T]);
    }

    q_hat[T + k] = neg_binomial_2_rng(state[21] * pop, eta_q);
    y_hat[T + k] = neg_binomial_2_rng(phi_y * pop * state[20], eta_y);

    state[20] = 0;
    state[21] = 0;
    
  }
}
