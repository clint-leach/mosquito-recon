functions {
  vector derivs(int t,
               vector y,
               real bv,
               real rov,
               real lambda,
               real ro,
               real gamma,
               real dv,
               real delta,
               real xi,
               real cap) {
    
    /**
    * documentation block
    */
    
    vector[9] dydt;
    
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
    Sv = y[7] - y[6] - y[5];
    R = 1 - y[1] - y[2] - y[3];

    loss = dv + cap;
    
    // Compute transition rates
    // Mosquito to human foi
    foi_vh = lambda * y[6] * y[1];

    // Human to mosquito foi
    foi_hv = lambda * y[3] * Sv;

    // Infectious humans
    infectious = ro * y[2];
    
    // Infectious mosq
    infect_mosq = rov * y[5];
    
    // Compute derivatives
    /*S*/  dydt[1] = b - b * y[1] - foi_vh + delta * R;
    
    /*E*/ dydt[2] = foi_vh - infectious - b * y[2];
    /*I*/ dydt[3] = infectious - (gamma + b) * y[3];
    
    /*VA*/ dydt[4] = xi * y[7] - bv * y[4];
    /*VE*/ dydt[5] = foi_hv - infect_mosq - loss * y[5];
    /*VI*/ dydt[6] = infect_mosq - loss * y[6];
    /*VN*/ dydt[7] = bv * y[4] - loss * y[7];
    
    /*VC*/ dydt[8] = cap * y[7];
    /*cases*/ dydt[9] = infectious;
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=0> T_pred;
  int<lower=0> D;
  int y[T];
  int q[T];
  vector[T] tau;
  real rov[T + T_pred];
  matrix[T + T_pred, D] covars;
  int week[T + T_pred];
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
  real<lower=0> beta0;
  vector[D] beta;
  vector[D] alpha;
  real<lower=0> sigmad;
  vector[T] eps_d;
}
transformed parameters {
  vector[9] y0;
  vector[T] dv;
  vector[T] bv;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;
  real<lower=0> xi;

  // initial conditions
  y0[1] = S0 * (pop - E0 - I0) / pop;
  y0[2] = E0 / pop;
  y0[3] = I0 / pop;
  y0[4] = exp(logNv);
  y0[5] = 0.0;
  y0[6] = 0.0;
  y0[7] = exp(logNv);
  y0[8] = 0.0;
  y0[9] = 0.0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // mosquito demographic parameters
  xi = 1.47 * beta0;
  bv = 1.47 * beta0 * exp(covars[1:T] * alpha);
  dv = 1.47 * beta0 * exp(covars[1:T] * beta + sigmad * eps_d);
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  vector[9] state[T * 7 + 1];
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
  beta0 ~ gamma(100, 100);
  beta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  
  sigmad ~ normal(0, 0.1);
  
  eps_d ~ normal(0, 1);

  // Process model
  
  state[1] = y0;

  for (t in 1:T){
    for(j in 1:7){
      
      state[idx + 1] = state[idx] + 1.0 / 7.0 * derivs(t, 
                                                       state[idx], 
                                                       bv[t], 
                                                       rov[t], 
                                                       lambda, 
                                                       ro, 
                                                       gamma, 
                                                       dv[t], 
                                                       delta, 
                                                       xi,
                                                       phi_q * tau[t]);
        
      idx = idx + 1;
    }
    
    q_hat[t] = state[idx, 8] - state[idx - 7, 8];
    y_hat[t] = state[idx, 9] - state[idx - 7, 9];
  }
  
  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {

  vector[T + T_pred] y_hat;
  vector[T + T_pred] q_hat;
  vector[T_pred] d_pred;
  vector[T_pred] b_pred;
  vector[9] state;
  
  state = y0;

  // Estimated trajectories
  for (t in 1:T){
    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(t, 
                                         state, 
                                         bv[t],
                                         rov[t], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         dv[t], 
                                         delta, 
                                         xi,
                                         phi_q * tau[t]);
    }

    q_hat[t] = neg_binomial_2_rng(state[8] * pop, eta_q);
    y_hat[t] = neg_binomial_2_rng(phi_y * pop * state[9], eta_y);

    state[8] = 0;
    state[9] = 0;
    
  }
  for (k in 1:T_pred){
    
    d_pred[k] = 1.47 * beta0 * exp(covars[T + k] * beta + normal_rng(0, sigmad));
    b_pred[k] = 1.47 * beta0 * exp(covars[T + k] * alpha);

    for(j in 1:7){
      
      state = state + 1.0 / 7.0 * derivs(T + k, 
                                         state, 
                                         b_pred[k],
                                         rov[T + k], 
                                         lambda, 
                                         ro, 
                                         gamma, 
                                         d_pred[k], 
                                         delta, 
                                         xi,
                                         phi_q * tau[T]);
    }

    q_hat[T + k] = neg_binomial_2_rng(state[8] * pop, eta_q);
    y_hat[T + k] = neg_binomial_2_rng(phi_y * pop * state[9], eta_y);

    state[8] = 0;
    state[9] = 0;
    
  }
}
