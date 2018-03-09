functions {
  vector derivs(int t,
               vector y,
               real rov,
               real lambda,
               real ro,
               real gamma,
               real eps_dv,
               real delta,
               real omega,
               real damp, 
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

    dv = exp(y[9] + 0.39);
    
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
    
    /*VE*/ dydt[4] = foi_hv - infect_mosq - (control[1] * dv) * y[4];
    /*VI*/ dydt[5] = infect_mosq - (control[1] * dv) * y[5];
    /*VN*/ dydt[6] = control[2] * (y[11] + dv) * y[6] - (control[1] * dv) * y[6];
    
    /*VC*/ dydt[7] = y[6];
    /*cases*/ dydt[8] = infectious;
    
    // Parameter processes
    /*logdv*/ dydt[9] = y[10];
    /*dlogdv*/ dydt[10] = -2 * damp * omega * y[10] - omega ^ 2 * y[9] + eps_dv;
    /*rv*/ dydt[11] = y[12];
    /*drv*/ dydt[12] = -2 * damp * omega * y[12] - omega ^ 2 * y[11];
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=1> steps;
  int y[T];
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
  real<lower=0> eta_inv_y;             // overdispersion of case reports
  real<lower=0> ro_c;                  // human latenet period
  real<lower=0> gamma_c;               // human infectious period
  real<lower=0> delta_c;               // cross-immune period
  real logNv;                         // initial mosquito population size
  real dv0;
  real rv0;
  real<lower=0> sigmadv;
  vector[T] eps_dv;
}
transformed parameters {
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  vector[T] y_hat;
  vector[12] state[T + 1];
  
  // initial conditions
  state[1, 1] = S0 * (pop - E0 - I0) / pop;
  state[1, 2] = E0 / pop;
  state[1, 3] = I0 / pop;
  state[1, 4] = 0.0;
  state[1, 5] = 0.0;
  state[1, 6] = exp(logNv);
  state[1, 7] = 0.0;
  state[1, 8] = 0.0;
  state[1, 9] = dv0;
  state[1, 10] = 0;
  state[1, 11] = rv0;
  state[1, 12] = 0;
  
  // measurement parameters
  eta_y = 1 / eta_inv_y;
  
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
    
    for(j in 1:steps){
      
      state[t] = state[t] + 1.0 / steps * derivs(t, 
                                                 state[t], 
                                                 rov[t - 1], 
                                                 lambda, 
                                                 ro, 
                                                 gamma, 
                                                 sigmadv * eps_dv[t - 1], 
                                                 delta, 
                                                 omega, 
                                                 damp, 
                                                 control[t - 1]);
    }
                                                      
    y_hat[t - 1] = state[t, 8] - state[t - 1, 8];

  }
}
model {

  // Priors
  
  // Initial conditions
  S0 ~ beta(4, 6);
  E0 ~ gamma(10, 0.1);
  I0 ~ gamma(6, 0.1);

  // Measurement models
  eta_inv_y ~ normal(0, 5);

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
  eps_dv ~ normal(0, 1);
  
  // Variance parameters
  sigmadv ~ normal(0, 0.2);

  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
}
generated quantities {
  vector[T] y_meas;

  // Estimated trajectories
  for (t in 1:T){
    
    y_meas[t] = neg_binomial_2_rng(phi_y * y_hat[t] * pop, eta_y);
    
  }
}