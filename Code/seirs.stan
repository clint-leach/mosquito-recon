functions {
  real[] derivs(real t,
                real[] y,
                real[] theta,
                real[] x_r,
                int[] x_i) {
    
    /**
    * documentation block
    */
    
    real dydt[8];
    
    real b;
    real lambda;
    real omega;
    real Sv;
    real R;
    
    real foi_vh;
    real foi_hv;
    real infectious;
    real infect_mosq;
    
    real ro = theta[1];
    real delta = theta[2];
    real gamma = theta[3];
    real cap = theta[4] * x_i[1];

    real rov = x_r[1];
    real ac = x_r[2];
    real lc = x_r[3];
    
    real bw = 8;
    
    vector[60] times = to_vector(tail(x_r, 60));
    vector[60] z_dv = to_vector(theta[5:64]);
    vector[60] z_rv =  to_vector(tail(theta, 60));
    
    vector[60] weights = 1 / (bw * sqrt(2 * pi())) * exp(-0.5 * square((t - times) / bw));
    real dv = sum(weights .* z_dv);
    real rv = sum(weights .* z_rv);
    
    // Assigning data
    b = 7.0 / (76 * 365);
    lambda = 4.87;

    // Computing mosquito population size
    Sv = y[6] - y[5] - y[4];
    R = 1 - y[1] - y[2] - y[3];

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
    
    /*VE*/ dydt[4] = foi_hv - infect_mosq - (ac * dv + cap) * y[4];
    /*VI*/ dydt[5] = infect_mosq - (ac * dv + cap) * y[5];
    /*VN*/ dydt[6] = lc * (rv + dv) * y[6] - (ac * dv + cap) * y[6];
    
    /*VC*/ dydt[7] = cap * y[6];
    /*cases*/ dydt[8] = infectious;
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  int<lower=0> T_pred;
  real ts[T + 1, 1];
  real pc_times[60];
  int y[T];
  int q[T];
  int tau[T, 1];
  real rov[T, 1];
  real control[T, 2];
  int pop;
}
transformed data {
  real data_real[T, 63];
  real phi_y = 1.0 / 12.0;
  
  data_real[, 1] = rov[, 1];
  data_real[, 2] = control[, 1];
  data_real[, 3] = control[, 2];
  for(i in 4:63){
    data_real[, i] = rep_array(pc_times[i-3], T);
  }
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
  real<lower=0> z_dv[60];
  real z_rv[60];
}
transformed parameters {
  real y0[8];
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
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  real state[T + 1, 8];
  real theta[124];

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
  z_dv ~ gamma(1, 1);
  z_rv ~ normal(0, 0.1);

  // Process model

  theta[1] = ro;
  theta[2] = delta;
  theta[3] = gamma;
  theta[4] = phi_q;
  theta[5:64] = z_dv;
  theta[65:124] = z_rv;

  state[1] = y0;

  for (i in 1:T){

    state[i + 1] = integrate_ode_rk45(derivs,
                                      state[i],
                                      ts[i, 1],
                                      ts[i + 1],
                                      theta,
                                      data_real[i],
                                      tau[i],
                                      1e-6,
                                      1e-6,
                                      200)[1];

    q_hat[i] = state[i + 1, 7] - state[i, 7];
    y_hat[i] = state[i + 1, 8] - state[i, 8];

  }

  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {

  vector[T] y_hat;
  vector[T] q_hat;
  real state[T + 1, 8];
  real theta[124];
  
// Process model
  
  theta[1] = ro;
  theta[2] = delta;
  theta[3] = gamma;
  theta[4] = phi_q;
  theta[5:64] = z_dv;
  theta[65:124] = z_rv;

  state[1] = y0;

  for (i in 1:T){

    state[i + 1] = integrate_ode_rk45(derivs, 
                                      state[i], 
                                      ts[i, 1], 
                                      ts[i + 1], 
                                      theta, 
                                      data_real[i], 
                                      tau[i], 
                                      1e-6, 
                                      1e-6, 
                                      200)[1];
    
    q_hat[i] = neg_binomial_2_rng(pop * (state[i + 1, 7] - state[i, 7]), eta_q);
    y_hat[i] = neg_binomial_2_rng(phi_y * pop * (state[i + 1, 8] - state[i, 8]), eta_y);

  }
}
