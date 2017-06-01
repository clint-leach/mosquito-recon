functions {
  real[] derivs(real t,
                real[] y,
                real[] theta,
                real[] x_r,
                int[] x_i) {
    
    /**
    * documentation block
    */
    
    real dydt[21];
    
    real b;
    real bv;
    real lambda;
    real Sv;
    real R;
    
    real foi_vh_1;
    real foi_vh_2;
    
    real foi_hv_1;
    real foi_hv_2;
    
    real alpha0 = theta[1];
    real alpha1 = theta[2];
    real alpha2 = theta[3];
    real ro = theta[4];
    real delta = theta[5];
    real gamma = theta[6];
    real cap = theta[7] * x_i[1];
    real dv = theta[8];
    
    real rov = x_r[1];
    
    // Assigning data
    b = 7.0 / (76 * 365);
    lambda = 4.87;
    bv = exp(alpha0 + alpha1 * sin(2 * pi() * t / 52) + alpha2 * cos(2 * pi() * t / 52));

    // Computing mosquito population size
    Sv = y[15] - y[16] - y[17] - y[18] - y[19];
    
    // Compute forces of infection
    // Mosquito to human
    foi_vh_1 = lambda * y[17];
    foi_vh_2 = lambda * y[19];
    
    // Human to mosquito
    foi_hv_1 = lambda * (y[3] + y[13]);
    foi_hv_2 = lambda * (y[7] + y[9]);
    
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
  real ts[T + 1, 1];
  int<lower=1> D;
  int y[T];
  int q[T];
  int tau[T, 1];
  real rov[T + T_pred, 1];
  matrix[T + T_pred, D] covars;
  int pop;
}
transformed data {
  matrix[T, D] Q_ast;
  matrix[D, D] R_ast;
  matrix[D, D] R_ast_inv;
  real phi_y = 1.0 / 12.0;

  Q_ast = qr_Q(covars[1:T, ])[, 1:D] * sqrt(T - 1);
  R_ast = qr_R(covars[1:T, ])[1:D, ] / sqrt(T - 1);
  R_ast_inv = inverse(R_ast);
}
parameters {
  vector[13] p0_raw;                   // untransformed initial conditions
  real log_phi_q;                      // log per-trap capture rate
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
  vector[T] logdv;
}
transformed parameters {
  vector[14] p0;
  real y0[21];
  vector[T] dv;
  real<lower=0> ro;
  real<lower=0> gamma;
  real<lower=0> delta;
  real<lower=0> eta_y;
  real<lower=0> eta_q;
  real<lower=0> phi_q;

  // initial conditions

  p0 = softmax_id(p0_raw);

  for(i in 1:14){
    y0[i] = p0[i];
  }
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
  dv = 1.47 * beta0 * exp(logdv);
}
model {
  vector[T] y_hat;
  vector[T] q_hat;
  real state[T + 1, 21];
  real theta[8];
  
  // Priors
  
  // Initial conditions
  p0_raw[1] ~ normal(0.8, 0.2);
  p0_raw[2] ~ normal(-9, 0.6);
  p0_raw[3] ~ normal(-9, 0.6);
  p0_raw[4] ~ normal(-11, 0.6); // Cross-immune to 2
  p0_raw[5] ~ normal(0.3, 0.2); // susceptible to 2, immune to 1
  p0_raw[6] ~ normal(-11, 0.6);
  p0_raw[7] ~ normal(-11, 0.6);
  p0_raw[8] ~ normal(-11, 0.6);
  p0_raw[9] ~ normal(-11, 0.6);
  p0_raw[10] ~ normal(-2, 0.6); // cross-immune to 1
  p0_raw[11] ~ normal(0.3, 0.2); // susceptible to 1, immune to 2
  p0_raw[12] ~ normal(-9, 0.6);
  p0_raw[13] ~ normal(-9, 0.6);


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
  beta ~ normal(0, 5);

  logdv ~ student_t(10, Q_ast * beta, sigmad);
  
  // Process model
  state[1] = y0;
  
  theta[1] = alpha0;
  theta[2] = alpha1;
  theta[3] = alpha2;
  theta[4] = ro;
  theta[5] = delta;
  theta[6] = gamma;
  theta[7] = phi_q;
  
  for(i in 1:T){
  
    theta[8] = dv[i];
    state[i, 20] = 0;
    state[i, 21] = 0;

    state[i + 1] = integrate_ode_rk45(derivs, state[i], ts[i, 1], ts[i + 1], theta, rov[i], tau[i], 
                                      1e-3, 1e-3, 7)[1];
      
    q_hat[i] = state[i + 1, 21];
    y_hat[i] = state[i + 1, 20];
    
  }
  
  // Measurement models
  y ~ neg_binomial_2(phi_y * y_hat * pop, eta_y);
  q ~ neg_binomial_2(q_hat * pop, eta_q);
}
generated quantities {

  vector[T + T_pred] y_hat;
  vector[T + T_pred] q_hat;
  vector[T + T_pred] mu_logdv;
  vector[T_pred] d_pred;

  mu_logdv = covars * R_ast_inv * beta;

  {
    real state[T + 1, 21];
    real theta[8];
    
    state[1] = y0;

    theta[1] = alpha0;
    theta[2] = alpha1;
    theta[3] = alpha2;
    theta[4] = ro;
    theta[5] = delta;
    theta[6] = gamma;
    theta[7] = phi_q;
    
    for(i in 1:T){
      
      theta[8] = dv[i];

      state[i + 1] = integrate_ode_rk45(derivs, state[i], ts[i, 1], ts[i + 1], theta, rov[i], tau[i], 
                                        1e-3, 1e-3, 7)[1];
      
      q_hat[i] = neg_binomial_2_rng(pop * (state[i + 1, 21] - state[i, 21]), eta_q);
      y_hat[i] = neg_binomial_2_rng(phi_y * pop * (state[i + 1, 20]- state[i, 20]), eta_y);
    }
  }
}
