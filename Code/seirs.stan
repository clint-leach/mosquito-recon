/**
Implements model of dengue transmission defined by SEIRS structure in the human
population and SEI structure in the mosquito population, with an exponentially-distributed
extrinsic incubation period.

Estimates time series of mosquito mortality and population growth rates from time
series of weekly moquito trap counts and reports of "dengue-like illness".
*/
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
    Computes the RHS of the SERS-SEI dengue transmission model.
    Inputs:
      t: integer giving the week index
      y: vector of length 12 giving the values of the state variables in week t
      eps_rv: real giving the week t forcing term for the rv oscillator
      rov: real giving the inverse of the extrinsic incubation period for week t
      lambda: real giving the transmission rate
      ro: real giving the inverse of the human latent period
      gamma: real giving the rate of loss of infectiousness in humans
      eps_dv: real giving the week t forcing term for the log(dv) oscillator
      delta: real giving the inverse of the period of cross-immunity
      cap: real giving the week t mosquito capture rate
      omega: real giving the period of the harmonic oscillators
      damp: real giving the damping coefficient of the harmonic oscillators
      control: vector of length 3, giving effects of control efforts in week t
               [1] the multiplicative change in mosquito mortality, 
               [2] the multiplicative change in mosquito birth rate,
               [3] the multiplicative change in adult mosquito abundance
               
    Outputs:
      dydt: a vector of length 12 giving the RHS of the ODE at time t
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
    
    /*VE*/ dydt[4] = foi_hv - infect_mosq - (control[1] * dv + cap) * y[4];
    /*VI*/ dydt[5] = infect_mosq - (control[1] * dv + cap) * y[5];
    /*VN*/ dydt[6] = control[2] * (y[11] + dv) * y[6] - (control[1] * dv + cap) * y[6];
    
    /*VC*/ dydt[7] = cap * y[6];
    /*cases*/ dydt[8] = infectious;
    
    // Parameter processes
    /*logdv*/ dydt[9] = y[10];
    /*dlogdv*/ dydt[10] = -2 * damp * omega * y[10] - omega ^ 2 * y[9] + eps_dv;
    /*rv*/ dydt[11] = y[12];
    /*drv*/ dydt[12] = -2 * damp * omega * y[12] - omega ^ 2 * y[11] + eps_rv;
    
    return dydt;
  }
}
data {
  int<lower=1> T;          // number of time steps
  int<lower=1> steps;      // number of Euler steps per week
  int y[T];                // length T vector of case reports data
  int q[T];                // length T vector of mosquito trap counts
  vector[T] tau;           // length T vector of number of mosquito traps inspected
  vector[T] rov;           // length T vector of weekly average extrinsic incubation period
  vector[3] control[T];    // T by 3 array of mosquito control actions
  int pop;                 // integer giving the human population size
}
transformed data {
  real lambda = 4.87;                              // DENV transmission rate
  real phi_y = 1.0 / 12.0;                         // proportion of human cases reported
  real damp = 0.0;                                 // damping coefficient on harmonic oscillators
  real omega = (pi() / 26) / sqrt(1 - damp ^ 2);   // period of harmonic oscillators
}
parameters {
  real<lower=0,upper=1> S0;            // initial proportion of human population susceptible
  real<lower=0> E0;                    // initial number of exposed humans
  real<lower=0> I0;                    // initial number of infectious humans
  real<upper=0> log_phi_q;             // log per-trap capture rate
  real<lower=0> eta_inv_y;             // overdispersion of case reports
  real<lower=0> eta_inv_q;             // overdispersion of mosquito capture
  real<lower=0> ro_c;                  // centered human latenet period
  real<lower=0> gamma_c;               // centered human infectious period
  real<lower=0> delta_c;               // centered cross-immune period
  real logNv;                          // initial log mosquito population size
  real dv0;                            // initial log mosquito mortality rate
  real rv0;                            // initial mosquito population growth rate
  real<lower=0> sigmadv;               // std dev of mortality rate perturbation process
  real<lower=0> sigmarv;               // std dev of growth rate perturbation process
  vector[T] eps_dv;                    // perturbations to log(dv) oscillator
  vector[T] eps_rv;                    // perturbations to rv oscillator
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
  eta_q = 1 / eta_inv_q;
  phi_q = exp(log_phi_q);
  
  // rate parameters
  ro = 1 / (0.87 * ro_c);
  gamma = 3.5 * gamma_c;
  delta = 1 / (97 * delta_c);
  
  // Process model
  for (t in 2:(T + 1)){
    
    state[t] = state[t - 1];
    
    // Implementing control actions
    state[t, 4] = state[t, 4] * control[t - 1, 3];
    state[t, 5] = state[t, 5] * control[t - 1, 3];
    state[t, 6] = state[t, 6] * control[t - 1, 3];
    
    // Implementing Euler integration steps
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
                                                      
    q_hat[t - 1] = state[t, 7] - state[t - 1, 7];
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
