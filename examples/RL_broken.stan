// matstanlib/examples/RL_broken.stan
// 
// this Stan model specification implements the original flawed version of 
// a hierarchical Bayesian reinforcement learning model of performance on 
// an N-armed bandit task.  
// 
// (for the reparameterized version, see matstanlib/examples/RL_fixed.stan)
// 
// this code is part of the matstanlib library, which is available from 
// github.com/baribault/matstanlib.
// 
// (c) beth baribault 2021 --- 

data { 
  int<lower=1> S;                      //number of subjects
  int<lower=1> T;                      //number of trials per subject
  int<lower=1> A;                      //number of bandit arms
  int<lower=1> N;                      //number of data points
  int<lower=1,upper=S> Subject[N];     //subject number
  int<lower=1,upper=T> Trial[N];       //trial number
  int<lower=1,upper=A> Action[N];      //action selected
  int<lower=0,upper=1> Reward[N];      //reinforcement given
  vector[A] Q0;                        //inital value of actions
}

transformed data { 
  //fixed parameters
  real epsilon = 0.005;
}

parameters { 
  //group-level parameters
  real<lower=0> mu_beta;
  real<lower=0> sigma_beta;
  real<lower=0,upper=1> mu_alpha;
  real<lower=0,upper=1> sigma_alpha;
  real<lower=0,upper=1> mu_phi;
  real<lower=0,upper=1> sigma_phi;

  //subject-level parameters
  real<lower=0> beta[S];               //inverse temperature
  real<lower=0,upper=1> alpha[S];      //learning rate
  real<lower=0,upper=1> phi[S];        //decay rate
}

model { 
  //local variables
  vector[A] Q;
  vector[A] pi;
  real delta;

  //group-level priors
  mu_beta ~ normal(10,5);
  sigma_beta ~ normal(0,5);
  mu_alpha ~ uniform(0,1);
  sigma_alpha ~ normal(0,0.5);
  mu_phi ~ uniform(0,1);
  sigma_phi ~ normal(0,0.5);

  //individual-level priors
  for (s in 1:S) { 
    beta[s] ~ normal(mu_beta,sigma_beta);
    alpha[s] ~ normal(mu_alpha,sigma_alpha);
    phi[s] ~ normal(mu_phi,sigma_phi);
  }
  
  //likelihood
  for (n in 1:N) { 
    //initialize value
    if (Trial[n]==1) { 
      Q = Q0;
    } else {
      //... or forget values from previous trial
      Q = (1-phi[Subject[n]])*Q + phi[Subject[n]]*Q0;
    }
    
    //choice behavior
    pi = (1-epsilon)*softmax(beta[Subject[n]]*Q) + epsilon/A;
    Action[n] ~ categorical(pi);
    
    //reward
    delta = Reward[n] - Q[Action[n]];
    Q[Action[n]] = Q[Action[n]] + alpha[Subject[n]]*delta;
  }
}