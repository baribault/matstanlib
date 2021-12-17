// matstanlib/examples/RL_broken.stan
// 
// this Stan model specification implements the reparameterized version of 
// a hierarchical Bayesian reinforcement learning model of performance on 
// an N-armed bandit task.  
// 
// (for the original version, see matstanlib/examples/RL_broken.stan)
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
  real<lower=0> b1;
  real<lower=0> mu_beta;
  real<lower=0> a1;
  real<lower=0> a2;
  real<lower=0> p1;
  real<lower=0> p2;
  
  //subject-level parameters
  real<lower=0> beta[S];               //inverse temperature
  real<lower=0,upper=1> alpha[S];      //learning rate
  real<lower=0,upper=1> phi[S];        //decay rate
} 

transformed parameters { 
  //derived parameters
  real b2 = (1+b1)/mu_beta;
  real sigma_beta = sqrt((1+b1) / (b2^2));
  real mu_alpha = (a1+1)/(a1+1 + a2+1);
  real sigma_alpha = sqrt((a1+1)*(a2+1)/( ((a1+1 + a2+1)^2)*(a1+1 + a2+1 + 1) ));
  real mu_phi = (p1+1)/(p1+1 + p2+1);
  real sigma_phi = sqrt((p1+1)*(p2+1)/( ((p1+1 + p2+1)^2)*(p1+1 + p2+1 + 1) ));
}

model { 
  //local variables
  vector[A] Q;
  vector[A] pi;
  real delta;
  
  //group-level priors
  b1 ~ gamma(5,1);
  mu_beta ~ normal(7.5,2.5);
  a1 ~ gamma(2,1);
  a2 ~ gamma(5,1);
  p1 ~ gamma(1,1);
  p2 ~ gamma(5,1);
  
  //individual-level priors
  for (s in 1:S) { 
    beta[s] ~ gamma(1+b1,b2);
    alpha[s] ~ beta(1+a1,1+a2);
    phi[s] ~ beta(1+p1,1+p2);
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