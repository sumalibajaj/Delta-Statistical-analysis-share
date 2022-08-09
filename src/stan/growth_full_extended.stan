
data {
  int<lower=1> T;  // total number of timepoints
  int<lower=1> I;  // total number of regions (utla)
  int<lower=1> J;  // total number of NUTS1 regions
  int<lower=0> nuts[I]; // individual ID
  int Y_delta[T, I];  // response variable = number of delta variant cases
  int Y_total[T, I];  // offset variable = total cases sampled
  int<lower=0,upper=1> include_lik[T,I]; //include in likelihood or not 
  int<lower=0> ncovs; // number of covariates
  int<lower=0> ncovs2; // number of covariates which remain the same as observed in counterfactuals
  matrix[T*I, ncovs] X;
  matrix[T*I, ncovs] X_90;
  matrix[T*I, ncovs2] X_time;
}

parameters {
  real rho_ov_reg[T,J];   //  NUTS1 and time specific transmission advantage
  real dev_rho_raw[T,I]; // for non-centred parametrisation
  vector[ncovs] beta;
  real phi_1[1, I];
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real a[J];
  real<lower=0> b[J];
  vector[ncovs2] beta_time;
}

transformed parameters{
  real dev_rho[T,I];
  real phi[T, I];
  real rho[T,I];
  real p[T, I];

  for(i in 1:I){
    for(t in 1:T){
      dev_rho[t,i] = dev_rho_raw[t,i]*sigma1;
    }
  }
  
  {int c = 1;
  for(i in 1:I){
    for(t in 1:T){
        rho[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i] + X[c]*beta + X_time[c]*beta_time;
      c = c+1;
    }
  }
  }
  
  for(i in 1:I){
    for(t in 1:T){
      if(t==1){
        phi[t,i] = phi_1[t,i];
      }
      else{
        phi[t,i] = phi[t-1,i] + rho[t-1,i];
      }
    }
  }
  
  p = inv_logit(phi);
  
}

model {
  for(t in 1:T){
    for(i in 1:I){
      if(include_lik[t,i] == 1)
      Y_delta[t,i] ~ binomial_logit(Y_total[t,i], phi[t,i]);
    }
  }
  
  for(j in 1:J){
    for(t in 2:T){
      rho_ov_reg[t,j] ~ normal(rho_ov_reg[t-1,j], sigma2);
    }
  }

  for(i in 1:I){
    for(t in 1:T){
      dev_rho_raw[t,i] ~ normal(0, 1);
    }
  }

  for(i in 1:I)
    phi_1[1,i] ~ normal(a[nuts[i]],b[nuts[i]]);
    
  for(j in 1:J)
    rho_ov_reg[1,j] ~ normal(0,1);

  sigma1 ~ normal(0,5);
  sigma2 ~ normal(0,5);
  beta ~ normal(0,5);
  beta_time ~ normal(0,5);
  
  a ~ normal(-3,1);
  b ~ normal(1.5,1);
}

generated quantities {
  int Y_delta_hat_post[T,I];
  int Y_delta_hat_prior[T, I];
  real p_wovac[T,I];
  real p_90vac[T,I];
  {
  real a_prior[J];
  real b_prior[J];
  real phi_1_prior[T,I];
  real rho_wovac[T,I];
  real phi_wovac[T,I];
  real rho_90vac[T,I];
  real phi_90vac[T,I];
  real rho_ov_reg_prior[T,J];
  real sigma1_prior;
  real sigma2_prior;
  vector[ncovs] beta_prior;
  real dev_rho_raw_prior[T,I];
  real rho_prior[T,I];
  real dev_rho_prior[T,I];
  real phi_prior[T,I];

  
  // State level p counterfactuals
  {int c = 1;
    for(i in 1:I){
    for(t in 1:T){
        rho_wovac[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i] + X_time[c]*beta_time;
        rho_90vac[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i] + X_90[c]*beta + X_time[c]*beta_time;
        c = c+1;
    }
    }
    
    for(i in 1:I){
    for(t in 1:T){
      if(t==1){
        phi_wovac[t,i] = phi_1[1,i];
        phi_90vac[t,i] = phi_1[1,i];
      }
      else{
        phi_wovac[t,i] = phi_wovac[t-1,i] + rho_wovac[t-1,i];
        phi_90vac[t,i] = phi_90vac[t-1,i] + rho_90vac[t-1,i];
      }
    }
  }
  p_wovac = inv_logit(phi_wovac);
  p_90vac = inv_logit(phi_90vac);
  }

// Posterior predictive checks
  for(i in 1:I){
    for(t in 1:T){
      if(include_lik[t,i] == 1)
      Y_delta_hat_post[t,i] = binomial_rng(Y_total[t,i], inv_logit(phi[t,i]));
    }
  }

// Prior predictive checks
  for(j in 1:J)
    rho_ov_reg_prior[1,j] = normal_rng(0,1);

  sigma1_prior = normal_rng(0,5);
  while(sigma1_prior<0){
      sigma1_prior = normal_rng(0,5);
  }
  
  sigma2_prior = normal_rng(0,5);
  while(sigma2_prior<0){
      sigma2_prior = normal_rng(0,5);
  }
  
  for(n in 1:ncovs){
    beta_prior[n] = normal_rng(0,5);
    while(beta_prior[n]<0){
      beta_prior[n] = normal_rng(0,5);
  }
  }
  
  for(j in 1:J){
    a_prior[j] = normal_rng(-3,1);
    b_prior[j] = normal_rng(1.5,1);
    while(b_prior[j]<0){
      b_prior[j] = normal_rng(1.5,1);
    }
  }
  for(i in 1:I){
    phi_1_prior[1,i] = normal_rng(a_prior[nuts[i]],b_prior[nuts[i]]);
  }

    

  for(i in 1:I){
    for(t in 1:T){
      dev_rho_raw_prior[t,i] = normal_rng(0, 1);
    }
  }
  
  
  for(j in 1:J){
    for(t in 2:T){
      rho_ov_reg_prior[t,j] = normal_rng(rho_ov_reg_prior[t-1,j], sigma2_prior);
    }
  }
  
  
  for(i in 1:I){
    for(t in 1:T){
      dev_rho_prior[t,i] = dev_rho_raw_prior[t,i]*sigma1_prior;
    }
  }
  
  
  {int c = 1;
  for(i in 1:I){
    for(t in 1:T){
        rho_prior[t,i] = rho_ov_reg_prior[t,nuts[i]] + dev_rho_prior[t,i] + X[c]*beta_prior;
      c = c+1;
    }
  }
  }
  
  for(i in 1:I){
    for(t in 1:T){
      if(t==1){
        phi_prior[t,i] = phi_1_prior[t,i];
      }
      else{
        phi_prior[t,i] = phi_prior[t-1,i] + rho_prior[t-1,i];
      }
    }
  }
  
  
  for(t in 1:T){
    for(i in 1:I){
      if(include_lik[t,i] == 1)
      Y_delta_hat_prior[t,i] = binomial_rng(Y_total[t,i], inv_logit(phi_prior[t,i]));
    }
  }


}
}
