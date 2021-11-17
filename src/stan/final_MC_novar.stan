
data {
  int<lower=1> T;  // total number of timepoints
  int<lower=1> I;  // total number of regions (utla)
  int<lower=1> J;  // total number of NUTS1 regions
  int<lower=0> nuts[I]; // individual ID
  int<lower=0> Y_delta[T, I];  // response variable = number of delta variant cases
  int<lower=0> Y_total[T, I];  // offset variable = total cases sampled
  // int<lower=0> Y_delta_test[T, I];  // response variable = number of delta variant cases
  // int<lower=0> Y_total_test[T, I];  // offset variable = total cases sampled
  int<lower=0,upper=1> include_lik[T,I]; //include in likelihood or not 
  int<lower=1> T1;  // number of extra weeks
  int<lower=0> Y_delta_extra[T1, I];  // extra = number of delta variant cases
  int<lower=0> Y_total_extra[T1, I];  // extra = offset variable = total cases sampled
  int<lower=0,upper=1> include_lik_extra[T1,I]; //include in likelihood or not 
}

parameters {
  real phi_1[1, I];
  real rho_ov_reg[T,J];   //  NUTS1 and time specific transmission advantage
  real dev_rho_raw[T,I];
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real a[J];
  real<lower=0> b[J];
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

  for(i in 1:I){
    for(t in 1:T){
        rho[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i];
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
  a ~ normal(-3,1);
  b ~ normal(1.5,1);
}

generated quantities {
  real loglikelihood_old[T,I];
  real dev_rho_extra[T1,I]; //extra dev_rho
  real rho_ov_reg_extra[T1,J];   //  extra = NUTS1 and time specific transmission advantage
  real rho_extra[T1,I];
  real phi_extra[T1, I];
  real loglikelihood_extra[T1,I];
  int Y_delta_hat_post_extra[T1,I];
  
  for(i in 1:I){
    for(t in 1:T){
      if(include_lik[t,i] == 1)
      loglikelihood_old[t,i] = binomial_logit_lpmf(Y_delta[t,i]|Y_total[t,i], phi[t,i]);
      else
      loglikelihood_old[t,i] = -99;
    }
  }
  
  // Get rho_ov_reg for the xtra time points and regions/NUTS1
  // First extra time point depends on the last time of data used in estimation
  for(j in 1:J){
    for(t1 in 1:T1){
      if(t1 == 1){
        rho_ov_reg_extra[t1,j] = normal_rng(rho_ov_reg[T,j], sigma2);
      }
      else{
        rho_ov_reg_extra[t1,j] = normal_rng(rho_ov_reg_extra[t1-1,j], sigma2);
      }
    }
  }
  
  // Get dev_rhos for all extra time points and all states/UTLAs
  for(i in 1:I){
    for(t1 in 1:T1){
      dev_rho_extra[t1,i] = normal_rng(0, sigma1);
      if(t1 == 1){
        rho_extra[t1,i] = rho_ov_reg_extra[t1,nuts[i]] + dev_rho_extra[t1,i];
        phi_extra[t1,i] = phi[T,i] + rho[T,i];
      }
      else{
        rho_extra[t1,i] = rho_ov_reg_extra[t1,nuts[i]] + dev_rho_extra[t1,i];
        phi_extra[t1,i] = phi_extra[t1-1,i] + rho_extra[t1-1,i];
      }
    }
  }

  for(i in 1:I){
    for(t1 in 1:T1){
      if(include_lik_extra[t1,i] == 1)
      loglikelihood_extra[t1,i] = binomial_logit_lpmf(Y_delta_extra[t1,i]|Y_total_extra[t1,i], phi_extra[t1,i]);
      else
      loglikelihood_extra[t1,i] = -99;
    }
  }

// Posterior
  for(i in 1:I){
    for(t1 in 1:T1){
      Y_delta_hat_post_extra[t1,i] = binomial_rng(Y_total_extra[t1,i], inv_logit(phi_extra[t1,i]));
    }
  }
}
