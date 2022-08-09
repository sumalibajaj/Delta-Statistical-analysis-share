
data {
  int<lower=1> T;  // total number of timepoints
  int<lower=1> I;  // total number of regions (utla)
  int<lower=1> J;  // total number of NUTS1 regions
  int<lower=0> nuts[I]; // individual ID
  int Y_delta[T, I];  // response variable = number of delta variant cases
  int Y_total[T, I];  // offset variable = total cases sampled
  int<lower=0,upper=1> include_lik[T,I]; //include in likelihood or not 
  int<lower=0> ncovs; // number of covariates
  matrix[T*I, ncovs] X;
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
        rho[t,i] = rho_ov_reg[t,nuts[i]] + dev_rho[t,i] + X[c]*beta;
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
  
  a ~ normal(-3,1);
  b ~ normal(1.5,1);
}
