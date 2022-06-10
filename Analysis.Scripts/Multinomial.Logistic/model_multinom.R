# Description of the Bayesian model fitted in this file
# Notation
# sample i = 1:N.sample
# cell type j = 1:N.ct
# y[i,j] = matrix of counts with each sample row i having set row sum S = 20 on j = 1:N.ct cell types
# p[i,j] = probability of sample i being in category j
# b1[age.group[i],j] = age group effect parameter for each sample i in cell type j
# b2[sex[i],j] = sex effect parameter for each sample i in cell type j


multinom.model <- "
model
{
  # Likelihood
  for (i in 1:N.sample) { # loop through samples
    y[i,1:N.ct] ~ dmulti(p[i,1:N.ct], N.total[i])
    for(j in 1:N.ct) { # loop through cell types

      log(q[i,j]) <- b0[j] +
                      b1[age.group[i],j] + #age.group  
                      b2[sex[i], j] #sex
                      
      p[i,j] <- q[i,j]/sum(q[i,1:N.ct])                  
    }
  }
  
  
  # Priors
   
  for(j in 1:N.ct){
      b1[1,j] <- 0 
      for(k in 2:N.age){
        b1[k,j] ~ dnorm(0, 0.001)
      }
      b2[1,j] <- 0 
      b2[2,j] ~ dnorm(0, 0.001)
      b0[j]  ~ dnorm(0, 0.001)
    
    }

}
"

