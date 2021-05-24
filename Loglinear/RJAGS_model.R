model <- "model{
	#n.sample = total number of samples
	#n.ct = total number of cell types
	#N = total number of cells per sample

	for (s in 1:n.sample){
		for (c in 1:n.ct){	
			
			#sample rate is offset by total number of cells for the sample
			#X = matrix of cell counts with sample s as rows and cell type c as columns
			X[s,c] ~ dpois(N[s]*mu[s,c])
			
			#mu is the rates of cell types
			#log mu dependent on covariates
			log(mu[s,c]) <- b0[c] + b1[age.group[s],c] + b2[sex[s],c] + b3[batch[s],c] 	
			
			#probability of the cell types for each sample
			fit[s,c] <- mu[s,c]*exp(-mu[s,c])	


		}
	}
	
	for (c in 1:n.ct){
		#four groups: Younger Males, Younger Females, Older Males, Older Females in the same batch
		#log rate
		log(mu.YM[c]) <- b0[c] 
		log(mu.YF[c]) <- b0[c] + b2[2,c]
		log(mu.OM[c]) <- b0[c] + b1[2,c]
		log(mu.OF[c]) <- b0[c] + b1[2,c] + b2[2,c]

		#probability of the cell types for each sample
		fit.YM[c] <- mu.YM[c]*exp(-mu.YM[c])
		fit.YF[c] <- mu.YF[c]*exp(-mu.YF[c])
		fit.OM[c] <- mu.OM[c]*exp(-mu.OM[c])	
		fit.OF[c] <- mu.OF[c]*exp(-mu.OF[c])
	}

	#priors
	
  for (c in 1:n.ct){b0[c] ~ dnorm(0, 0.00001)}

  for (c in 1:n.ct){
      b1[1,c] <- 0
      b2[1,c] <- 0
      b1[2, c] ~ dnorm(0, 0.00001)
      b2[2, c] ~ dnorm(0, 0.00001)

    b3[1,c] <- 0
    for(i in 2:4){ #loop on # of batch
      b3[i, c] ~ dnorm(0, 0.00001)
    }  
  }
  


}
"
