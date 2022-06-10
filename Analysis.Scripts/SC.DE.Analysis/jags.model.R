#RJAGS bayesian gaussian model for the analysis of cell type specific differential genes of aging

model <- "
	model  
{    
#the expression of each gene follows a log normal distribution
#can fit normal disribution to log tranformed gene expression data


for(c in 1:N.cell){

	y.gene[c] ~ dnorm(mu[c],tau)

	#age group reference: younger age group
	mu[c] <- b1 * age.group.middle[c] + b2 * age.group.old[c] + b3 * age.group.EL[c] + b4 * sex[c] + b5 * (batch[c] - mean(batch[])) + b6 * ethnic[c] + bs[subj[c]]           

}

#priors for subject random effect
for(j in 1:N.subj){
	bs[j] ~ dnorm(b0,tau.s)
}

#priors for fixed effects
b0 ~ dnorm(0,0.0001)
b1 ~ dnorm(0,0.0001)
b2 ~ dnorm(0,0.0001)
b3 ~ dnorm(0,0.0001)
b4 ~ dnorm(0,0.0001)
b5 ~ dnorm(0,0.0001)
b6 ~ dnorm(0,0.0001)

### variance components
tau ~ dgamma(0.0001,0.0001)
tau.s ~ dgamma(0.00001,0.00001)



}
"
