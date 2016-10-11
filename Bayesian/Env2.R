sink("Bayesian/Env2.jags")
cat("
    model{
    
    #Liklihood
    
    for (cell in 1:cells){

      #Behavior
      beh[cell] ~ dbern(rho[cell] * X[cell]+0.000001)
  
      #Occurrence
      X[cell] ~ dbern(phi[cell])
      
      #Occ function
      logit(phi[cell]) = alpha + beta * env[cell]  
  
      #Behavior function
      logit(rho[cell]) = alpha2 + beta2 * env2[cell]  

    }
    
    #Priors
    alpha ~ dnorm(0,0.386)
    beta ~ dnorm(0,0.386)
    
    alpha2 ~ dnorm(0,0.386)
    beta2 ~ dnorm(0,0.386)
    
    }"
    ,fill=TRUE)
sink()