sink("Bayesian/Env.jags")
cat("
    model{
    
    #Liklihood
    
    for (cell in 1:cells){
      X[cell] ~ dbern(z[cell])
      z[cell] = rho * phi[cell]
  
      #Env function
      logit(phi[cell]) = alpha + beta * env[cell]  
      }
    

    #Priors
    rho ~ dbeta(1,1)
    alpha ~ dnorm(0,0.386)
    beta ~ dnorm(0,0.386)

    
    }"
    ,fill=TRUE)
sink()