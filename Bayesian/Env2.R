sink("Bayesian/Env2.jags")
cat("
    model{
    
    #Liklihood
    
    for (cell in 1:cells){
    X[cell] ~ dbern(z[cell])
    z[cell] = rho[cell] * phi[cell]
    
    #Env function
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