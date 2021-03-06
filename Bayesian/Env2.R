sink("Bayesian/Env2.jags")
cat("
    model{
    
    #Liklihood
    
    for (cell in 1:cells){
    
    #Conditional Behavior
    v[cell] ~ dbern(z[cell])
    
    z[cell] <- rho[cell] * X[cell]
    #Occurrence
    X[cell] ~ dbern(phi[cell])
    
    #Occ function
    logit(phi[cell]) = alpha + beta * env[cell]  
    
    #Behavior 
    logit(rho[cell]) = alpha2 + beta2 * env2[cell]  

    }
    
    #marginal probability
    m=mean(z[])
    
    #Priors
    alpha ~ dnorm(0,0.386)
    beta ~ dnorm(0,0.386)
    
    alpha2 ~ dnorm(0,0.386)
    beta2 ~ dnorm(0,0.386)    
    }"
    ,fill=TRUE)
sink()