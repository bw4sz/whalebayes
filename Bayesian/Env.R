sink("Bayesian/Env.jags")
cat("
    model{
    
    #Liklihood
    
    for (cell in 1:cells){
    
    #Conditional Behavior
    v[cell] ~ dbern(z[cell])
    
    z[cell] <- rho * X[cell]
    #Occurrence
    X[cell] ~ dbern(phi[cell])
    
    #Occ function
    logit(phi[cell]) = alpha + beta * env[cell]  

    }
    
    #marginal probability
    m=mean(z[])

    #Priors
    alpha ~ dnorm(0,0.386)
    beta ~ dnorm(0,0.386)
    
    rho ~ dbeta(1,1)
    
    }"
    ,fill=TRUE)
sink()