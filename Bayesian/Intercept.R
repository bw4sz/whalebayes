sink("Bayesian/Intercept.jags")
cat("
model{

    #Liklihood

    for (cell in 1:cells){
      X[cell] ~ dbern(z)
    }

    #conditional
    z=rho * phi

    #Priors
    rho ~ dbeta(1,1)
    phi ~ dbeta(1,1)

    }"
    ,fill=TRUE)
sink()