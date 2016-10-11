sink("Bayesian/Intercept.jags")
cat("
model{

    #Liklihood

    for (cell in 1:cells){
      
      beh[cell] ~ dbern(rho*X[cell]+0.000001)

      X[cell] ~ dbern(phi)
    }

    #conditional
    #z=rho * phi

    rho ~ dbeta(1,1)
    phi ~ dbeta(1,1)

    }"
    ,fill=TRUE)
sink()