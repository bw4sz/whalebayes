sink("Bayesian/Intercept.jags")
cat("
model{

    #Liklihood

    for (cell in 1:cells){
      
      beh[cell] ~ dbern(X[cell] * rho + 0.0001)

      X[cell] ~ dbern(phi)
    }

    z= rho * phi

    rho ~ dbeta(1,1)
    phi ~ dbeta(1,1)

    }"
    ,fill=TRUE)
sink()