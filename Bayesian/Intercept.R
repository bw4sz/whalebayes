sink("Bayesian/Intercept.jags")
cat("
model{

    #Liklihood

    for (cell in 1:cells){
      
      beh[cell] ~ dbern(z[cell])

      z[cell]=X[cell] * rho + 0.00001

      X[cell] ~ dbern(phi)
    }

    m=mean(z[])

    rho ~ dbeta(1,1)
    phi ~ dbeta(1,1)

    }"
    ,fill=TRUE)
sink()