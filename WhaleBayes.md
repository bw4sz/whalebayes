# Whale Bayes: On seperating occurrence and behavior in movement ecology
Ben Weinstein  
`r Sys.time()`  



A fundamental goal in ecology is to attribute the movement of animals across space and time to ecological mechanisms. Animals travel to find food, mates, shelter, and predator free space. Using individual data to parameterize movement models, ecologist have gained insight into animal distribution, migration, and behavior.Pinning down the causes and predictors of animal movement remains a challenging task for ecologists and conservation managers in a changing world. 

While there has been immense focus on describing movement processes, mechanisms, autocorrelation structures, and behavioral phases, it remains difficult to distinguish movement mechanisms across space and time. It is natural to assume that movement phases ('foraging', 'traveling', 'resting') are functions of the environment and life-history strategies. One glaring challenge is partitioning  species presence and behavior. Given that species must be present to participate in behavior, it is natural to assume that two processes are invariably linked. However, when we begin to think about predicting species behavior, we may spuriously conflate the predictors of animal behavior with the predictor of species presence. To date, all movement models assume that presence is given, and seek to extract the environment signatures of behavior based on logistic funtions and markov-models of state dependence. Our aim to step back and relate the huge body of movement ecology literature with the equally well-developed literature of species distribution modeling. We were inspired by the recent paper (Gravel) seeking to biotic interactions and distributions at scale. Our conceptual framework is based on Bayes Rule and statements of conditional probability.

Consider a grid of cells

<img src="WhaleBayes_files/figure-html/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

Species presence at a given cell X is

$$ P(X=1) = 0.5 $$

The probability of species X existing in one of two behavior phases is:

$$ P(B=1) = 0.5 $$

Since we assume species movement stem from just two behavioral phases, the converse probability is:

$$ P(B=1) = 0.5 $$

Therefore the conditional probability of observing a species in cell X = [x,y] existing in behavioral phase B=1 is:

$$ P(B=1|X=1) = P(B=1) * B(X=1)$$

We can therefore model the joint probability of occurrence and behavior as arriving from seperate functions.

### Liklihood Analysis

Since we have a mixture of two processes, we can model the outcome occurrence and behavior state as a binomial and multinomial mixture model. For the sake of simplicity, we model just the probability of state = 1.


$$ X = Bernoulli(z)  $$
$$ z= \phi * \rho$$
$$ Occurrence \sim Bernoulli(\phi)$$
$$ Behavior \sim Bernoulli(\rho)$$


## Example 1: Random occurrence with random behavior

<img src="WhaleBayes_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />


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


```
##    user  system elapsed 
##    0.06    0.02    3.08
```

<img src="WhaleBayes_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-5-2.png" style="display: block; margin: auto;" />

While this will work for simple cases, clearly as we see more complex functions, the probability of occurrence (phi) and the probability of behavior == 1 (rho) will become unidentifiable, and the liklihood landscape will not converge. 

The key value of interest is z, the conditional probability of foraging given occurrence. Here the mean estimate is 0.76, very close to the true known value of 0.25. As we increase the grid size (more data), we would converge on the true answer.

## Example 2: Environmentally dependent occurrence with random behavior

Moving towards the aim of analysis, let's continue with environmentally dependent occurrence, but environmentally independent behavior. We do not claim that behavior is itself random, but that the marginal probabilities of behavior are 0.5 with respect to the environmental conditions of the cell.

$$ X = Bernoulli(z)  $$
$$ z= \phi * \rho$$
$$ logit(\phi) = \alpha + \beta * environment $$

<img src="WhaleBayes_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-6-2.png" style="display: block; margin: auto;" />


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


```
##    user  system elapsed 
##    0.70    0.08   31.03
```

<img src="WhaleBayes_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

## Example 3: Environmentally dependent occurrence with environmentally dependent behavior

Moving towards the aim of analysis, let's continue with environmentally dependent occurrence, but environmentally independent behavior. We do not claim that behavior is itself random, but that the marginal probabilities of behavior are 0.5 with respect to the environmental conditions of the cell.

$$ X = Bernoulli(z)  $$
$$ z= \phi * \rho$$
$$ logit(\phi) = \alpha + \beta * environment $$
$$ logit(\rho) = \alpha_2 + \beta_2 * environment_2 $$

Where the environmental predictors for occurrence and behavior are different variables (eg. temperature, bathymetry).

<img src="WhaleBayes_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-9-2.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-9-3.png" style="display: block; margin: auto;" />


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


```
##    user  system elapsed 
##    0.03    0.00   44.16
```

<img src="WhaleBayes_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" /><img src="WhaleBayes_files/figure-html/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

The critical thing to notice here is that while were able to parameterize a occurrence function (true state in dashed red lines), the probability of behavior is broadly centered on 0.5, with wide confidence intervals. This alerts us that the behavior itself has little environmental influence. 

The critical thing to notice here is that while were able to parameterize a occurrence function (true state in dashed red lines), the probability of behavior is broadly centered on 0.5, with wide confidence intervals. This alerts us that the behavior itself has little environmental influence. 

We want to stress that in no way are we degraded the very excellent work by many movement ecologists, and especially the designers of the moveHMM package. The model being tested was not developed with such a mixture in mind. Nor are we claiming that the divide among occurrence and behavior was unknown. Instead we aim to link two parts of ecological thought, which we hope will motivate continued development within these communities. 

