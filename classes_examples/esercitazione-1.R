library(R2jags)

 library(car)
data(Angell)
angell.1 <- Angell[, -4] ## take off the fourth
 ## column (remember, the order is (row, column))
 
 moral <- angell.1$moral
hetero <- angell.1$hetero
mobility <- angell.1$mobility
 N <- length(angell.1$moral)
angell.data <- list("moral", "hetero", "mobility", "N")

angell.params <- c("alpha", "beta1", "beta2","tau")
inits1 <- list("alpha"=0, "beta1"=0, "beta2"=0)
inits2 <- list("alpha"=1, "beta1"=1, "beta2"=1)
angell.inits <- list(inits1, inits2)
set.seed(123)
angellfit <- jags(data=angell.data, 
                  inits=angell.inits, 
                  angell.params, 
                  n.chains=2, 
                  n.iter=9000, 
                  n.burnin=1000, 
                  model.file="angell.model.jags")
traceplot(angellfit)

angellfit<-update(angellfit,n.iter=10000)
round(angellfit$BUGS$summary,2)
alpha<-angellfit$BUGS$sims.matrix[,grep("alpha",colnames(angellfit$BUGS$sims.matrix))]
require(coda)
alpha.m<-mcmc.list(alpha1=mcmc(alpha[1:1000],start=1,end=10000,thin=10), alpha2=mcmc(alpha[1001:2000],start=1,end=10000,thin=10))
#two chains 
plot(alpha.m)
#check convergence
gelman.diag(alpha.m)#computes Potential scale reduction factor for the object
