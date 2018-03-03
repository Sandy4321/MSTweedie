#library(Matrix)
#library(tweedie)
#library(mvtnorm)
#library(HDtweedie)
#dyn.load("R/MTtweedie.so")


gen.data.inf <- function(K,p,df,n){
   ##generate model
   ids <- sample(1:p, df)
   tasks <- sample(1:K, df, replace = TRUE)
   model <- matrix(0, p, K)
   for(i in seq(df)){
      model[ids[i], tasks[i]] <- 1
   }
   #beta <- t(matrix(rnorm(p*K),p,K)*model)
   beta <- t(matrix(3,p,K)*model)
   ##matrix of covariances
   s1 <- matrix(rnorm(p*p),p)
   sigma <- s1 %*% t(s1)
   sigma <- sigma / p
   ##generate covariates
   x <- lapply(n, function(n){
      rmvnorm(n, sigma = sigma)
   })
   ##generate responses
   y <- lapply(1:K, function(k){
      mu <- t(beta[k,]) %*% t(x[[k]])
      rtweedie(n[k], power = 1.5, mu = exp(mu[1,]),phi = 1)
   })
   list(x = x, y = y, beta = beta, K=K, p=p, df = df, n = n, sigma = sigma,
        model = cbind(tasks = tasks, var = ids))
}

gen.data.2 <- function(K,p,df,n){
   ##generate model
   ids <- sample(1:p, df)
   model <- matrix(0, p, K)
   for(i in seq(df)){
      model[ids[i], ] <- 1
   }
   #beta <- t(matrix(rnorm(p*K),p,K)*model)
   beta <- t(matrix(2,p,K)*model)
   beta0 <- rnorm(K)
   ##matrix of covariances
   s1 <- matrix(rnorm(df*df),df)
   sigma <- s1 %*% t(s1)
   sigma <- sigma / df
   ##generate covariates
   x <- lapply(n, function(n){
      good <- rmvnorm(n, sigma = sigma)
      mat <- matrix(runif(n*p, -1, 1), nrow = n, ncol = p)
      mat[,ids] <- good
      mat
   })
   ##generate responses
   y <- lapply(1:K, function(k){
      mu <- t(beta[k,]) %*% t(x[[k]]) + beta0[K]
      rtweedie(n[k], power = 1.5, mu = exp(mu[1,]), phi = 1)
   })
   list(x = x, y = y, beta = beta, K=K, p=p, df = df, n = n, sigma = sigma,
        beta0 = beta0, model = ids)
}


##############################
# tests
##############################

##parameters
K <- 5
n <- rep(60,K)
p <- 10
df <- 3

data <- gen.data.inf(K, p, df, n)
fit <- MTtweedie(data$x, data$y, reg = 'Linf')
fit <- MTtweedie(data$x, data$y, reg = 'L2')

kkt.check(fit, cond=c(1:4))

par(mfrow=c(2,2), mar = c(2,2,5,1))

data <- gen.data.inf(K, p, df, n)
fit <- MTtweedie(data$x, data$y, reg = 'Linf')
descent.plot(fit, main = 'Linf')
fit$time
fit <- MTtweedie(data$x, data$y, reg = 'L2')
descent.plot(fit, main= 'L2')
fit$time

fit

u <- ceiling(sqrt(K+1))
v <- ceiling((K+1)/u)
par(mfrow=c(u,v), mar = c(2,4,1,1))

plot.MTtweedie(x = fit, type.coef = 'coef', log.lambda = TRUE)
plot.MTtweedie(x = fit, type.coef = 'norm', log.lambda = TRUE)
data$model

######################
## CROSS VALIDATION
######################

data <- gen.data.2(K, p, df, n)

cv.inf <- cv.MTtweedie(data$x, data$y, kktstop = F,
                       reg = 'Linf', lambda.min=1e-4)
cv.2 <- cv.MTtweedie(data$x, data$y, kktstop = F,
                     reg = 'L2', lambda.min=1e-4)

par(mfrow=c(2,2), mar = c(2,2,2,1))
plot.cv.MTtweedie(cv.inf, main = 'Linf')
plot.MTtweedie(x = cv.inf$MTtweedie.fit, type.coef = 'norm', log.lambda = TRUE)
plot.MTtweedie(x = cv.2$MTtweedie.fit, type.coef = 'coef', log.lambda = TRUE)
plot.cv.MTtweedie(cv.2, main = 'L2')
plot.MTtweedie(x = cv.2$MTtweedie.fit, type.coef = 'norm', log.lambda = TRUE)

######################
## INPUT TYPES
######################

data <- dat[seq(1000),]

fit <- MTtweedie(x=data, y=1, source=4, reg = 'Linf')

plot.MTtweedie(x = cv$MTtweedie.fit, type.coef = 'coef', log.lambda = TRUE)
plot.cv.MTtweedie(cv)
plot.MTtweedie(x = cv$MTtweedie.fit, type.coef = 'norm', log.lambda = TRUE)

