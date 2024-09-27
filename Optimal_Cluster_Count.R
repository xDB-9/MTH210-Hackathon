gamma_ick <- function(x, mu, sig2, pis, C)
{
  gamma_prob <- matrix(0, nrow = length(x), ncol = C)
  for(c in 1:C)
  {
    gamma_prob[ ,c]  <- dnorm(x, mean = mu[c], sd = sqrt(sig2[c]))* pis[c]
  }
  
  gamma_prob <- gamma_prob/(rowSums(gamma_prob))
  return(gamma_prob)
}

GMMoneDim <- function(x, C )
{
  pis <- numeric(length = C)
  for(i in 1:C){
    pis[i] <- runif(1)
  }
  pis <- pis/sum(pis)
  mu <- seq(min(x), max(x), length.out = C)
  sig2 <- rep(var(x), C)
  diff <- 100
  tol <- 1e-5
  iter <- 0
  
  current <- c(pis, mu, sig2)
  store <- current
  
  while(diff > tol)
  {
    previous <- current
    iter <- iter + 1
    
    # E step: find gamma_{i,c,k} for just c = 1, since for c = 2 is just 1-Ep
    # Ep <- current[1]*dnorm(x, current[2], sqrt(current[4]))/
    #   (current[1]*dnorm(x, current[2], sqrt(current[4])) + (1 - current[1])*dnorm(x, current[3], sqrt(current[5])))
    # 
    Ep <- gamma_ick(x, mu, sig2, pis, C)
    
    # M-step
    pis <- colMeans(Ep)
    mu <- colSums(Ep*x) / colSums(Ep)
    for(c in 1:C)
    {
      sig2[c] <- sum(Ep[,c]*(x - mu[c])^2) / sum(Ep[,c])
    }
    current <- c(pis, mu, sig2)
    
    diff <- norm(previous - current, "2")
    store <- rbind(store, current)
  }
  
  print(current)# final estimates
  
  
  # Final estimates of the probability
  # that each observation is in Class C.
  Prob.Z <- gamma_ick(x, mu = mu, sig2 = sig2, pis = pis, C)
  Prob.Z <- apply(Prob.Z, 1, which.max)
  rtn <- list(Prob.Z, mu, sig2, pis)
  return(rtn)
}

y <- assign4_train$y
erup2 <- GMMoneDim(y, C = 2)
erup3 <- GMMoneDim(y, C = 3)

log.like <- function(y, mu, sig2, pis) {
  n <- length(y)
  logL <- 0
  for (i in 1:n) {
    logL <- logL + log(sum(pis * dnorm(y[i], mu, sqrt(sig2))))
  }
  return(logL)
}

n <- length(y)
loglike2 <- log.like(y, mu = erup2[[2]], sig2 = erup2[[3]], pis = erup2[[4]])
loglike3 <- log.like(y, mu = erup3[[2]], sig2 = erup3[[3]], pis = erup3[[4]])
# number of parameters
C <- 2
K2 <- 3*C - 1
C <- 3
K3 <- 3*C - 1
BIC.2 <- -2*loglike2 + K2*log(n)
BIC.3 <- -2*loglike3 + K3*log(n)

BIC_scores <- c(BIC.2,BIC.3)
iter <- length(BIC_scores)
for(i in 1:length(BIC_scores)){
  if(BIC_scores[i]<BIC_scores[iter])
    iter = i
}

cluster_classif <- list(erup2[[1]], erup3[[1]])
optimal_cluster <- data.frame(cluster_classif[[iter]])

n <- unique(optimal_cluster$cluster_classif..iter..)[1]
mat <- c()
vec1 <- which(optimal_cluster==1)
vec2 <- which(optimal_cluster==2)

cluster1 <- assign4_train[c(vec1),]
cluster2 <- assign4_train[c(vec2),]
