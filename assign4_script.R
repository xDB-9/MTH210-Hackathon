### Paste all your codes for model building
### and cross-validation here

### Paste all your codes for model building
### and cross-validation here
library(readr)
assign4_train <- read_csv("assign4_train.csv")

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


y <- assign4_train$y
x <- data.matrix(assign4_train[,-c(1)])

p <- dim(x)[2]

##############################
#cluster 1

y1 <- cluster1$y
x1 <- data.matrix(cluster1[,-c(1)])
n1 <- dim(x1)[1]
p1 <- dim(x1)[2]
lam.vec1 <- c(10^(seq(-8, -6, by = .01)))
# Will store CV error in this.
CV.error1 <- numeric(length = length(lam.vec1))
for(l in 1:length(lam.vec1))
{
  track.cv1 <- 0
  lam1 <- lam.vec1[l]
  for(i in 1:n1)
  {
    # Making training data
    x1.train <- x1[-i,] # removing ith X
    y1.train <- y1[-i]  #removing ith y
    # fitting model for training data
    beta.train1 <- solve(t(x1.train) %*% x1.train  + diag(lam1, p1)) %*% t(x1.train) %*% y1.train
    # test error
    track.cv1 <- track.cv1 + sum((y1-x1%*%beta.train1)^2)/p1
  }
  CV.error1[l] <- track.cv1/n1
}
chosen.lam1 <- lam.vec1[which.min(CV.error1)]
beta.final1 <- solve(t(x1) %*% x1 + chosen.lam1*diag(p1)) %*% t(x1) %*% y1

# min_value1 <- min(norm_pos1)
# min_index1 <- which.min(norm_pos1)
# norm_pos1_df <- data.frame(norm_pos1)
# beta1_min <- beta1[min_index1,]
y1_pred <- x1%*%beta.final1
mse1 <- sum((y1-y1_pred)^2)/p1

##############################
#cluster 2

y2 <- cluster2$y
x2 <- data.matrix(cluster2[,-c(1)])

n2 <- dim(x2)[1]
p2 <- dim(x2)[2]

lam.vec2 <- c(10^(seq(-8, -6, by = .01)))
# Will store CV error in this.
CV.error2 <- numeric(length = length(lam.vec2))
for(l in 1:length(lam.vec2))
{
  track.cv2 <- 0
  lam2 <- lam.vec2[l]
  for(i in 1:n2)
  {
    # Making training data
    x2.train <- x2[-i,] # removing ith X
    y2.train <- y2[-i]  #removing ith y
    # fitting model for training data
    beta.train2 <- solve(t(x2.train) %*% x2.train  + diag(lam2, p2)) %*% t(x2.train) %*% y2.train
    # test error
    track.cv2 <- track.cv2 + sum((y2-x2%*%beta.train2)^2)/p2
  }
  CV.error2[l] <- track.cv2/n2
}
chosen.lam2 <- lam.vec2[which.min(CV.error2)]
beta.final2 <- solve(t(x2) %*% x2 + chosen.lam2*diag(p2)) %*% t(x2) %*% y2

# min_value2 <- min(norm_pos2)
# min_index2 <- which.min(norm_pos2)
# norm_pos2_df <- data.frame(norm_pos2)
# beta2_min <- beta2[min_index2,]
y2_pred <- x2%*%beta.final2
mse2 <- sum((y2-y2_pred)^2)/p2

####################################
#Conditional Distribution
####################################

p_pos <- seq(0,1,by=1e-4)
beta <- matrix(0,ncol = dim(x)[2], nrow = length(p_pos))
y_pred <- matrix(0,ncol = dim(x)[1], nrow = length(p_pos))
mse <- numeric(length = length(p_pos))

min_mse = 1e5
min_p_pos_index <- 0
beta_min_p_pos <- numeric(length = dim(x)[2])
for(i in 1:length(p_pos)){
  beta[i,] <- (1-p_pos[i])*beta.final1 + p_pos[i]*beta.final2
  y_pred[i,] <- x%*%beta[i,]
  mse[i] <- sum((y-y_pred[i,])^2)/p
  if(mse[i]<min_mse) {
    min_mse = mse[i]
    min_p_pos_index = i
    beta_min_p_pos <- beta[i,]
  }
}

####################################
#Test



y_prediction <- x%*%beta_min_p_pos
mse_final <- sum((y-y_prediction)^2)/500
save(beta_min_p_pos,file = "fit_params.RData")





