

# par(mfrow=c(1,2))
# hist(y1)
# hist(y2)


y <- assign4_train$y
x <- data.matrix(assign4_train[,-c(1)])

p <- dim(x)[2]

##############################
#cluster 1

y1 <- cluster1$y
x1 <- data.matrix(cluster1[,-c(1)])

p1 <- dim(x1)[2]
k_pos1 <- seq(1e-7,1e-4,by=1e-8)

beta1 <- matrix(c(rep(0,length(k_pos1)*dim(x1)[2])), nrow=length(k_pos1),ncol = dim(x1)[2])
y_pred1 <- matrix(c(rep(0,length(k_pos1)*dim(x1)[1])), nrow=length(k_pos1),ncol = dim(x1)[1])
norm_pos1 <- numeric(length = length(y1))
for(i in 1:length(k_pos1)){
  beta1[i,] <- solve(t(x1) %*% x1  + diag(k_pos1[i], p1)) %*% t(x1) %*% y1
  y_pred1[i,] <- x1%*%beta1[i,]
  norm_pos1[i] <- norm(y1-y_pred1[i,],"2")
}

min_value1 <- min(norm_pos1)
min_index1 <- which.min(norm_pos1)
norm_pos1_df <- data.frame(norm_pos1)
beta1_min <- beta1[min_index1,]
y1_pred <- x1%*%beta1_min
mse1 <- sum((y1-y1_pred)^2)/p1

##############################
#cluster 2

y2 <- cluster2$y
x2 <- data.matrix(cluster2[,-c(1)])

p2 <- dim(x2)[2]
k_pos2 <- seq(1e-4,1e-1,by=1e-5)

beta2 <- matrix(c(rep(0,length(k_pos2)*dim(x2)[2])), nrow=length(k_pos2),ncol = dim(x2)[2])
y_pred2 <- matrix(c(rep(0,length(k_pos2)*dim(x2)[1])), nrow=length(k_pos2),ncol = dim(x2)[1])
norm_pos2 <- numeric(length = length(y2))
for(i in 1:length(k_pos2)){
  beta2[i,] <- solve(t(x2) %*% x2  + diag(k_pos2[i], p2)) %*% t(x2) %*% y2
  y_pred2[i,] <- x2%*%beta2[i,]
  norm_pos2[i] <- norm(y2-y_pred2[i,],"2")
}

min_value2 <- min(norm_pos2)
min_index2 <- which.min(norm_pos2)
norm_pos2_df <- data.frame(norm_pos2)
beta2_min <- beta2[min_index2,]
y2_pred <- x2%*%beta2_min
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
  beta[i,] <- (1-p_pos[i])*beta1_min + p_pos[i]*beta2_min
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


save(beta,file = "fit_params.RData")




