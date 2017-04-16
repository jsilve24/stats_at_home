# Calculate with FTheta and GTheta with (optional) Time Varying Matricies -----

# Given a model this function calculates G_t*Theta_t-1, note that 
# G_t can have a time-varying matrix, also ignores last entry of theta assuming
# this will have 1 more row than mod$X
# model is dlm model object (just a list) with the following components
#     GG, (Optional: JGG, X)
# theta is a matrix with 1 more row than observations in the time-series (first
#     row corresponds to unobserbved initial state)
gtheta <- function(mod, theta){
  if (is.null(mod$JGG)) gt <- theta[-(nrow(theta)), ] %*% t(mod$GG)
  else {
    nz <- mod$JGG != 0
    JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
    gt <- array(0, dim=c(nrow(theta)-1, ncol(theta))) # this is the result
    GG <- mod$GG
    for (i in 1:(nrow(theta)-1)){
      GG[JGG[,-3, drop=FALSE]] <- mod$X[i, JGG[,3]]
      gt[i,] <- theta[i,] %*% t(GG)
    }
  }
  return(gt) 
}

# Given a model this function calculates F_t*Theta_t, note that 
# F_t can have a time-varying matrix, also ignores the first entry of theta 
# assuming this will have 1 more row than mod$X
# model is dlm model object (just a list) with the following components
#     FF, (Optional: JFF, X)
# theta is a matrix with 1 more row than observations in the time-series (first
#     row corresponds to unobserbved initial state)
ftheta <- function(mod, theta){
  if (is.null(mod$JFF)) ft <- theta[-1, ] %*% t(mod$FF)
  else {
    nz <- mod$JFF != 0
    JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
    ft <- array(0, dim=c(nrow(theta)-1, nrow(mod$FF))) # this is the result
    FF <- mod$FF
    for (i in 1:(nrow(theta)-1)){
      FF[JFF[,-3, drop=FALSE]] <- mod$X[i, JFF[,3]]
      ft[i,] <- theta[i+1,] %*% t(FF)
    }
  }
  return(ft) 
}

# function to tidy a multi-dimentional array
tidy_array <- function(a){
  d <- dim(a)
  l <- list()
  for (i in 1:length(d)){
    l[[paste0('dim_',i)]] <- 1:d[i]
  }
  tidy <- expand.grid(l)
  tidy[["var"]] <- a[as.matrix(tidy)]
  tidy[["parameter"]] <- lazyeval::expr_text(a)
  
  return(tidy)
}

# Random draws from Logistic Normal distribution 
# with parameters given in PhILR/ILR space. 
# m is the mean
# Sigma is the covariance matrix
# Tr is an object of class phylo
rlogisticnormal <- function(n, m, Sigma, tr){
  sbp <- phylo2sbp(tr)
  V <- buildilrBasep(sbp, p=c(1,1,1)) # Build contrast matrix witout tayon-weights
  
  r <- MASS::mvrnorm(n, m, Sigma)
  r <- ilrInv(r, V)
  return(unclass(r))
}


