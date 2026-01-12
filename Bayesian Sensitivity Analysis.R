# Simulation,
# Sample size to consider 200, 500, 1000 patients;
# y, an end-of-study binary outcome;
# A, a binary treatment assignment;
# U, a binary unmeasured confounder
#    i) independent of measured confounders;
#    ii)* interactions with two measured confounders;
#    to consider degrees of unmeasured confounding low, mid, high;
# X, a vector of measured confounders, 1 binary and 1 continuous;

if( !file.exists("simulation_3_visit_info_prior_jags.txt") ){
  
  cat( " model {
  for (i in 1:n){
  Y[i] ~ dnorm(py[i], tau)
  py[i] <- b30 + b31 * X1[i] + b32 * X2[i] + b33* X3[i] + b34 * U1[i] +  b35 * U2[i] + b36 * U3[i] + b37 * A1[i] + b38 * A2[i] + b39 * A3[i]
  
  X1[i] ~ dbern(theta_x1[i])
  logit(theta_x1[i]) <- t10
    
  U1[i] ~ dbern(theta_u1[i])
  logit(theta_u1[i]) <- r10  
  
  A1[i] ~ dbern(theta_a1[i])
  logit(theta_a1[i]) <- a10 + a11 * X1[i] + a12 * U1[i]
  
  X2[i] ~ dbern(theta_x2[i])
  logit(theta_x2[i]) <- t20 + t21 * X1[i] + t22 * A1[i]
  
  U2[i] ~ dbern(theta_u2[i])
  logit(theta_u2[i]) <- r20 + r21 * U1[i] + r22 * A1[i]
  
  A2[i] ~ dbern(theta_a2[i])
  logit(theta_a2[i]) <- a20 + a21 * X1[i] + a22 * X2[i] + a23 * U1[i] + a24 * U2[i] + a25 * A1[i]
  
  X3[i] ~ dbern(theta_x3[i])
  logit(theta_x3[i]) <- t30 + t31 * X1[i] + t32 * X2[i] + t33 * A1[i] + t34 * A2[i]
  
  U3[i] ~ dbern(theta_u3[i])
  logit(theta_u3[i]) <- r30 + r31 * U1[i] + r32 * U2[i] + r33 * A1[i] + r34 * A2[i] 
  
  A3[i] ~ dbern(theta_a3[i])
  logit(theta_a3[i]) <- a30 + a31 * X1[i] + a32 * X2[i] + a33 * X3[i] + a34 * U1[i] + a35 *U2[i] + a36 * U3[i] + a37 * A2[i] + a38 * A1[i]
  
 }
 ### Priors
    b30 ~ dnorm(0, 0.01);
    b31 ~ dnorm(0, 0.1);
    b32 ~ dnorm(0, 0.1);
    b33 ~ dnorm(0, 0.1);
    b37 ~ dnorm(0, 0.1);
    b38 ~ dnorm(0, 0.1);
    b39 ~ dnorm(0, 0.1);
    

    t10 ~ dnorm(0, 0.01);
    t20 ~ dnorm(0, 0.01);
    t21 ~ dnorm(0, 0.1);
    t22 ~ dnorm(0, 0.1);
    t30 ~ dnorm(0, 0.01);
    t31 ~ dnorm(0, 0.1);
    t32 ~ dnorm(0, 0.1);
    t33 ~ dnorm(0, 0.1);
    t34 ~ dnorm(0, 0.1);
    
    
    a10 ~ dnorm(0, 0.1);
    a11 ~ dnorm(0, 0.1);
    a20 ~ dnorm(0, 0.01);
    a21 ~ dnorm(0, 0.1);
    a22 ~ dnorm(0, 0.1);
    a25 ~ dnorm(0, 0.1);
    a30 ~ dnorm(0, 0.01);
    a31 ~ dnorm(0, 0.1);
    a32 ~ dnorm(0, 0.1);
    a33 ~ dnorm(0, 0.1);
    a37 ~ dnorm(0, 0.1);
    a38 ~ dnorm(0, 0.1);
    
    tau ~ dgamma(0.1,0.1);
    
    # bias parameter
    b34 ~ dunif(-10, 10);
    b35 ~ dunif(-10, 10);
    b36 ~ dunif(-10, 10);
    
    a12 ~ dunif(-2, 2);
    a23 ~ dunif(-2, 2);
    a24 ~ dunif(-2, 2);
    a34 ~ dunif(-2, 2);
    a35 ~ dunif(-2, 2);
    a36 ~ dunif(-2, 2);
    
    
    r10 ~ dunif(-5, 5);
    r20 ~ dunif(-5, 5);
    r21 ~ dunif(-2, 2);
    r22 ~ dunif(-2, 2);
    r30 ~ dunif(-5, 5);
    r31 ~ dunif(-2, 2);
    r32 ~ dunif(-2, 2);    
    r33 ~ dunif(-2, 2);
    r34 ~ dunif(-2, 2);
} ", file = "simulation_3_visit_info_prior_jags.txt")
  
}
#parallel;
library(parallel)
library(foreach)
library(doParallel)

registerDoParallel(40)


mysim.bin.cens <- function(outfile, from=1, to=10, U=1, #parameters defining the simulation iterations and outcome;
                           samplesize=500, nboot = 500,#parameters for data simulation, either 500/ 1000
                           #1. no unmeasured confounder; 2. time-constant unmeasured confounder; 3. time varying unmeasured confounder
                           B=500, pdraw = 4000, nburnin = 3000,#parameters for analysis;
                           Px11 = .5, Px12 = .4, #simulate measured covariates;
                           coef_u = -.3, coef_a = .1
) {

  expit <- function(x){
    x <- exp(x)/(exp(x)+1)
    return(x)
}

  results_run<-foreach(i=from:to, .combine='rbind',.inorder=T, .verbose=T) %dopar% {
  library(tidyverse)
  library(R2jags)
  # data generation seed;

  set.seed(i+123)

  results.it <- matrix(NA, 1, 9) #to change number of parameters;
  colnames(results.it)<- c("true value","MSM withU -est","MSM withU -sd","MSM noU -est","MSM noU -sd",
                           "BSA -est","BSA -sd","BSA -CIL","BSA -CLU")
  X1_1 <- rbinom(samplesize,1, prob = 0.5)
  X2_1 <- rbinom(samplesize,1, prob = 0.5)

  # you can specify coef_UZ1 = coef_UZ visit;
  A1 <- rbinom(samplesize, 1, expit(1 - 0.5*X1_1+0.5*X2_1))

  ## Second Visit
  X1_2 = rbinom(samplesize,1, prob = expit(0.5*X1_1-0.5*A1))
  X2_2 = rbinom(samplesize,1, prob = expit(0.5*X2_1+0.5*A1)) #does not dependent on X1_2
  A2 <- rbinom(samplesize, 1, expit(1 - 0.3*A1-0.3*X1_1+0.4*X2_1-0.4*X1_2+0.5*X2_2))

  ## Third Visit
  X1_3 = rbinom(samplesize,1, prob = expit(0.5*X1_2 - 0.5*A2))
  X2_3 = rbinom(samplesize,1, prob = expit(0.5*X2_2 + 0.5*A2))

  A3 <- rbinom(samplesize, 1, expit(1 - 0.3*A2-0.3*X1_1+0.4*X2_1-0.4*X1_2+0.5*X2_2-0.5*X1_3+0.6*X2_3))

  Y111 = rnorm(samplesize, mean = 10 - 4-3-2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a2,a1
  Y110 = rnorm(samplesize, mean = 10 - 4-3  +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a2
  Y100 = rnorm(samplesize, mean = 10 - 4    +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3
  Y101 = rnorm(samplesize, mean = 10 - 4  -2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a1
  Y011 = rnorm(samplesize, mean = 10 -   3-2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #    a2 a1
  Y010 = rnorm(samplesize, mean = 10 -   3  +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #    a2
  Y000 = rnorm(samplesize, mean = 10        +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #
  Y001 = rnorm(samplesize, mean = 10 -     2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #       a1
  
  Y = Y111*(A1)*(A2)*(A3) + Y110*(1-A1)*(A2)*(A3)+ Y100*(1-A1)*(1-A2)*(A3) + Y101*(A1)*(1-A2)*(A3)+ 
    Y011*(A1)*(A2)*(1-A3) + Y010*(1-A1)*(A2)*(1-A3)+ Y000*(1-A1)*(1-A2)*(1-A3) + Y001*(A1)*(1-A2)*(1-A3)
  
  # Y <- rnorm(samplesize, mean = 10 -4*A3-3*A2-2*A1+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2)

  mytrue <- function(A1s=1, A2s=1, A3s=1){
    #visit 1;
    X1_1s <- rbinom(samplesize*10,1, prob = 0.5)
    X2_1s <- rbinom(samplesize*10,1, prob = 0.5)
    ## Second Visit
    X1_2s = rbinom(samplesize*10,1, prob = expit(0.5*X1_1s-0.5*A1s))
    X2_2s = rbinom(samplesize*10,1, prob = expit(0.5*X2_1s+0.5*A1s))
    ## Third Visit
    X1_3s = rbinom(samplesize*10,1, prob = expit(0.5*X1_2s-0.5*A2s))
    X2_3s = rbinom(samplesize*10,1, prob = expit(0.5*X2_2s+0.5*A2s))
    #Visit 3, simulate potential outcomes;
    truemean = (10-4*A3s-3*A2s-2*A1s+1*X1_1s+2*X1_2s+3*X1_3s-1*X2_1s-2*X2_2s-3*X2_3s)
    return(mean(truemean))
  }
  
  generate <- function(samplesize){
    X1_1 <- rbinom(samplesize,1, prob = 0.5)
    X2_1 <- rbinom(samplesize,1, prob = 0.5)
    
    # you can specify coef_UZ1 = coef_UZ visit;
    A1 <- rbinom(samplesize, 1, expit(1 - 0.5*X1_1+0.5*X2_1))
    
    ## Second Visit
    X1_2 = rbinom(samplesize,1, prob = expit(0.5*X1_1-0.5*A1))
    X2_2 = rbinom(samplesize,1, prob = expit(0.5*X2_1+0.5*A1)) #does not dependent on X1_2
    A2 <- rbinom(samplesize, 1, expit(1 - 0.3*A1-0.3*X1_1+0.4*X2_1-0.4*X1_2+0.5*X2_2))
    
    ## Third Visit
    X1_3 = rbinom(samplesize,1, prob = expit(0.5*X1_2 - 0.5*A2))
    X2_3 = rbinom(samplesize,1, prob = expit(0.5*X2_2 + 0.5*A2))
    
    A3 <- rbinom(samplesize, 1, expit(1 - 0.3*A2-0.3*X1_1+0.4*X2_1-0.4*X1_2+0.5*X2_2-0.5*X1_3+0.6*X2_3))
    
    Y111 = rnorm(samplesize, mean = 10 - 4-3-2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a2,a1
    Y110 = rnorm(samplesize, mean = 10 - 4-3  +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a2
    Y100 = rnorm(samplesize, mean = 10 - 4    +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3
    Y101 = rnorm(samplesize, mean = 10 - 4  -2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) # a3,a1
    Y011 = rnorm(samplesize, mean = 10 -   3-2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #    a2 a1
    Y010 = rnorm(samplesize, mean = 10 -   3  +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #    a2
    Y000 = rnorm(samplesize, mean = 10        +1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #
    Y001 = rnorm(samplesize, mean = 10 -     2+1*X1_1+2*X1_2+3*X1_3-1*X2_1-2*X2_2-3*X2_3, sd = 2) #       a1
    
    Y = Y111*(A1)*(A2)*(A3) + Y110*(1-A1)*(A2)*(A3)+ Y100*(1-A1)*(1-A2)*(A3) + Y101*(A1)*(1-A2)*(A3)+ 
      Y011*(A1)*(A2)*(1-A3) + Y010*(1-A1)*(A2)*(1-A3)+ Y000*(1-A1)*(1-A2)*(1-A3) + Y001*(A1)*(1-A2)*(1-A3)
    full_dat <- data.frame(X1_1, X2_1, X1_2, X2_2, X1_3, X2_3, A1,A2,A3,Y111,Y110,Y101,Y100,Y011,Y010,Y001,Y000,Y,id=c(1:samplesize))
    return(full_dat)
  }
  
     generateMatrix <- function(n) {
       combos <- expand.grid(replicate(n, c(0, 1), simplify = FALSE))
       combos <- as.data.frame(combos) %>%
         # intercept
         mutate(int = 1, .before = Var1)
       return(as.matrix(combos))
     }

  # saving final data
  # Step 6: put together the observed;
  full_dat <- data.frame(X1_1, X2_1, X1_2, X2_2, X1_3, X2_3, A1,A2,A3,Y111,Y110,Y101,Y100,Y011,Y010,Y001,Y000,Y,id=c(1:samplesize))

  
  results.it[1,1] <- mytrue(A1s=1,A2s=1,A3s=1) - mytrue(A1s=0,A2s=0,A3s=0)
# # 5. Analysis: Full Bayesian Sensitivity analysis;
#       # Create a list of data
   data_model <- list(
     Y = Y,
     n = length(Y),
     A1 = A1,
     A2 = A2,
     A3 = A3,
     X1 = X2_1,
     X2 = X2_2,
     X3 = X2_3
   )
    parameters <- c("t10", "t20", "t21","t22", "t30","t31","t32","t33","t34","r10", "r20", "r21","r22","r30","r31","r32","r33","r34","b30","b31","b32","b33","b34","b35","b36","b37","b38","b39","tau")
    model<- R2jags::jags(model.file = "simulation_3_visit_info_prior_jags.txt",
                         parameters.to.save = parameters,
                         data= data_model,
                         # jags.seed = i+123,
                         n.chains = 4,
                         n.iter = pdraw,
                         n.burnin = nburnin,
                         n.thin = 3)


   #### Probability calculation
   coef_x1_1 = cbind(model$BUGSoutput$sims.list$r10)
   coef_x2_1 = cbind(model$BUGSoutput$sims.list$t10)
   # coef_a1 = cbind(model$BUGSoutput$sims.list$a10, model$BUGSoutput$sims.list$a11, model$BUGSoutput$sims.list$a12)
   coef_x1_2 = cbind(model$BUGSoutput$sims.list$r20, model$BUGSoutput$sims.list$r21, model$BUGSoutput$sims.list$r22)
   coef_x2_2 = cbind(model$BUGSoutput$sims.list$t20, model$BUGSoutput$sims.list$t21, model$BUGSoutput$sims.list$t22)
   # coef_a2 = cbind(model$BUGSoutput$sims.list$a20, model$BUGSoutput$sims.list$a21, model$BUGSoutput$sims.list$a22, model$BUGSoutput$sims.list$a23, model$BUGSoutput$sims.list$a24, model$BUGSoutput$sims.list$a25)
   coef_x1_3 = cbind(model$BUGSoutput$sims.list$r30,model$BUGSoutput$sims.list$r31,model$BUGSoutput$sims.list$r32, model$BUGSoutput$sims.list$r33, model$BUGSoutput$sims.list$r34)
   coef_x2_3 = cbind(model$BUGSoutput$sims.list$t30,model$BUGSoutput$sims.list$t31,model$BUGSoutput$sims.list$t32, model$BUGSoutput$sims.list$t33, model$BUGSoutput$sims.list$t34)
   # coef_a3 = cbind(model$BUGSoutput$sims.list$a30, model$BUGSoutput$sims.list$a31, model$BUGSoutput$sims.list$a32, model$BUGSoutput$sims.list$a33, model$BUGSoutput$sims.list$a34, model$BUGSoutput$sims.list$a35,model$BUGSoutput$sims.list$a36, model$BUGSoutput$sims.list$a37, model$BUGSoutput$sims.list$a38)
   coef_y = cbind(model$BUGSoutput$sims.list$b30, model$BUGSoutput$sims.list$b31, model$BUGSoutput$sims.list$b32, model$BUGSoutput$sims.list$b33, model$BUGSoutput$sims.list$b34, model$BUGSoutput$sims.list$b35,model$BUGSoutput$sims.list$b36, model$BUGSoutput$sims.list$b37, model$BUGSoutput$sims.list$b38, model$BUGSoutput$sims.list$b39)
#
   generateMatrix <- function(n) {
     combos <- expand.grid(replicate(n, c(0, 1), simplify = FALSE))
     combos <- as.data.frame(combos) %>%
       # intercept
       mutate(int = 1, .before = Var1)
     return(as.matrix(combos))
   }
   # calculate probability
   #first visit
   p_U1 = cbind(1 - expit(coef_x1_1), expit(coef_x1_1)) # P(U1)
   p_X1 = cbind(1 - expit(coef_x2_1), expit(coef_x2_1)) # P(X1)
   # second visit
   iden2 = generateMatrix(ncol(coef_x1_2)-2)
   iden_2_a1 = as.matrix(cbind(iden2,matrix(1, nrow = nrow(iden2), ncol = 1)))
   p_X2_1 =  cbind((1 - expit(coef_x2_2 %*% t(iden_2_a1))) * p_X1, expit(coef_x2_2 %*% t(iden_2_a1)) * p_X1) # P(X2,X1)
   p_U2_1 =  cbind((1 - expit(coef_x1_2 %*% t(iden_2_a1))) * p_U1, expit(coef_x1_2 %*% t(iden_2_a1)) * p_U1) # P(U2,U1)
   iden_2_a0= as.matrix(cbind(iden2,matrix(0, nrow = nrow(iden2), ncol = 1)))
   p_X2_0 =  cbind((1 - expit(coef_x2_2 %*% t(iden_2_a0))) * p_X1, expit(coef_x2_2 %*% t(iden_2_a0)) * p_X1) # P(X2,X1)
   p_U2_0 =  cbind((1 - expit(coef_x1_2 %*% t(iden_2_a0))) * p_U1, expit(coef_x1_2 %*% t(iden_2_a0)) * p_U1) # P(U2,U1)
   # third visit
   iden3 = generateMatrix(ncol(coef_x1_3)-3)
   iden_3_a11 = as.matrix(cbind(iden3,matrix(1, nrow = nrow(iden3), ncol = 2)))
   p_X3_1 = cbind((1 - expit(coef_x2_3 %*% t(iden_3_a11)))* p_X2_1, expit(coef_x2_3 %*% t(iden_3_a11))* p_X2_1)
   p_U3_1 = cbind((1 - expit(coef_x1_3 %*% t(iden_3_a11)))* p_U2_1, expit(coef_x1_3 %*% t(iden_3_a11))* p_U2_1)
   iden_3_a00 = as.matrix(cbind(iden3,matrix(0, nrow = nrow(iden3), ncol = 2)))
   p_X3_0 = cbind((1 - expit(coef_x2_3 %*% t(iden_3_a00)))* p_X2_0, expit(coef_x2_3 %*% t(iden_3_a00))* p_X2_0)
   p_U3_0 = cbind((1 - expit(coef_x1_3 %*% t(iden_3_a00)))* p_U2_0, expit(coef_x1_3 %*% t(iden_3_a00))* p_U2_0)
   # marginal probability
   PX0_expanded <- p_X3_0[, rep(1:ncol(p_X3_0), each=ncol(p_U3_0))]
   PU0_expanded <- p_U3_0[, rep(1:ncol(p_U3_0), ncol(p_X3_0))]
   P_0 <- PX0_expanded * PU0_expanded
   PX1_expanded <- p_X3_1[, rep(1:ncol(p_X3_1), each=ncol(p_U3_1))]
   PU1_expanded <- p_U3_1[, rep(1:ncol(p_U3_1), ncol(p_X3_1))]
   P_1 <- PX1_expanded * PU1_expanded
#
   # identity matrix for Y
   iden_y = generateMatrix(ncol(coef_y) - 4)
   iden_muy0 <- as.matrix(cbind(iden_y, matrix(0, nrow = nrow(iden_y), ncol = 3)))
   iden_muy1 <- as.matrix(cbind(iden_y, matrix(1, nrow = nrow(iden_y), ncol = 3)))
#
   muy1 = rowSums(coef_y %*% t(iden_muy1) * P_1)
   muy0 = rowSums(coef_y %*% t(iden_muy0) * P_0)
#
   results.it[1,6] <- mean(muy1 - muy0)
   results.it[1,7] <- sd(muy1 - muy0)
   results.it[1,8:9] <- quantile(muy1 - muy0,prob=c(0.025, 0.975))

#

#    # combining parallel results;
   cbind(i,results.it)
#   
  }
  save(results_run,file=paste0(outfile,".RData"))
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = TRUE)
  
}

mysim.bin.cens("simulation result for bsa", from=1, to= 1000, U=1, #parameters defining the simulation iterations and outcome;
         samplesize=500,#parameters for data simulation
         B=500, pdraw = 8000, nburnin = 6000,#parameters for analysis;
         Px11 = .5, Px12 = .4, coef_u = -0.2, coef_a = .1)#simulate measured covariates
