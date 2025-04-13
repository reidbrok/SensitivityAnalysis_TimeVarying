# Simulation,
# Sample size to consider 200, 500, 1000 patients;
# y, an end-of-study binary outcome;
# A, a binary treatment assignment;
# U, a binary unmeasured confounder
#    i) independent of measured confounders;
#    ii)* interactions with two measured confounders;
#    to consider degrees of unmeasured confounding low, mid, high;
# X, a vector of measured confounders, 1 binary and 1 continuous;
#

options(warn=-1)
#parallel;
library(parallel)
library(foreach)
library(doParallel)

registerDoParallel(10)

samplesize = 500
expit <- function(x){
  x <- exp(x)/(exp(x)+1)
  return(x)
}
mytrue <- function(A1s=1, A2s=1, A3s=1){
  #visit 1;
  X1_1s <- rbinom(samplesize*1000,1, prob = 0.5)
  X2_1s <- rbinom(samplesize*1000,1, prob = 0.5)
  ## Second Visit
  X1_2s = rbinom(samplesize*1000,1, prob = expit(0.5*X1_1s-0.5*A1s))
  X2_2s = rbinom(samplesize*1000,1, prob = expit(0.5*X2_1s+0.5*A1s))
  ## Third Visit
  X1_3s = rbinom(samplesize*1000,1, prob = expit(0.5*X1_2s-0.5*A2s))
  X2_3s = rbinom(samplesize*1000,1, prob = expit(0.5*X2_2s+0.5*A2s))
  #Visit 3, simulate potential outcomes;
  trueY <- rnorm(samplesize*1000,
                 mean = 10-4*A3s-3*A2s-2*A1s+1*X1_1s+2*X1_2s+3*X1_3s-1*X2_1s-2*X2_2s-3*X2_3s, sd = 2)
  return(mean(trueY))
}

mytrue(1,1,1) - mytrue(0,0,0)
# -10.27

mysim.bin.cens <- function(outfile, from=1, to=100, U=1, #parameters defining the simulation iterations and outcome;
                           samplesize=500, nboot = 500,#parameters for data simulation, either 500/ 1000
                           #1. no unmeasured confounder; 2. time-constant unmeasured confounder; 3. time varying unmeasured confounder
                           B=500, pdraw = 4000, nburnin = 3000,#parameters for analysis;
                           Px11 = .5, Px12 = .4, #simulate measured covariates;
                           coef_u = -.3, coef_a = .1
) {


  #
  #
  # from=1; to=4;
  # samplesize=500;
  # B=450;  pdraw = 8000; nburnin = 5000; i= 1; nboot = 500

  expit <- function(x){
    x <- exp(x)/(exp(x)+1)
    return(x)
  }

  results_run<-foreach(i=from:to, .combine='rbind',.inorder=T, .verbose=T) %dopar% {

    library(tidyverse)
    library(coda)
    library(WeightIt) #package to assess covariates balance by treatment;
    library(marginaleffects)
    library(data.table)
    # data generation seed;

    set.seed(i+123)

    results.it <- matrix(NA, 1, 17) #to change number of parameters;
    colnames(results.it)<- c("true value","MSM withU -est","MSM withU -sd","MSM noU -est","MSM noU -sd",
                             "G_compute sf -est","G_compute sf -sd","G_compute sf error -est",
                             "G_compute sf error -sd","G_compute sf error -naive sd",
                             "G_compute sf error h2 -est",
                             "G_compute sf error h2 -sd","G_compute sf h2 error -naive sd",
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

    # saving final data
    # Step 6: put together the observed;
    full_dat <- data.frame(X1_1, X2_1, X1_2, X2_2, X1_3, X2_3,
                           A1,A2,A3,
                           Y111,Y110,Y101,Y100,Y011,Y010,Y001,Y000,Y,id=c(1:samplesize))

    # xtabs(~A1+A2+A3, data = full_dat)
    # write.csv(full_dat, "full_dat.csv")

    # data approximate where the true value is provided at the top at -10.27;
    # results.it[1,1] <- mean(full_dat$Y111) - mean(full_dat$Y000)
    mytrue <- function(A1s=1, A2s=1, A3s=1){
      #visit 1;
      X1_1s <- rbinom(samplesize*1000,1, prob = 0.5)
      X2_1s <- rbinom(samplesize*1000,1, prob = 0.5)
      ## Second Visit
      X1_2s = rbinom(samplesize*1000,1, prob = expit(0.5*X1_1s-0.5*A1s))
      X2_2s = rbinom(samplesize*1000,1, prob = expit(0.5*X2_1s+0.5*A1s))
      ## Third Visit
      X1_3s = rbinom(samplesize*1000,1, prob = expit(0.5*X1_2s-0.5*A2s))
      X2_3s = rbinom(samplesize*1000,1, prob = expit(0.5*X2_2s+0.5*A2s))
      #Visit 3, simulate potential outcomes;
      trueY <- rnorm(samplesize*1000,
                     mean = 10-4*A3s-3*A2s-2*A1s+1*X1_1s+2*X1_2s+3*X1_3s-1*X2_1s-2*X2_2s-3*X2_3s, sd = 2)
      return(mean(trueY))
    }
    results.it[1,1] <- mytrue(A1s=1,A2s=1,A3s=1) - mytrue(A1s=0,A2s=0,A3s=0)

    WMSM_withU <- weightitMSM(
      list(A1 ~ X1_1 + X2_1,
           A2 ~ A1 + X1_1 + X2_1+ X1_2 + X2_2,
           A3 ~ A2 + A1 + X1_1 + X2_1+ X1_2 + X2_2+ X1_3 + X2_3),
      method = "glm",
      data = full_dat,
      stabilize = T)
    full_dat$weights_withU <- WMSM_withU$weights

    MSM_withU <- lm(Y ~ (A1 + A2 + A3), weights = weights_withU, data = full_dat)

    p<-avg_predictions(MSM_withU,
                       vcov = "HC3",
                       newdata = datagrid(A1 = 0:1, A2 = 0:1, A3 = 0:1),
                       by = c("A1", "A2", "A3"),
                       wts = "weights_withU",
                       type = "response")
    results.it[1,2]<- hypotheses(p, "b8 - b1 = 0")$estimate
    results.it[1,3]<- hypotheses(p, "b8 - b1 = 0")$std.error
    #remove all the X1_
    WMSM_noU <- weightitMSM(
      list(A1 ~ X2_1,
           A2 ~ A1+X2_1+X2_2,
           A3 ~ A2+A1+X2_1+X2_2+X2_3),
      method = "glm",
      data = full_dat,
      stabilize = T)

    full_dat$weights_withoutU <- WMSM_noU$weights
    MSM_noU <- lm(Y ~ A1 + A2 + A3, weights = weights_withoutU, data = full_dat)
    p<-avg_predictions(MSM_noU,
                       vcov = "HC3",
                       newdata = datagrid(A1 = 0:1, A2 = 0:1, A3 = 0:1),
                       by = c("A1", "A2", "A3"),
                       wts = "weights_withoutU",
                       type = "response")

    results.it[1,4]<- hypotheses(p, "b8 - b1 = 0")$estimate
    results.it[1,5]<- hypotheses(p, "b8 - b1 = 0")$std.error

    #
    # 4. Marginal Sturtural Model with sensitivity function

    trtmodel_1 = glm(A1 ~ X2_1, family = "binomial", data = full_dat)
    trt_1 <- predict(trtmodel_1,  newdata = full_dat, type = "response")  #f(A1 = 1|X1)
    trtmodel_2 = glm(A2 ~ A1 + X2_1 + X2_2, family = "binomial", data = full_dat)
    trt_2 <- predict(trtmodel_2, newdata = full_dat %>% mutate(A1 = 1), type = "response")  #f(A2 = 1|X,A1 = 1)
    con_2 <- predict(trtmodel_2, newdata = full_dat %>% mutate(A1 = 0), type = "response")  #f(A2 = 1|X,A1 = 0)
    trtmodel_3 = glm(A3 ~ A2 + A1 + X2_1 + X2_2 + X2_3, family = "binomial", data = full_dat)
    trt_3 <- predict(trtmodel_3,  newdata = full_dat%>% mutate(A2 = 1, A1 = 1), type = "response") # f(A3 = 1|A2 = A1 =1, X)
    con_3_a1 <- predict(trtmodel_3,  newdata = full_dat%>% mutate(A2 = 0, A1 = 1), type = "response") # f(A3 = 1|A2 = 0, A1 =1 ,X)
    con_3_a2 <- predict(trtmodel_3,  newdata = full_dat%>% mutate(A2 = 1, A1 = 0), type = "response") # f(A3 = 1|A2 = 1, A1 =0 ,X)
    con_3 <- predict(trtmodel_3,  newdata = full_dat%>% mutate(A2 = 0, A1 = 0), type = "response") # f(A3 = 1|A1= A2 = 0,X)

    full_dat = full_dat %>% mutate(ptrt_1 = trt_1,
                                   ptrt_2 = trt_2,
                                   ptrt_3 = trt_3,
                                   pcon_2 = con_2,
                                   pcon_3_a1 = con_3_a1,
                                   pcon_3_a2 = con_3_a2,
                                   pcon_3 = con_3)

    #SF function;
    #Visit 3;
    #c(3,a1,a2,a3,x1,x2,x3) = E(APO|a3,a2,a1,x3,x2,x1) - E(APO|1-a3,a2,a1,x3,x2,x1)
    SF_calculation3 <- function(df = full_dat, A1=1, A2=1, A3=1){
      varName <- paste0("Y", A1, A2, A3)
      mx1<- glm(X2_1~1, family = "binomial", data= df)
      mx2<- glm(X2_2~X2_1+A1, family = "binomial", data= df)
      mx3<- glm(X2_3~X2_2+X2_1+A1+A2, family = "binomial", data= df)
      exp1 <- df[df$A1 == A1 & df$A2 == A2 & df$A3 == A3, varName]*predict(mx3,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == A3,], type="response")*predict(mx2,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == A3,], type="response")*predict(mx1,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == A3,], type="response")
      exp2 <- df[df$A1 == A1 & df$A2 == A2 & df$A3 == (1-A3), varName]*predict(mx3,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == (1-A3),], type="response")*predict(mx2,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == (1-A3),], type="response")*predict(mx1,newdata = df[df$A1 == A1 & df$A2 == A2 & df$A3 == (1-A3),], type="response")
      return(mean(exp1) - mean(exp2))
    }
    SF_calculation3(A1=1, A2=1, A3=1)
    # -1.108315
    # -0.4463642

    #Visit 2;
    #c(3,a1,a2,a3,x1,x2,x3) = E(APO|a2,a1,x2,x1) - E(APO|1-a2,a1,x2,x1)
    SF_calculation2 <- function(df = full_dat, A1=1, A2=1, A3=1){
      varName <- paste0("Y", A1, A2, A3)
      mx1<- glm(X2_1~1, family = "binomial", data= df)
      mx2<- glm(X2_2~X2_1+A1, family = "binomial", data= df)
      mx3<- glm(X2_3~X2_2+X2_1+A1+A2, family = "binomial", data= df)
      exp1 <- df[df$A1 == A1 & df$A2 == A2, varName]*predict(mx2,newdata = df[df$A1 == A1 & df$A2 == A2,], type="response")*predict(mx1,newdata = df[df$A1 == A1 & df$A2 == A2,], type="response")
      exp2 <- df[df$A1 == A1 & df$A2 == (1-A2), varName]*predict(mx2,newdata = df[df$A1 == A1 & df$A2 == (1-A2),], type="response")*predict(mx1,newdata = df[df$A1 == A1 & df$A2 == (1-A2),], type="response")
      return(mean(exp1) - mean(exp2))
    }
    SF_calculation2(A1=1, A2=1, A3=1)
    # -0.6633976
    # -0.4544505
    #Visit 1;
    #c(3,a1,a2,a3,x1,x2,x3) = E(APO|a1,x1) - E(APO|1-a1,x1)
    SF_calculation1 <- function(df = full_dat, A1=1, A2=1, A3=1){
      varName <- paste0("Y", A1, A2, A3)
      mx1<- glm(X2_1~1, family = "binomial", data= df)
      mx2<- glm(X2_2~X2_1+A1, family = "binomial", data= df)
      mx3<- glm(X2_3~X2_2+X2_1+A1+A2, family = "binomial", data= df)
      exp1 <- df[df$A1 == A1  , varName]*predict(mx1,newdata = df[df$A1 == 1,], type="response")
      exp2 <- df[df$A1 == (1-A1), varName]*predict(mx1,newdata = df[df$A1 == (1-A1),], type="response")
      return(mean(exp1) - mean(exp2))
    }
    # SF_calculation1(A1=1, A2=1, A3=1)


    full_dat = full_dat  %>% mutate(bias_1 = case_when(
      A3 == 1& A2 == 1& A1 == 1 ~ SF_calculation1(A1=1, A2=1, A3=1)*(1-ptrt_1) ,
      A3 == 1& A2 == 1& A1 == 0 ~ SF_calculation1(A1=1, A2=1, A3=0)*ptrt_1,
      A3 == 1& A2 == 0& A1 == 0 ~ SF_calculation1(A1=1, A2=0, A3=0)*ptrt_1 ,
      A3 == 1& A2 == 0& A1 == 1 ~ SF_calculation1(A1=1, A2=0, A3=1)*(1-ptrt_1),
      A3 == 0& A2 == 1& A1 == 1 ~ SF_calculation1(A1=0, A2=1, A3=1)*(1-ptrt_1),
      A3 == 0& A2 == 1& A1 == 0 ~ SF_calculation1(A1=0, A2=1, A3=0)*ptrt_1 ,
      A3 == 0& A2 == 0& A1 == 1 ~ SF_calculation1(A1=0, A2=0, A3=1)*(1-ptrt_1) ,
      A3 == 0& A2 == 0& A1 == 0 ~ SF_calculation1(A1=0, A2=0, A3=0)*ptrt_1,
      TRUE ~ NA_real_)) %>% mutate(bias_2 = case_when(
        A3 == 1& A2 == 1& A1 == 1 ~ SF_calculation2(A1=1, A2=1, A3=1)*(1-ptrt_2) ,
        A3 == 1& A2 == 1& A1 == 0 ~ SF_calculation2(A1=1, A2=1, A3=0)*(1-pcon_2),
        A3 == 1& A2 == 0& A1 == 0 ~ SF_calculation2(A1=1, A2=0, A3=0)*pcon_2,
        A3 == 1& A2 == 0& A1 == 1 ~ SF_calculation2(A1=1, A2=0, A3=1)*ptrt_2,
        A3 == 0& A2 == 1& A1 == 1 ~ SF_calculation2(A1=0, A2=1, A3=1)*(1-ptrt_2),
        A3 == 0& A2 == 1& A1 == 0 ~ SF_calculation2(A1=0, A2=1, A3=0)*(1-pcon_2),
        A3 == 0& A2 == 0& A1 == 1 ~ SF_calculation2(A1=0, A2=0, A3=1)*ptrt_2,
        A3 == 0& A2 == 0& A1 == 0 ~ SF_calculation2(A1=0, A2=0, A3=0)*pcon_2,
        TRUE ~ NA_real_)) %>% mutate(bias_3 = case_when(
          A3 == 1& A2 == 1& A1 == 1 ~ SF_calculation3(A1=1, A2=1, A3=1)*(1-ptrt_3) ,
          A3 == 1& A2 == 1& A1 == 0 ~ SF_calculation3(A1=1, A2=1, A3=0)*(1 - pcon_3_a2),
          A3 == 1& A2 == 0& A1 == 0 ~ SF_calculation3(A1=1, A2=0, A3=0)*(1- pcon_3),
          A3 == 1& A2 == 0& A1 == 1 ~ SF_calculation3(A1=1, A2=0, A3=1)*(1 - pcon_3_a1),
          A3 == 0& A2 == 1& A1 == 1 ~ SF_calculation3(A1=0, A2=1, A3=1)*ptrt_3,
          A3 == 0& A2 == 1& A1 == 0 ~ SF_calculation3(A1=0, A2=1, A3=0)*pcon_3_a2,
          A3 == 0& A2 == 0& A1 == 1 ~ SF_calculation3(A1=0, A2=0, A3=1)*pcon_3_a1,
          A3 == 0& A2 == 0& A1 == 0 ~ SF_calculation3(A1=0, A2=0, A3=0)*pcon_3,
          TRUE ~ NA_real_))

    full_dat = full_dat %>% mutate(Y_sf = Y - bias_1 - bias_2 - bias_3)

    MSM_noU_sf <- lm(Y_sf ~ A1 + A2 + A3, weights = weights_withoutU, data = full_dat)
    p<-avg_predictions(MSM_noU_sf,
                       vcov = "HC3",
                       newdata = datagrid(A1 = 0:1, A2 = 0:1, A3 = 0:1),
                       by = c("A1", "A2", "A3"),
                       wts = "weights_withoutU",
                       type = "response")

    results.it[1,6]<- hypotheses(p, "b8 - b1 = 0")$estimate
    results.it[1,7]<- hypotheses(p, "b8 - b1 = 0")$std.error


    ### Variance Calculation For MSM
    # bootstrap variance;
    Bootest_withU<-rep(NA, nboot)

    Bootest_noU<-rep(NA, nboot)

    Bootest_SF<-rep(NA, nboot)
    for (j1 in 1:nboot){
      bootidx <- sample(1:dim(full_dat)[1], replace=TRUE)
      bootdata <- full_dat[bootidx,]

      MSM_withU <- lm(Y ~ (A1 + A2 + A3), weights = weights_withU, data = bootdata)
      Bootest_withU[j1] <- as.numeric(MSM_withU$coefficients[2]+ MSM_withU$coefficients[3] + MSM_withU$coefficients[4])
      MSM_withoutU <- lm(Y ~ (A1 + A2 + A3), weights = weights_withoutU, data = bootdata)
      Bootest_noU[j1] <- as.numeric(MSM_withoutU$coefficients[2]+ MSM_withoutU$coefficients[3] + MSM_withoutU$coefficients[4])
      MSM_withoutUSF <- lm(Y_sf ~ (A1 + A2 + A3), weights = weights_withoutU, data = bootdata)
      Bootest_SF[j1] <- as.numeric(MSM_withoutUSF$coefficients[2]+ MSM_withoutUSF$coefficients[3] + MSM_withoutUSF$coefficients[4])
    }
    #
    #
    #
    results.it[1,3]<- sd(Bootest_withU)
    results.it[1,5]<- sd(Bootest_noU)
    results.it[1,7]<- sd(Bootest_SF)

    #    # combining parallel results;
    cbind(i,results.it)
    #
  }
  save(results_run,file=paste0(outfile,".RData"))
  write.table(results_run, file = paste0(outfile,".txt"), row.names = FALSE,col.names = TRUE)

}

mysim.bin.cens("500binaryU", from=1, to= 1000, U=1, #parameters defining the simulation iterations and outcome;
               samplesize=500,#parameters for data simulation
               B=500, pdraw = 8000, nburnin = 6000,#parameters for analysis;
               Px11 = .5, Px12 = .4, coef_u = -0.2, coef_a = .1)#simulate measured covariates

