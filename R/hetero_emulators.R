library(lhs)



targetGen_ML<- function(data){
  
  target<- function(params){
    sigma <- exp(params[1])
    nug   <- exp(params[2])
    omega <- giveOmegas(giveDeltas(params[3:length(params)]))
    K <- sigma * get_corr_matrix(data$d1,omega) 
    R<-tryCatch(chol(K + diag(nug, data$N)),error=function(e) return(-999999))
    if (is.matrix(R)==F){ return (R)}
    # add prior for omegas? and sigma2? 
    LL<-  ( - 0.5 * t(data$fofD) %*% backsolve(R ,backsolve(R, data$fofD, transpose=T))
            - 0.5 * sum(log(diag(R))) ) 
    return (LL )
  }
  return(target)
}

estimate_params_ML<-function(data){
  target<-targetGen_ML(data)
  control<-list(fnscale=-1)
  # sigma, nugget, omegas
  guess<-c(log(1),log(1), rep(0,data$p+data$q))
  answer<-optim(guess, target, control=control)
  params<-answer$par
  sigma <- exp(params[1])
  nug   <- exp(params[2])
  omega <- giveOmegas(giveDeltas(params[3:length(params)]))
  deltas <- giveDeltas(params[3:length(params)])
  K <- sigma * get_corr_matrix(data$d1,omega) 
  R<-tryCatch(chol(K + diag(nug, data$N)),error=function(e) return(-999999))
  return(list(sigma2_hat=sigma, 
              nugget_hat=nug, alpha_hat=sigma/(sigma+nug), 
              v_hat=sigma+nug, delta_hat=deltas, 
              omega_hat=omega, R=R, A=K))
}

#estims<-estimate_params_ML(data)

targetGen_ML_fixed_noise<- function(data){
  
  target<- function(params){
    sigma <- exp(params[1])
    omega <- giveOmegas(giveDeltas(params[2:length(params)]))
    K <- sigma * get_corr_matrix(data$d1,omega) 
    R<-tryCatch(chol(K + diag(as.vector(data$noise), data$N)),error=function(e) return(-999999))
    if (is.matrix(R)==F){ return (R)}
    # add prior for omegas? and sigma2? 
    #taus<- params[2:length(params)]
    #tau_lo <- -10.6; tau_hi <- 9.2
    #tau_prior<- -2* sum(exp(-2*(taus- tau_lo )) + exp(2*(taus-tau_hi)))
    LL<-  ( - 0.5 * t(data$fofD) %*% backsolve(R ,backsolve(R, data$fofD, transpose=T))
            - 0.5 * sum(log(diag(R))) ) 
    return (LL)# + tau_prior )
  }
  return(target)
}

estimate_params_ML_fixed_noise<-function(data,reps){
  target<-targetGen_ML_fixed_noise(data)
  control<-list(fnscale=-1)
  # sigma, nugget, omegas
  
  x<-rep(dim(data$d1)[2] +1,reps)
  guess<-t(mapply(runif,x,-6,5))
  
  # run optimisation routine from each of these points
  lower<-rep(-6,dim(data$d1)[2])
  upper<-rep(10,dim(data$d1)[2])
  lower<-c(-Inf, lower)
  upper<-c(Inf, upper)
  
  
  test<- apply(guess,1, optim, fn=target, control=control,
               hessian=T,method="L-BFGS-B", lower=lower, upper= upper)
  outs<-sapply(test, "[[", "par")
  vals<-sapply(test, "[[", "value")
  hess_neg_def<-sapply(test, 
                       function(optim_out){all(eigen(optim_out$hessian)$values<0)})
  
  results<-cbind(vals,t(outs), hess_neg_def)
  
  results<-results[order(results[,1], decreasing=T),]
  
  answer<-results[1,]
  
  
  params<-answer[2:(length(answer)-1)]
  sigma <- exp(params[1])
  
  omega <- giveOmegas(giveDeltas(params[2:length(params)]))
  deltas <- giveDeltas(params[2:length(params)])
  K <- sigma * get_corr_matrix(data$d1,omega) 
  R<-tryCatch(chol(K + diag(as.vector(data$noise), data$N)),error=function(e) return(-999999))
  return(list(sigma2_hat=sigma, 
              nugget_hat=data$noise, alpha_hat=sigma/(sigma+data$noise), 
              v_hat=sigma+data$noise, delta_hat=deltas, 
              omega_hat=omega, R=R, A=K))
}


estim_point<-function(data,estims, newdata){
  newdata <-as.matrix(newdata)
  tmat<-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  return(tmat %*%  backsolve(estims$R ,backsolve(estims$R, data$fofD, transpose=T)))
}

estim_var_var <- function(data,estims,newdata){
  newdata <-as.matrix(newdata)
  tmat <-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  cmat <- estims$sigma2_hat*get_corr_matrix(newdata,estims$omega_hat) #+ 
  #    diag(data$noise, dim(newdata)[1] )
  return (  cmat - (tmat) %*%  
              backsolve(estims$R ,backsolve(estims$R, t(tmat), transpose=T)))
}


estim_mean_var <- function(data,estims,newdata,noise){
  newdata <-as.matrix(newdata)
  tmat <-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  cmat <- estims$sigma2_hat*get_corr_matrix(newdata,estims$omega_hat) + 
    diag(as.vector(noise), dim(newdata)[1] )
  return (  cmat - (tmat) %*%  
              backsolve(estims$R ,backsolve(estims$R, t(tmat), transpose=T)))
}

create_var_data_object<-function(input,output,obs_per_point){
  # input must be UNTRANSFORMED VARIANCE
  # function to create data object for input into estimate pars.
  # from an lhs and a vector of outputs. 
  # sample must be [n,p+q]
  # output must be [n,1]
  d1 = as.matrix(input)
  hmat = as.matrix(cbind(rep(1,length(output)),d1))
  n <- obs_per_point
  r = dim(d1)[2] + 1  # number of basis functions
  p = dim(d1)[2]      # number of input variables
  q = 0                  # only used in calibration case
  N = dim(d1)[1]      # number of observations
  s = dim(d1)[2]  +1    # number of correlation hyperparameters
  fofD = as.matrix(log(output))- digamma((n-1)/2) - log(2) + log(n-1)
  noise = trigamma((n-1)/2)
  data <- list(d1 = d1, hmat = hmat, fofD = fofD, # data objects
               r = r, p = p, q = q, s = s, N = N, # size objects
               obs_per_point = obs_per_point, # number of observations used to calc each log sample var
               noise = noise # uncertainty due to log sample estimate of variance
  ) 
  return (data)
}

create_mean_data_object<-function(input,output,log_var_preds, obs_per_point){
  # function to create data object for input into estimate pars.
  # from an lhs and a vector of outputs. 
  # sample must be [N,p+q]
  # output must be [N,1] 
  
  d1 = as.matrix(input)
  hmat = as.matrix(cbind(rep(1,length(output)),d1))
  fofD = as.matrix(output)
  r = dim(d1)[2] + 1  # number of basis functions
  p = dim(d1)[2]      # number of input variables
  q = 0                  # only used in calibration case
  N = dim(d1)[1]      # number of observations
  s = dim(d1)[2]  +1    # number of correlation hyperparameters
  noise = exp(log_var_preds)/obs_per_point
  data <- list(d1 = d1, hmat = hmat, fofD = fofD, # data objects
               r = r, p = p, q = q, s = s, N = N, # size objects
               obs_per_point = obs_per_point, # number of observations used to calc each log sample var
               noise = noise # uncertainty in mean estimate due to sampling variance
  ) 
  return (data)
}


get_func_var_no_mean<-function(data,estims,data_var,estims_var,newdata){
  pred_noise<-exp(estim_point(data_var,estims_var,newdata))
  newdata <-as.matrix(newdata)
  tmat <-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  cmat <- estims$sigma2_hat*get_corr_matrix(newdata,estims$omega_hat) + 
    diag(as.vector(pred_noise))
  return (  cmat - (tmat) %*%  
              backsolve(estims$R ,backsolve(estims$R, t(tmat), transpose=T)))
}


gen_target_fixed_noise<-function(data){
  target<-function(params){
    sigma2 <- exp(params[1])
    omega <- giveOmegas(giveDeltas(params[2:length(params)]))
    A <- get_corr_matrix(data$d1, omega)
    ## calc LL using Ks -----------------------------
    K <- sigma2 *A + diag(as.vector(data$noise), data$N)
    R<-tryCatch(chol(K),error=function(e) return(-999999))
    if (is.matrix(R)==F){ return (R)}
    HKH <- t(data$hmat) %*% backsolve(R, backsolve(R, data$hmat, transpose=T))
    RR <- chol(HKH)
    HKfofD <- t(data$hmat) %*% backsolve(R, backsolve(R, data$fofD, transpose=T))
    beta_hat <- backsolve(RR, backsolve( RR, t(data$hmat), transpose=T)) %*% backsolve(R, backsolve(R, data$fofD, transpose=T))
    e <- data$fofD - data$hmat %*% beta_hat
    Q <- backsolve( R, e, transpose=T)
    eKe <- t(Q) %*% (Q)
    P <-  backsolve(RR, backsolve( RR, t(data$hmat), transpose=T))
    RRP <- backsolve(RR, P, transpose=T)
    C <-  t(RRP) %*% RRP
    eCe <- t(e) %*% C %*% e 
    LL <- ( - sum(log(diag(R))) -  sum(log(diag(RR))) # determinants
            - 0.5 * eKe # + 0.5 * eCe  # error terms
            - log(sigma2) ) # prior
    return(LL)
  }
  return(target)
}



estimate_par<-function(data,reps=10){
  target  <- gen_target(data)
  control <- list(fnscale=-1)
  
  x<-rep(dim(data$d1)[2] +2,reps)
  guess<-t(mapply(runif,x,-6,5))
  
  # run optimisation routine from each of these points
  lower<-rep(-6,dim(data$d1)[2])
  upper<-rep(5,dim(data$d1)[2])
  if (data$det==F){
    lower<-c(-Inf, -Inf, lower)
    upper<-c(Inf, Inf, upper)
  }
  control=list(fnscale=-1)
  test<- apply(guess,1, optim, fn=target, control=control,
               hessian=T,method="L-BFGS-B", lower=lower, upper= upper)
  
  # for each result extract the parameters, 
  # log-like value and an indicator as to whether the hessian was -ve def
  outs<-sapply(test, "[[", "par")
  vals<-sapply(test, "[[", "value")
  hess_neg_def<-sapply(test, 
                       function(optim_out){all(eigen(optim_out$hessian)$values<0)})
  
  results<-cbind(vals,t(outs), hess_neg_def)
  
  results<-results[order(results[,1], decreasing=T),]
  
  answer<-results[1,]
  
  # get omegas and alpha 
  delta_hat<-giveDeltas(answer[4:(data$p+data$q+3)])
  omegas_hat<- giveOmegas(delta_hat)
  
  sigma2_hat<-exp(answer[2])
  nugget<-exp(answer[3])
  
  A <- get_corr_matrix(data$d1,  omegas_hat)
  K <- sigma2_hat*A + diag(nugget, data$N)
  # get chol of covar matrix A and estimate parameters
  R<-chol(K)
  
  beta_hat<-estimBeta(data$hmat,data$fofD,R)
  
  #return as a list
  return(list(beta_hat=beta_hat,sigma2_hat=sigma2_hat, 
              nugget_hat=nugget, omega_hat=giveOmegas(delta_hat),
              R=R,A=A, modesFound=results))
}

estimate_par_fixed_noise<-function(data,reps=10){
  target  <- gen_target_fixed_noise(data)
  control <- list(fnscale=-1)
  
  x<-rep(dim(data$d1)[2] +1,reps)
  guess<-t(mapply(runif,x,-6,5))
  
  # run optimisation routine from each of these points
  lower<-rep(-6,dim(data$d1)[2])
  upper<-rep(5,dim(data$d1)[2])
  
  lower<-c(-Inf, lower)
  upper<-c(Inf, upper)
  
  control=list(fnscale=-1)
  test<- apply(guess,1, optim, fn=target, control=control,
               hessian=T,method="L-BFGS-B", lower=lower, upper= upper)
  
  # for each result extract the parameters, 
  # log-like value and an indicator as to whether the hessian was -ve def
  outs<-sapply(test, "[[", "par")
  vals<-sapply(test, "[[", "value")
  hess_neg_def<-sapply(test, 
                       function(optim_out){all(eigen(optim_out$hessian)$values<0)})
  
  results<-cbind(vals,t(outs), hess_neg_def)
  
  results<-results[order(results[,1], decreasing=T),]
  
  answer<-results[1,]
  
  # get omegas and alpha 
  delta_hat  <- giveDeltas(answer[3:(data$p+data$q+2)])
  omegas_hat <- giveOmegas(delta_hat)
  
  sigma2_hat<-exp(answer[2])
  nugget<-data$noise
  
  A <- get_corr_matrix(data$d1,  omegas_hat)
  K <- sigma2_hat*A + diag(as.vector(nugget), data$N)
  # get chol of covar matrix A and estimate parameters
  R<-chol(K)
  
  beta_hat<-estimBeta(data$hmat,data$fofD,R)
  
  #return as a list
  return(list(beta_hat=beta_hat,sigma2_hat=sigma2_hat, 
              nugget_hat=nugget, omega_hat=giveOmegas(delta_hat),
              R=R,A=A, modesFound=results))
}



estim_non_zero_mean_point<-function(data,estims, newdata, hfunc=createH){
  newdata <-as.matrix(newdata)
  part1 <- hfunc(newdata) %*% estims$beta_hat
  tmat<-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  part2 <- tmat %*%  backsolve(estims$R ,backsolve(estims$R, data$fofD- data$hmat %*%estims$beta_hat, transpose=T))
  return(part1 + part2)
}

get_func_var<-function(data,estims,data_var,estims_var,newdata){
  pred_noise<-exp(estim_point(data_var,estims_var,newdata))
  newdata <-as.matrix(newdata)
  tmat <-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  cmat <- estims$sigma2_hat*get_corr_matrix(newdata,estims$omega_hat) + 
    diag(as.vector(pred_noise))
  RH<-createH(newdata) - tmat %*% backsolve(estims$R, backsolve(estims$R,data$hmat, transpose=T))
  w <- backsolve( estims$R,  data$hmat, transpose=T)
  Q <- t(w) %*% w
  P <- chol(Q)
  K <- backsolve(P, t(RH), transpose=T)
  # function is ok - maybe more stable implmentations available. 
  return (  cmat - (tmat) %*%  
              backsolve(estims$R ,backsolve(estims$R, t(tmat), transpose=T))
            + t(K) %*% K 
  )
}

get_func_var_in_mean <- function(data,estims,newdata){
  # same as above, just excluding additional variance from variance emulator.
  # note that heteroskedastic does allow better mean prediction.
  newdata <-as.matrix(newdata)
  tmat <-estims$sigma2_hat*get_t_matrix(data$d1,newdata,estims$omega_hat)  
  cmat <- estims$sigma2_hat*get_corr_matrix(newdata,estims$omega_hat) 
  RH<-createH(newdata) - tmat %*% backsolve(estims$R, backsolve(estims$R,data$hmat, transpose=T))
  w <- backsolve( estims$R,  data$hmat, transpose=T)
  Q <- t(w) %*% w
  P <- chol(Q)
  K <- backsolve(P, t(RH), transpose=T)
  
  return (  cmat - (tmat) %*%  
              backsolve(estims$R ,backsolve(estims$R, t(tmat), transpose=T))
            + t(K) %*% K 
  )
}

fit_hetero_emulator <- function(design, mean, variance, reps){
  var_em <- fit_variance_emulator(design, variance, reps)
  mean_em <- fit_mean_emulator(design, mean, reps, var_em)
  
  outob <- list(mean_em=mean_em, var_em=var_em)
  class(outob) <- "het_em"
  return(outob)
}

predict.het_em <- function(object, new_data=NULL, variance=c("func","mean"),
                           percentile=F){
  variance <- match.arg(variance)
  het_em <- object
  if (!is.null(new_data)){
    # TODO check new data size
    d1 <- as.matrix(new_data)
  }  
  else {
    d1 <- het_em$mean_em$data$d1
  }
  
  if (variance == "func"){
    pred_variance <- diag(get_func_var(het_em$mean_em$data,
                                       het_em$mean_em$estims,
                                       het_em$var_em$data,
                                       het_em$var_em$estims, 
                                       d1))
  }
  else {
    pred_variance <- diag(get_func_var_in_mean(het_em$mean_em$data,
                                               het_em$mean_em$estims,
                                               d1))
  }
  pred_mean <- estim_non_zero_mean_point(het_em$mean_em$data,
                                         het_em$mean_em$estims,
                                         d1)
  if (percentile){
    low <- (pred_mean - sqrt(pred_variance) * qnorm(1 - (1 - percentile)/2))
    high <- (pred_mean + sqrt(pred_variance) * qnorm(1 - (1 - percentile)/2))
    outlist <- list(mean = pred_mean, variance = pred_variance, 
                    var_type=variance, low,high)
    names(outlist)[4:5] <- c(paste0("q", 1 - (1 - percentile)/2),
                             paste0("q", (1 - percentile)/2))
    return(outlist)
  }
  return(list(mean = pred_mean, variance = pred_variance, var_type=variance))
}

fit_variance_emulator <- function(design, variance, reps){
  data_var <- create_var_data_object(design, variance, obs_per_point = reps) 
  estims_var <- estimate_params_ML_fixed_noise(data_var, 5)
  # predict at points
  logvars <- estim_point(data_var, estims_var, data_var$d1)
  var_covar     <- estim_var_var(data_var, estims_var, data_var$d1)
  var_variance <- diag(var_covar)
  return(list(data=data_var, estims=estims_var, preds=logvars, vars=var_variance))
}

fit_mean_emulator <- function(design, mean, reps, var_em){
  data_mean <- create_mean_data_object(design, mean,
                                       estim_point(var_em$data, var_em$estims,
                                                   var_em$data$d1),
                                       var_em$data$obs_per_point)
  estims_mean <- estimate_par_fixed_noise(data_mean)
  return(list(data=data_mean, estims=estims_mean))
}

plot.het_em <- function(object, percentile=0.95, variance="func"){
  het_em <- object
  preds <- predict(het_em, percentile = percentile,variance = variance)
  plot(preds$mean[order(preds$mean)], type="l")
  points(preds$q0.025[order(preds$mean)], type="l", lty=3)
  points(preds$q0.975[order(preds$mean)], type="l", lty=3)
}


