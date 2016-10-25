estimate_params<-function(data){
  # estimates parameters from data and returns them

  # generate likelihood func to be optimised
  target<-targetGen(data)
  # set arguments to optim function
  control<-list(fnscale=-1)
  uppers<-rep(5,data$s)
  lowers<-rep(-5,data$s)
  guess<-pickStart(target,uppers,lowers,10,data$s)


  # find max likelihood estimate for scale and alpha params
  answer<-optim(guess,target,control=control, hessian=T)

  # get omegas and alpha
  delta_hat<-giveDeltas(answer$par[1:(data$p+data$q)])

  if (data$det==F){
    alpha_hat<- boot::inv.logit(answer$par[data$s])
  }
  else{
    #set alpha to 1 for deterministic
    alpha_hat<-1
  }

  A<-createAmat(data$d1, giveOmegas(delta_hat), alpha=alpha_hat)

  # get chol of covar matrix A and estimate parameters
  R<-chol(A)
  beta_hat<-estimBeta(data$hmat,data$fofD,R)
  V_hat<-estimV(beta_hat,data$hmat,data$fofD,data$N,data$r,R)
  sigma2_hat<- alpha_hat*V_hat
  nugget<-V_hat-sigma2_hat

  #return as a list
  return(list(beta_hat=beta_hat,sigma2_hat=sigma2_hat,
              nugget_hat=nugget, alpha_hat=alpha_hat,
              v_hat=V_hat, delta_hat=delta_hat,
              omega_hat=giveOmegas(delta_hat), hessian=answer$hessian,
              R=R, A=A))
}




estimate_params2<-function(data, reps=10){
  # estimates parameters from data and returns them
  # picks best of 'reps' attempts

  # generate likelihood func to be optimised
  target<-targetGen(data)

  #create a vector of length reps, of random optimisation starting points
  x<-rep(data$s,reps)
  guess<-t(mapply(runif,x,-5,1))

  # set arguments to optim function
  if (data$det==F){
    lower<-rep(-6,data$s-1)
    upper<-rep(1,data$s-1)
    lower<-c(lower,-Inf)
    upper<-c(upper, Inf)
  } else {
    lower <- rep(-6,data$s)
    upper <- rep(1, data$s)
  }
  control=list(fnscale=-1)
  if (data$s == 1){
    guess <- t(guess) # annoying behaviour of mapply.
    method <- "Brent" #better in 1d?
  } else {
    method <- "L-BFGS-B"
  }

  # run optimisation routine from each of these points
  test<- apply(guess, 1, optim, fn=target, control=control,
               hessian=T, method=method, lower=lower, upper= upper)

  # for each result extract the parameters,
  # log-like value and an indicator as to whether the hessian was -ve def
  outs<-sapply(test, "[[", "par")
  vals<-sapply(test, "[[", "value")
  convergence <- sapply(test, "[[", "convergence")
  message <- sapply(test, "[[", "message")

  hess_neg_def<-sapply(test,
                       function(thing){all(eigen(thing$hessian)$values<0)})
  if (data$s>1){
    outs<-t(outs)
  }

  results<-cbind(vals,outs, hess_neg_def, convergence)

  results<-results[order(results[,1], decreasing=T),]
  answer<-NULL
  for (i in 1:dim(results)[1]){
    if (results[i,"hess_neg_def"] & sum(results[i,2:data$s])<1){
      answer<-results[i,]
      break
    }
  }
  if (is.null(answer)){
    answer<-results[1,]
    warning("no reasonable mode-likelihood estimates of the correlation
            hyperparameters found - best estimate used ")
  }

  # find max likelihood estimate for scale and alpha params
  #answer<-optim(guess,target,control=control, hessian=T)

  # get omegas and alpha
  delta_hat<-giveDeltas(answer[2:(data$p+data$q+1)])

  if (data$det==F){
    alpha_hat<- boot::inv.logit(answer[data$s+1])
  }
  else{
    #set alpha to 1 for deterministic
    alpha_hat=1
  }

  A<-createAmat(data$d1,giveOmegas(delta_hat),alpha=alpha_hat)

  # get chol of covar matrix A and estimate parameters
  R<-chol(A)
  beta_hat<-estimBeta(data$hmat,data$fofD,R)
  V_hat<-estimV(beta_hat,data$hmat,data$fofD,data$N,data$r,R)
  sigma2_hat<- alpha_hat*V_hat
  nugget<-V_hat-sigma2_hat

  #return as a list
  return(list(beta_hat=beta_hat,sigma2_hat=sigma2_hat,
              nugget_hat=nugget, alpha_hat=alpha_hat,
              v_hat=V_hat, delta_hat=delta_hat,
              omega_hat=giveOmegas(delta_hat),
              R=R,A=A, modesFound=results))
}

targetGen<-function(data){
  #spits out function to pass to optim
  target<-function(params){
    # transform from form sent to optim function

    if (data$det==F){
      alpha<-boot::inv.logit(params[data$s])
      omegas<- giveOmegas(deltas=giveDeltas(params[1:(data$s-1)]))
      #create covar matrix
      amat<-createAmat(data$d1,omegas,alpha)
    }
    else {
      omegas<- giveOmegas(deltas=giveDeltas(params))

      #create covar matrix
      amat<-get_corr_matrix(data$d1,omegas=omegas[1:(data$s)])
    }
    #if possible, get chol. decomposition
    R<-tryCatch(chol(amat),error=function(e) return(-999999))
    if (is.matrix(R)==F){ return (R)}

    #estimate parameters
    outs<-estimBeta(data$hmat,data$fofD,R, giveK=T)
    vhat<-estimV(outs$beta_hat,data$hmat,data$fofD,data$N,data$r,R)

    #calc likelihood
    term1<- -((data$N-data$r)/2) * log(vhat)
    term2<- - sum(log(diag(R)))
    term3<- -sum(log(diag(outs$K)))

    return(term1+term2 +term3)
  }
  return(target)
}

estimCovarMatrix<-function(data,newdata,estimates=NULL, hfunc=createH){
  if (is.null(estimates)){
    # no estimates provided - use 'real' parameter values (for testing)
    A <- createAmat(data$d1,omegas = data$reality$omegas,alpha=data$reality$alpha)
    R <- chol(A) # note errors with data$reality$R
    sigma2<-data$reality$sigma2
    omegas<-data$reality$omegas
    nugget<-data$reality$nugget
    alpha <- data$reality$alpha
    v <- data$reality$v
  }
  else{
    #A<- get_corr_matrix(data$d1,omegas=estimates$omega_hat)
    #R<- chol(A)
    R <- estimates$R
    sigma2<-estimates$sigma2_hat
    omegas<-estimates$omega_hat
    nugget<-estimates$nugget_hat
    alpha <- estimates$alpha
    v <-estimates$v
  }
  hmatdash<- hfunc(newdata)
  newdata<-as.matrix(newdata)

  cmat <- createAmat(newdata,omegas,alpha=alpha)

  tmat <- alpha *  get_t_matrix(data$d1,newdata,omegas=omegas)

  w<-backsolve(R, t(tmat), transpose=T)
  covar<-cmat - t(w) %*% w
  I<-hmatdash - tmat %*%
    backsolve( R, backsolve(R,  data$hmat, transpose=T))
  w<-backsolve(R, data$hmat, transpose = TRUE)
  Q<-t(w)%*%w
  K<-chol(Q)
  part3<- I %*% backsolve(K, backsolve(K, t(I), transpose = TRUE))
  covar<-v[1] * (covar +part3)
  #   if (length(covar)==1){
  #     covar<-covar + estimates$nugget
  #   }
  #   else{
  #     covar<- covar + diag(rep(nugget, dim(covar)[1]))
  #   }
  return(covar)
}

estimBeta<-function(hmat,fofD,R,giveK=F){
  #estimate betas from data, response
  # and cholesky decomposition of the correl matrix
  # K optionally return if needed
  w<-backsolve(R, hmat, transpose = TRUE)
  Q<-t(w)%*%w
  K<-chol(Q)
  part1<-backsolve(K, backsolve(K, t(hmat), transpose = TRUE))
  part2<-backsolve(R, backsolve(R, fofD, transpose = TRUE))
  beta_hat <- part1%*%part2
  if (giveK) return(list(beta_hat=beta_hat,K=K))
  else return (beta_hat)
}

estimV<-function(betahat, hmat, fofD,n,r,R){
  # estimate total  variance given beta hat and data
  e<-(fofD - hmat%*%betahat)
  w2<-backsolve(R, e, transpose = TRUE)
  eae<- t(w2)%*%w2
  return( (1/(n-r-2)) * eae )
}


pickStart<-function(target,uppers, lowers, g, r){
  # choose starting point for optim
  # pick g random starting points of dimension r
  guess<-c()
  for (i in 1:r){
    guess<-cbind(guess,runif(g,min=lowers[i],max=uppers[i]))
  }
  # calculate their scores
  score<-rep(0,g)
  for (i in 1:dim(guess)[1]){
    score[i]<-target(guess[i,])
  }
  # return the max
  start_index<-which.max(score)
  return(guess[start_index,])
}

estimateMean<-function(data,newPoints,estimates=NULL, hfunc=createH){
  #estimate mean surface

  if (is.null(estimates)){
    A <- createAmat(data$d1,omegas = data$reality$omegas,alpha=data$reality$alpha)
    R <- chol(A) # note errors with data$reality$R
    e<-data$reality$e
    betas<-data$reality$betas
    sigma2<-data$reality$sigma2
    omegas<-data$reality$omegas
    alpha<-data$reality$alpha
  }
  else{
    R<-estimates$R
    e <- data$fofD - data$hmat%*%estimates$beta_hat
    betas <- estimates$beta_hat
    sigma2 <-estimates$sigma2_hat
    omegas <-estimates$omega_hat
    alpha <- estimates$alpha
  }
  hmat <- hfunc(newPoints)
  tmat <- alpha * get_t_matrix(data$d1,newPoints,omegas=omegas)
  part2 <- backsolve(R, backsolve(R, e, transpose = TRUE))
  mean <- hmat %*% betas + ((tmat) %*% part2)
  return(mean)
}
