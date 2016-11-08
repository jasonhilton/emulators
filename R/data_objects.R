# functions to fit an GP emulator - as per the MUCM framework.
#    wwww.mucm.ac.uk
# author Jason Hilton : jason_hilton@yahoo.com
# Many thanks to BACCO package author Robin K. S. Hankin - particularly for the
# fast corr.matrix implementation

create_data_object<-function(input,output,det,hfunc=createH){
  # function to create data object for input into estimate pars.
  # from an lhs and a vector of outputs.
  # sample must be [n,p+q]
  # output must be [n,1]
  # det must be boolean indicating whether the data if deterministic (True)
  # or stocastic (False)
  d1 = as.matrix(input)
  #hmat = as.matrix(cbind(rep(1,length(output)),d1))
  hmat = hfunc(input)
  fofD = as.matrix(output)
  #r = dim(d1)[2] + 1  # number of basis functions
  r = dim(hmat)[2]    # number of basis functions
  p = dim(d1)[2]      # number of input variables
  q = 0                  # only used in calibration case
  N = dim(d1)[1]      # number of observations
  s = dim(d1)[2]      # number of correlation hyperparameters
  if (det == F){
    s = s + 1            # extra nugget hyperparameter in stochastic case
  }
  outob <- list(d1 = d1, hmat = hmat, fofD = fofD, # data objects
               r = r, p = p, q = q, s = s, N = N, # size objects
               det = det # deterministic flag
  )
  class(outob) <- "em_data"
  return (outob)
}


createH<- function(data){
  #create the Hmat with a constant intercept column
  return(as.matrix(cbind(rep(1, dim(data)[1]), data)))
}

createHzero <- function(data){
  return (as.matrix(rep(1, dim(data)[1])))
}


createAmat<-function(dmat,omegas,alpha){
  # creates correl matrix from design data, scale params omegas,
  # and adds nugget to diag based on ratio alpha
  amat <- get_corr_matrix(dmat, omegas=omegas)
  amat<-amat * alpha
  amat<- amat + diag( 1-alpha,dim(amat)[1])
  return(amat)
}


createDMatrix<-function(data,n,p){
  # takes data, number of obs n, and number of variables p
  # returns n x p scaled matrix within inputs between 0 and 1
  dmat<-matrix(data,n,p)
  for (i in 1:p){
    dmat[,i]<-scale(dmat[,i])
  }
  return  (dmat)
}


get_sq_dist<-function(index,design){
  dsq<- replicate(length(design[,index]), design[,index]**2)
  dcross<-design[,index] %*% t(design[,index])
  return( dsq+t(dsq)-dcross-t(dcross))
}

get_corr_matrix<-function(design, omegas){
  # based on hankin bacco implementation
  # reproduced here to avoid import
  if (!is.matrix(omegas)){
    omegas <- diag(omegas)
  }
  quadmat <- design %*% omegas %*% t(design)
  dsq <- kronecker(diag(quadmat), t(rep(1,dim(design)[1])))
  return(exp(quadmat + t(quadmat) -dsq -t(dsq)))
}

get_t_matrix<-function(design, newdata, omegas){
  # based on hankin bacco implementation
  # reproduced here to avoid import
  if (!is.matrix(omegas)){
    omegas <- diag(omegas)
  }

  quad1 <- design %*% omegas %*% t(newdata)
  quad2 <- newdata %*% omegas %*% t(design)

  quadmat <- design %*% omegas %*% t(design)
  dsq <- kronecker(diag(quadmat), t(rep(1,dim(newdata)[1])))

  n_quadmat <- newdata %*% omegas %*% t(newdata)
  n_dsq <- kronecker(diag(n_quadmat), t(rep(1,dim(design)[1])))

  return(exp(t(quad1) + quad2 - n_dsq- t(dsq)))
}


createHInteractionFunction<-function( indexPairs){
  createHInteraction<-function(data){
    interactions<-c()
    for (i in seq_along(indexPairs)){
      interactions<-cbind(interactions, data[,indexPairs[[i]][1]]*
                            data[,indexPairs[[i]][2]])
    }
    hmat<-cbind(rep(1,dim(data)[1]),data, interactions)
    return(hmat)
  }
  return(createHInteraction)
}


create_validation_sample<-function(design,M, hiRange, loRange,flag_scaleUp=T){
  # Create a validation sample.
  # half are pertubations of the existing design
  # the other half are space-filling additions.
  # this gives a mix of short and long 'distances'

  valid_d<-design[sample(nrow(design),size=floor(M/2),replace=F),]
  valid_d<-apply(valid_d, 2,jitter,factor=3)
  if (max(design)>1 || min(design)<0){
    Ssample<-scaleDownLHS(design,hiRange, loRange)
  }
  else{
    valid_d[valid_d<0]<-0.0001
    valid_d[valid_d>1]<-0.9999
    Ssample<-design
  }
  Ssample<-augmentLHS(Ssample,m=ceiling(M/2))


  if (flag_scaleUp==T){
    valid_d<-scaleUpLHS(valid_d,hiRange,loRange)
    Ssample<-scaleUpLHS(Ssample,hiRange,loRange)
  }
  valid_d<-rbind(valid_d,Ssample[(dim(design)[1]+1):dim(Ssample)[1],])

  return(valid_d)
}

createNewPred<-function(train, new_data, n,new_n,p){
  # takes original training design, new points to be estimated
  # the number of new obs new_n and the number variables p
  # returns a input matrix of the new data scaled to the same input
  # space as the original data
  new_data<-matrix(new_data,new_n,p)
  train<-matrix(train,n,p)
  for (i in 1:p){
    scale1<-min(train[,i])
    scale2<-max(train[,i]-scale1)
    new_data[,i]<-(new_data[,i]-scale1)/scale2
  }
  return(new_data)
}


