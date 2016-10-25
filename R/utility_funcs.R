scale<-function(xs){
  return((xs- min(xs))/ max(xs- min(xs)))
}

inv_scale<-function(scl_xs, xs){
  return(scl_xs * max(xs-min(xs)) + min(xs))
}
correl<-function(x, xdash, delta){
  return(exp(-sum(((x-xdash)/delta)**2)))
}

correl_tau<-function(x, xdash, tau){
  return(exp(-sum(((x-xdash)**2/exp(tau)))))
}


giveDeltas<-function(taus){
  # convert taus to deltas
  return(exp(taus/2))
}

giveTaus<-function(deltas){
  #convert deltas to taus
  return ( 2* log (deltas))
}

giveOmegas<-function(deltas){
  #convert deltas to omegas
  return(1/(deltas**2))
}


giveDeltas_oms<-function(omegas){
  return( 1/sqrt(omegas))
}

scaleUpLHS<-function(sample_D, hiRange, loRange){
  hiRange<-as.numeric(hiRange)
  loRange<-as.numeric(loRange)
  for (i in 1:dim(sample_D)[2]){
    sample_D[,i]<-(sample_D[,i]*(hiRange[i]-loRange[i]))+loRange[i]
  }
  return(sample_D)
}

scaleDownLHS<-function(sample, hiRange, loRange){
  for (i in 1:dim(sample)[2]){
    sample[,i]<-(sample[,i]-loRange[i] )/(hiRange[i]-loRange[i])
  }
  return(sample)
}
