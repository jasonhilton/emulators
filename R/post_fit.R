cross_validate_em<-function(data, vdesign,vout, estims){
  if  (data$det==F){
    print("stochastic code normality")
    d1ests<-estimateMean(data,newPoints=data$d1,estimates=estims)
    d1covar<-estimCovarMatrix(data,data$d1,estimates=estims)
    d1resids<-(d1ests-data$fofD)/sqrt(estims$nugget_hat[1])
    qqnorm(d1resids)
    qqline(d1resids)
    mtext("Residuals at original design points")
    plot(d1resids)
    mtext("Residuals at original design points")
    print(shapiro.test(d1resids))
  }
  print("cross-validation")
  dvests<-estimateMean(data,newPoints=vdesign,estimates=estims)
  dvcovar<-estimCovarMatrix(data,newdata=vdesign,estimates=estims)
  #dvresids<-(dvests-log(poptotalV))/sqrt(diag(dvcovar))
  dvresids<-(dvests-vout)/sqrt(diag(dvcovar))
  qqnorm(dvresids)
  qqline(dvresids)
  mtext("Residuals at cross-validation points")
  plot(dvresids)
  mtext("Residuals at cross-validation points")
  # Malonhobis distance
  ER<- chol(dvcovar)
  Q<-backsolve(ER,dvests-vout, transpose=T)
  MD<-t(Q) %*% Q
  multiplier<-(data$N-data$r)/(dim(vdesign)[1] * ( data$N - data$r-2))
  print ("Expectation of MD:")
  print(dim(vdesign)[1])
  print ("Expected 90% intervals:")
  print(qf(c(0.05,0.95), dim(vdesign)[1],data$N-data$r)/multiplier)
  print ("Actual value:")
  print(MD)
  return(MD)
}


cross_validate_em_ret<-function(data, vdesign,vout, estims){
  if  (data$det==F){
    #print("stochastic code normality")
    d1ests<-estimateMean(data,newPoints=data$d1,estimates=estims)
    d1covar<-estimCovarMatrix(data,data$d1,estimates=estims)
    d1resids<-(d1ests-data$fofD)/sqrt(estims$nugget_hat[1])
    res_original<-qplot(1:length(d1resids),d1resids) + geom_hline(yintercept=0) +
      labs(title = "Residuals at original design points") +
      theme(title=element_text()) + xlab("Training Run Number") +
      ylab("Standardised Residual")
    #qqnorm(d1resids)
    #qqline(d1resids)
    #mtext("Residuals at original design points")
    #plot(d1resids)
    #mtext("Residuals at original design points")
    #print(shapiro.test(d1resids))
  }
  #print("cross-validation")
  dvests<-estimateMean(data,newPoints=vdesign,estimates=estims)
  dvcovar<-estimCovarMatrix(data,newdata=vdesign,estimates=estims)
  #dvresids<-(dvests-log(poptotalV))/sqrt(diag(dvcovar))
  dvresids<-(dvests-vout)/sqrt(diag(dvcovar))
  res_cv<-qplot(1:length(dvresids),dvresids) + geom_hline(yintercept=0) +
    labs(title = "Residuals at cross validation points") +
    theme(title=element_text()) + xlab("Cross-Validation Run Number") +
    ylab("Standardised Residual")
  #qqnorm(dvresids)
  #qqline(dvresids)
  #mtext("Residuals at cross-validation points")
  #plot(dvresids)
  #mtext("Residuals at cross-validation points")
  # Malonhobis distance
  ER<- chol(dvcovar)
  Q<-backsolve(ER,dvests-vout, transpose=T)
  MD<-t(Q) %*% Q
  multiplier<-(data$N-data$r)/(dim(vdesign)[1] * ( data$N - data$r-2))
  #print ("Expectation of MD:")
  #print(dim(vdesign)[1])
  E_MD<- dim(vdesign)[1]
  E_95<-qf(c(0.025,0.975), dim(vdesign)[1],data$N-data$r)/multiplier
  #print ("Expected 90% intervals:")
  #print(qf(c(0.025,0.975), dim(vdesign)[1],data$N-data$r)/multiplier)
  #print ("Actual value:")
  #print(MD)
  return(list(E_MD=E_MD,E_95=E_95,MD=MD,res_original=res_original,
              res_cv=res_cv,dvresids = dvresids,
              d1resids = d1resids))
}

plot_1d_preds<-function(data, estims, index, input_range = c(0,1),
                        input_name = "Input",

                        alpha=0.1, num_points=100){
  # plot predictions for all values of one input
  # at central points for all other dimensions
  # takes original data object
  # estims object
  # and the column index of the input to be ploted
  newData = matrix(0.5, num_points, data$p+data$q)
  newData[, index] = seq(0, 1, length.out = num_points)
  mean_hat  <- estimateMean(data, newData, estims)
  covar_hat <- estimCovarMatrix(data, newData, estims)
  var_hat   <- diag(covar_hat)

  df <- data$N - data$r
  pi     <- qt(1 - alpha / 2, df) * sqrt(var_hat)
  upper  <- mean_hat + pi
  lower  <- mean_hat - pi
  ylims  <- c(min(lower) * 0.95, max(upper) * 1.05)
  X_vals <- seq(input_range[1], input_range[2], length.out = num_points)
  plot(X_vals, mean_hat, type = "l", lty = 4,
       ylim=ylims, xlab=input_name)
  polygon(c(newData[, index], rev(newData[, index])),
          c(upper, rev(lower)), density=NA, border=NA, col="lightblue")
  points(X_vals, mean_hat, type = "l")
}

plot_2d_preds<-function(data, estims, indices,
                        input_range1 = c(0,1),
                        input_name1 = "Input1",
                        input_range2 = c(0,1),
                        input_name2 = "Input2",
                        transformation = function(x) {x},
                        num_points=100){
  # plot predictions for all values of one input
  # at central points for all other dimensions
  # takes original data object
  # estims object
  # and the column index of the input to be ploted
  newData = matrix(0.5, num_points*(num_points+1), data$p+data$q)
  newVals<-expand.grid(seq(0, 1, length.out = num_points),
                       seq(0, 1, length.out = num_points+1))

  newData[, indices[1]] = newVals[,1]
  newData[, indices[2]] = newVals[,2]
  mean_hat  <- estimateMean(data, newData, estims)
  mean_hat<-matrix(mean_hat,nrow=100,ncol=101)

  X_vals <- seq(input_range1[1], input_range1[2], length.out = num_points)
  Y_vals <- seq(input_range2[1], input_range2[2], length.out = num_points+1)
  filled.contour(X_vals,Y_vals, transformation(mean_hat), xlab=input_name1,
                 ylab=input_name2)
}

