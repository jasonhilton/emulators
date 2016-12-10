
get_m_dash_k<-function(C,B,m, xk){
  return( solve(2*C + B) %*% (2*C %*% xk + B %*% m) )
}

get_Q_Xk<-function(C,B, m,x, xk){
  part1<-2 * t(x - xk) %*% C %*% (x - xk)
  part2<-t(x - m) %*% B %*% (x - m)
  return(part1 + part2)
}


get_Rt_k<-function(C,B, xk,m, alpha=1){
  part1<-(alpha) %*% sqrt(det(B))%*%(1/sqrt(det(2*C+B)))
  part2<- exp(-get_Q_Xk(C,B,m,get_m_dash_k(C,B,m,xk),xk)/2)
  return(part1%*%part2)
}

get_Q_x_kl<-function(C,B,m,x,xk,xl){
  part1<-2 * t(x - xk) %*% C %*% (x - xk)
  part2<-2 * t(x - xl) %*% C %*% (x - xl)
  part3<-t(x - m) %*% B %*% (x - m)
  return ( part1+ part2 + part3)
}

get_m_dash_kl<-function(C,B,m,xk,xl){
  return( solve(4*C + B) %*% (2*C %*% xk +2*C %*% xl +   B %*% m) )
}


get_R_t<-function(data, C,B,m,alpha){
  R_t<-rep(0,dim(data$d1)[1])
  for (i in 1:dim(data$d1)[1]){
    R_t[i]<-get_Rt_k(C,B,data$d1[i,],m,alpha)
  }
  return(R_t)
}

get_Q<-function(C, B, m,x,xdash){
  part1<- 2*t(x-xdash) %*% C %*% (x-xdash)
  part2<- t(x-m) %*% B %*% (x-m)
  part3<- t(xdash- m) %*% B %*% (xdash- m)
  return (part1 + part2 +  part3)
}

get_R_ht<-function(data, C, B, m, R_t){
  R_ht<-matrix(0, nrow=dim(data$hmat)[2],ncol=dim(data$hmat)[1])
  for (i in 1:dim(data$d1)[1]){
    R_ht[,i]<- R_t[i] %*% t(c(1,get_m_dash_k(C,B,m,data$d1[i,])))
  }
  return(R_ht)
}

get_R_tt<- function(data, C, B, m,alpha){
  R_tt<- matrix(0, nrow=dim(data$d1)[1], ncol=dim(data$d1)[1])
  for (i in 1:dim(data$d1)[1]){
    for (j in 1:dim(data$d1)[1]){
      mkl<- get_m_dash_kl(C, B, m, data$d1[i,], data$d1[j,])
      R_tt[i, j]<- (((alpha)^2) * sqrt(det(B)) 
                    * (1 / sqrt(det(4 * C + B)))
                    * exp( - get_Q_x_kl(C, B, m, mkl, 
                                        data$d1[i, ], data$d1[j, ]) / 2))
    }
  }
  return(R_tt)
}

get_m_dash_k_U<-function(Bk,B,m,C,xk){
  return (solve(Bk) %*% c( B %*% m , 2 * C %*% xk + B %*% m ))
}


uncertainty<-function(mean_x, covar_x,data,estims,return_intermediates=F){
  m <- as.vector(mean_x)
  B <- matrix(covar_x, nrow = length(mean_x))
  C <- diag(estims$omega_hat)
  R_h <- c(1, m)
  R_hh <- diag(c(0, 1 / diag(B))) + R_h %*% t(R_h)
  R_t <- get_R_t(data,C,B,m,estims$alpha_hat)
  R_ht <- get_R_ht(data, C, B, m, R_t)
  e <- backsolve(estims$R, 
                 backsolve(estims$R,(data$fofD - data$hmat %*% estims$beta_hat), 
                           transpose=T))
  E <- t(R_h) %*% estims$beta_hat + R_t %*% e
  R_tt <- get_R_tt(data, C, B, m, estims$alpha_hat)
  
  G <- backsolve(estims$R, backsolve(estims$R,
                                     data$hmat, transpose=T))
  
  w <- backsolve(estims$R, data$hmat, transpose=T)
  Q <- t(w)%*%w
  K <- chol(Q)
  mm <- rbind(m,m)
  BB <- rbind(cbind( 2*C+B, -2*C),cbind(-2*C, 2*C+B )) 
  
  U <- (estims$alpha_hat) %*% det(B) %*% 1/sqrt(det(BB))
  U_h <- U[1]* R_h
  U_hh <- U[1]*R_hh
  
  
  o <- backsolve(estims$R, R_t, transpose=T)
  RAR <- t(o) %*% o
  Rh_GRt <- (R_h - t(G) %*% R_t)
  temp <- backsolve(K, Rh_GRt, transpose=T)
  part3 <- t(temp) %*% temp
  
  V <- estims$v_hat* ( U - RAR + part3)
  
  ainvR_tt <- backsolve(estims$R,backsolve(estims$R,
                                           R_tt,transpose=T))
  
  part3 <- R_hh - 2*R_ht %*% G + t(G) %*% R_tt %*% G 
  
  I1 <- (estims$v_hat * (1 - sum(diag(ainvR_tt)) + 
                           sum(diag(backsolve(K, backsolve(K, part3, transpose=T))))))                      
  
  
  I2 <- ( t(estims$beta_hat) %*% R_hh %*% estims$beta_hat 
          + 2 * t(estims$beta_hat) %*% R_ht %*% e + t(e) %*% R_tt %*% e)
  
  E_Var <- I1-V + I2 - E**2
  
  if (return_intermediates){
    return(list(E=E, V=V, E_Var=E_Var, I1=I1, I2=I2,e=e,R_t=R_t,R_h=R_h))
  }
  
  return(list(E=E, V=V, E_Var=E_Var, I1=I1, I2=I2))
}



uncertainty_MC<-function(mean,precision,data, estims,reps=100){
  mean<-as.vector(mean)
  MC_means<-rep(0,reps)
  MC_covar<- rep(0,reps)
  I1_MC<-I2_MC <-I3_MC <- I4_MC <- I5_MC <-I6_MC <-rep(0,reps)
  
  x<-rmvnorm(reps,t(mean),solve(precision))
  x_dash<-rmvnorm(reps,t(mean),solve(precision))
  x_dash2<-rmvnorm(reps,t(mean),solve(precision))
  
  for (i in 1:reps){
    M_MC<-estimateMean(data, matrix(x[i,], nrow=1), estims)
    CV<-estimCovarMatrix(data, rbind( x[i,], x_dash[i,]), estims)
    M_dash_MC<-estimateMean(data, matrix(x_dash[i,], nrow=1), estims)
    CV2<-estimCovarMatrix(data, rbind( x[i,], x_dash2[i,]), estims)
    MC_means[i]<-M_MC
    MC_covar[i] <- CV[2]
    I1_MC[i] <- CV[1]
    I2_MC[i] <- M_MC**2
    I3_MC[i] <- CV[2] **2
    I4_MC[i] <- M_MC  * M_dash_MC * CV[2] 
    I5_MC[i] <- CV[2] * CV2[2]
    I6_MC[i] <- M_MC *CV[2]
  }
  
  V_MC <- mean(MC_covar)
  E_MC <- mean(MC_means)
  E_V_MC<-mean(I1_MC)- V_MC + mean(I2_MC)-E_MC **2
  I3_MC<-mean(I3_MC)
  I4_MC<-mean(I4_MC)
  I5_MC<-mean(I5_MC)
  I6_MC<-mean(I6_MC)
  V_V_MC<- 2 * (I3_MC - 2 * I5_MC +V_MC **2) + 4 *  (I4_MC - 2 * I6_MC * E_MC + E_MC**2 * V_MC)
  return(list(E=E_MC, V_E=V_MC,E_V=E_V_MC, V_V_MC=V_V_MC, I1=I1_MC,I2=I2_MC, MC_means=MC_means,MC_covar=MC_covar))
}



compute_U_w<-function(w_i,w_bar,B,C,estims,UA){
  # reduces to U = alpha %*% det(B) %*% 1/sqrt(det(BB)) when w_i=0
  u_vec<-vector("numeric", length=length(w_bar))
  for (ii in 1:length(w_bar)){
    i <- w_bar[ii]
    if (i==0){
      u_vec<-1
      break
    }
    u_vec[ii]<-sqrt(B[i,i]/ (B[i,i] + 4 * C[i,i]))
  }
  U_w <- ( estims$alpha_hat) * prod(u_vec)
  return(U_w)
}


compute_P_w<-function(w_i,w_bar,B, C, m, data, estims){
  P_w <-matrix(0,nrow=data$N,ncol=data$N) 
  for (k in 1:data$N){
    for (l in 1:data$N){
      part1<-vector("numeric", length=length(w_bar))
      for(ii in 1:length(w_bar)){
        i<- w_bar[ii]
        if (i==0){ 
          part1<-1
          break
        }
        part1[ii]<-(B[i,i]/ (B[i,i] + 2 * C[i,i])) * 
          exp( -1/2 * ((2*C[i,i] *B[i,i])/ (B[i,i] + 2 * C[i,i]))* 
                 ((data$d1[k,i]- m[i])**2 + (data$d1[l,i]- m[i])**2))
      }
      part2<-vector("numeric", length=length(w_i))
      for(ii in 1:length(w_i)){
        i<- w_i[ii]
        if (i==0){
          part2<-1
          break
        }
        part2[ii]<-sqrt(B[i,i]/ (B[i,i] + 4 * C[i,i])) * 
          exp( -1/2 * (1/ (B[i,i] + 4 * C[i,i]))* 
                 (4* C[i,i]**2 * (data$d1[k,i]- data$d1[l,i])**2 
                  + 2*C[i,i]*B[i,i]*((data$d1[k,i]- m[i])**2 +(data$d1[l,i]- m[i])**2)) 
          )
      }
      
      P_w[k,l] <- (estims$alpha_hat**2) * prod(part1) * prod(part2)
    }
  }
  return(P_w)
}


get_S_w<-function(w_i,w_bar,B,C,m,data,estims){
  
  S_w<-matrix(0, nrow=dim(data$hmat)[2],ncol=dim(data$hmat)[1])
  for (k in 1:data$r){
    for (l in 1:data$N){
      part2<-vector("numeric", length(data$p+data$q))
      E_h <- get_E_h(w_i,w_bar,B,C,m,data,l)
      for (i in 1:dim(data$d1)[2]){
        part2[i] <- (sqrt(B[i,i]) /sqrt(2*C[i,i] + B[i,i])) * 
          exp(- (1/2) * (2*C[i,i]*B[i,i] / (2* C[i,i] + B[i,i])) * 
                (data$d1[l,i]- m[i])**2) 
      }
      S_w[k,l]<-estims$alpha_hat * E_h[k] *( prod(part2))
    }  
  }
  return(S_w)
}

get_Q_w<- function(data, m, B, w_i){
  Q_w<-matrix(0,nrow=data$r,ncol=data$r)  
  m1<-c(1,m)
  Q_w<-outer(m1,m1)
  if (w_i[1] != 0){
    for ( i in w_i){
      Q_w[i+1,i+1]<-Q_w[i+1,i+1] +1/B[i,i]
    }
  }
  return(Q_w)  
}

do_sensitivity_w<-function(w_i,data, estims, UA, m, B){
  #G<-backsolve(estims$R,backsolve(estims$R,
  #                                data$hmat,transpose=T))
  C<- diag(estims$omega_hat)
  e<-backsolve(estims$R,backsolve(estims$R,
                                  (data$fofD - data$hmat%*%estims$beta_hat),transpose=T))
  W<- solve(t(backsolve(estims$R,data$hmat,transpose=T)) %*% 
              (backsolve(estims$R,data$hmat,transpose=T)))
  
  indices<-seq_len(data$p + data$q) # all possible indexes
  if (test_all_indices(w_i, indices)){
    w_bar <- 0
  } else if (length(w_i) ==1 & w_i[1] == 0) {
    w_bar <- indices
  } else {
    w_bar <- indices[!indices %in% w_i] # the ones not in w_bar
  }
  
  Q_w <- get_Q_w(data,m,B,w_i)
  
  S_w <- get_S_w(w_i,w_bar, B,C, m ,data,estims)   
  P_w <- compute_P_w(w_i, w_bar,B,C,m,data, estims)
  U_w <- compute_U_w(w_i, w_bar, B, C, estims)
  
  E_E_fX_xw2<-get_E_E_fX_xw2(data, estims, U_w,Q_w,P_w, S_w,e,W)  
  E_fX2 <- UA$V + UA$E**2
  return(E_E_fX_xw2 - E_fX2)
  #return(list(E_E_fX_xw2=E_E_fX_xw2,E_fX2=E_fX2,SV=E_E_fX_xw2-E_fX2))
}


test_all_indices<-function(w_i, indices){
  if (! (length( w_i) == length(indices) ) ) {
    return(F)  
  } else { 
    return ( all.equal( sort(w_i), indices))
  }
}

get_E_E_fX_xw2<-function(data,estims, U_w, Q_w,P_w, S_w, e, W){
  part1 <- U_w - sum(diag(backsolve(estims$R, backsolve(estims$R, P_w,transpose=T))))
  temp  <- backsolve(estims$R, backsolve(estims$R, data$hmat, transpose=T))
  part2 <- Q_w - S_w %*%temp - t(S_w %*%temp) + t(temp) %*% P_w %*% temp
  partA <- estims$v_hat * (part1 + sum(diag(W %*% part2)))
  part3 <- t(e) %*% P_w %*% e
  part4 <- 2 * (t(estims$beta_hat) %*% S_w %*% e)
  part5 <- t(estims$beta_hat) %*% Q_w %*% estims$beta_hat
  partB<- part3+part4+part5
  return(partA + partB)
}



get_E_h<-function(w_i,w_bar,B,C,m, data, l){
  E_h<-vector("numeric", length=data$r)
  E_h[1]<-1
  if (w_bar[1]!=0){
    E_h[w_bar+1]<-m[w_bar]
  }
  if(w_i[1]!=0){
    for (ii in 1:length(w_i)){
      i<-w_i[ii]
      E_h[i+1]<- (2*C[i,i]* data$d1[l,i] + B[i,i] *m[i] )/ ( 2* C[i,i] + B[i,i])
    }
  }
  return(E_h)
}

do_sensitivity<-function(data,estims,UA,m,B){
  indices <- seq(data$p + data$q)
  oneways <- combn(indices,1)
  twoways <- combn(indices,2)
  twoways_list <- lapply(1:dim(twoways)[2], function(i) twoways[,i])
  main_effects <- sapply(oneways, do_sensitivity_w, 
                         data=data, estims=estims, UA=UA, m=m,B=B)
  twoway_effects <- lapply(twoways_list,do_sensitivity_w, 
                           data=data, estims=estims, UA=UA, m=m,B=B)
  
  interacts <- mapply(calc_interaction, twoway_effects, twoways_list,
                      MoreArgs= list(main_effects))
  total_var <- do_sensitivity_w(indices, data, estims, UA, m,B)
  #if( (UA$E_Var - total_var - estims$nugget_hat)/ UA$E_Var > 1e-9){
  #  stop("Expected Variance not equal to total sensitivity+nugget")
  #}
  resid_var<- total_var - sum(main_effects, interacts)
  
  Effect_names <- c(oneways, as.character(twoways_list),"resid")
  Sensitivity_variances <- c(main_effects, interacts, resid_var)
  Sensitivity_index <- Sensitivity_variances / total_var
  out_table <- data.frame(Effect_names,Sensitivity_variances, Sensitivity_index) 
  
  return(out_table)
  
}

calc_interaction<-function(two_way,indices, main_effects){
  return(two_way - sum(main_effects[indices]))
}




