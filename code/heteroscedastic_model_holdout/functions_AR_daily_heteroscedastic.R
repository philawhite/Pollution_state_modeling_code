################################################
################ coverage Function   ###########
################################################

in_function = function(x,val,alp=0.1){
  CI = quantile(x,c(alp/2,1 - alp/2))
  1*(val <= CI[2] & val >= CI[1])
}

hour_function = function(hours){
  out = hours %% 24
  ifelse(out != 0, out, 24)
}
 

is.stationary = function(coef,lags){
  n_poly = max(lags) + 1
  coef_vec = numeric(n_poly)
  coef_vec[1] = 1
  coef_vec[lags + 1] = -coef
  all(Mod(polyroot(coef_vec)) > 1)
}

box_tran <- function(dat,lam1,off = 1/2){
  if( lam1 != 0){             ## 1use power transformation in lam != 0
    ((dat + off)^lam1 - 1) / lam1
  }else{                     ## use log(data) if lam =0
    log(dat + off)
  }
}

inv_box_tran <- function(dat,lam1,off = 1/2){
  if( lam1 != 0){             ## 1use power transformation in lam != 0
    (dat*lam1 + 1 )^(1/lam1) - off
  }else{                     ## use log(data) if lam =0
    exp(dat) - off
  }
}

################################################
################ Gibbs Functions   #############
################################################

Sig10b_update = function( bet1, bet10, M_b1, nu_b1 , ns){
  
  BtB = Reduce('+',lapply(1:ns,function(j,m,v){(m[j,] - v)%*%t(m[j,] - v )},m=bet1,v=bet10) )
  return( rwish(ns + nu_b1, solve( M_b1 + BtB)) )
  
}

Sig20b_update = function( bet2, bet20, M_b2, nu_b2 , ns){
  BtB = Reduce('+',lapply(1:ns,function(j,m,v){(m[j,] - v)%*%t(m[j,] - v )},m=bet2,v=bet20) )
  return( rwish(ns + nu_b2, solve( M_b2 + BtB)) )
}

Sig10g_update = function( gam1, gam10, M_g1, nu_g1 , ns){
  BtB = Reduce('+',lapply(1:ns,function(j,m,v){(m[j,] - v)%*%t(m[j,] - v )},m=gam1,v=gam10) )
  return( rwish(ns + nu_g1, solve( M_g1 + BtB)) )
}

Sig20g_update = function( gam2, gam20, M_g2, nu_g2 , ns){
  BtB = Reduce('+',lapply(1:ns,function(j,m,v){(m[j,] - v)%*%t(m[j,] - v )},m=gam2,v=gam20) )
  return( rwish(ns + nu_g2, solve( M_g2 + BtB)) )
}

################################################
################ Data
################################################

################################################################
################   Mexico all data
################################################################

MCMC_mexico = function(dat,dat_dec ,lags1,lags2 ,reps ,burn ,ns , nt,nh = 24,
                       logY = FALSE,lam = 0,seed = 1,n_hold = n_hold){
  
  set.seed(seed)
  hold_ind = sample( nt * ns , n_hold )
  hold_ind = hold_ind[order(dat$loc_ind[hold_ind],dat$time_ind[hold_ind])]
  hold_dat = dat[hold_ind,c("obs_ind","Hour","time_ind","loc_ind","PM10","O3")]
  lag_times = vector(mode = "list",length = n_hold)
  dat$O3[hold_ind] = NA
  dat$PM10[hold_ind] = NA
  
  
  for(i in 1:n_hold){
    lag_poss = hold_dat$time_ind[i] + lags1
    idx = which(lag_poss <= nt)
    lag_times[[i]] = lag_poss[idx]
  }
  
  if(logY == FALSE){
    Z1_loc = lapply(1:ns,function(x){ dat$O3[dat$loc_ind == x]  } )
    Z1_time = lapply(1:nt,function(x){dat$O3[dat$time_ind == x] } )
    Z1_all =  dat$O3
    Z1_loc_old = lapply(1:ns,function(x){ dat_dec$O3[dat_dec$loc_ind == x]  } )
    Z1_old = dat_dec$O3
    
    Z2_loc = lapply(1:ns,function(x){ dat$PM10[dat$loc_ind == x] } )
    Z2_time = lapply(1:nt,function(x){ dat$PM10[dat$time_ind == x] } )
    Z2_all =  dat$PM10
    Z2_loc_old = lapply(1:ns,function(x){ dat_dec$PM10[dat_dec$loc_ind == x] } )
    Z2_old = dat_dec$PM10
    
  } else{
    Z1_loc = lapply(1:ns,function(x){ box_tran(dat$O3[dat$loc_ind == x],lam, 1/2)  } )
    Z1_time = lapply(1:nt,function(x){box_tran(dat$O3[dat$time_ind == x],lam, 1/2) } )
    Z1_all =  box_tran(dat$O3,lam, 1/2)
    Z1_loc_old = lapply(1:ns,function(x){ box_tran(dat_dec$O3[dat_dec$loc_ind == x],lam, 1/2)  } )
    Z1_old =  box_tran(dat_dec$O3,lam, 1/2)
    
    Z2_loc = lapply(1:ns,function(x){ box_tran(dat$PM10[dat$loc_ind == x],lam, 1/2) } )
    Z2_time = lapply(1:nt,function(x){ box_tran(dat$PM10[dat$time_ind == x],lam, 1/2) } )
    Z2_all =  box_tran(dat$PM10,lam, 1/2)
    Z2_loc_old = lapply(1:ns,function(x){ box_tran(dat_dec$PM10[dat_dec$loc_ind == x],lam, 1/2) } )
    Z2_old =  box_tran(dat_dec$PM10,lam, 1/2)
    
  }
  
  X_all = as.matrix(cbind(1,scale(dat[,c("RH","TMP")],scale=FALSE)))
  X_loc = lapply(1:ns,function(x){ X_all[dat$loc_ind == x,] } )
  X_time = lapply(1:nt,function(x){ X_all[dat$time_ind == x,] } )
  XtX_loc = lapply(X_loc,function(X) t(X) %*% X )

  p = ncol(X_all)
  n_lags1 = length(lags1)
  n_lags2 = length(lags2)
  
  L_all1 = matrix(0,ncol=length(lags1),nrow=(ns*nt))
  L_all2 = matrix(0,ncol=length(lags2),nrow=(ns*nt))
  
  for(i in 1:(nt*ns)){
    t_ind = dat$time_ind[i] 
    s_ind = dat$loc_ind[i]
    lag_ind1 = (t_ind - lags1)
    lag_ind2 = (t_ind - lags2)
    idx1_2017 = which(1:nt %in% lag_ind1)
    idx2_2017 = which(1:nt %in% lag_ind2)
    idx1_2016 = which((-nt_dec + 1):0 %in% lag_ind1)
    idx2_2016 = which((-nt_dec + 1):0 %in% lag_ind2)
    
    if(length(idx1_2016) == 0){
      L_all1[i,] = Z1_loc[[s_ind]][idx1_2017]
    } else if(length(idx1_2017) ==0 ) {
      L_all1[i,] = Z1_loc_old[[s_ind]][idx1_2016]
    } else{
      L_all1[i,] = c(Z1_loc_old[[s_ind]][idx1_2016],Z1_loc[[s_ind]][idx1_2017])
    }
    
    if(length(idx2_2016) == 0){
      L_all2[i,] = Z2_loc[[s_ind]][idx2_2017]
    } else if(length(idx2_2017) == 0 ) {
      L_all2[i,] = Z2_loc_old[[s_ind]][idx2_2016]
    } else{
      L_all2[i,] = c(Z2_loc_old[[s_ind]][idx2_2016],Z2_loc[[s_ind]][idx2_2017])
    }
    
  }
  
  L_all1 = L_all1[,n_lags1:1]
  L_all2 = L_all2[,n_lags2:1]
  
  L_loc1 = lapply(1:ns,function(x){ L_all1[dat$loc_ind == x,] } )
  L_time1 = lapply(1:nt,function(x){ L_all1[dat$time_ind == x,] } )
  
  L_loc2 = lapply(1:ns,function(x){ L_all2[dat$loc_ind == x,] } )
  L_time2 = lapply(1:nt,function(x){ L_all2[dat$time_ind == x,] } )
  
  S_b1_inv = solve(1e3 * diag(p))
  S_b2_inv = solve(1e3 * diag(p))
  
  m_b1 = rep(0,p)
  m_b2 = rep(0,p) 
  
  Smb1 = S_b1_inv %*% m_b1
  Smb2 = S_b2_inv %*% m_b2
  
  S_g1_inv = solve(1e3 * diag(n_lags1))
  S_g2_inv = solve(1e3 * diag(n_lags2))
  
  m_g1 = rep(0,n_lags1)
  m_g2 = rep(0,n_lags2)
  
  Smg1 = S_g1_inv %*% m_g1
  Smg2 = S_g2_inv %*% m_g2
  
  as1 = 1
  as2 = 1
  bs1 = 1
  bs2 = 1
  preds1 = matrix(0,ncol = n_hold,nrow = (reps + burn))
  preds2 = matrix(0,ncol = n_hold,nrow = (reps + burn)) 
  
  bet1 = vector(mode= "list",length=(reps+burn)) ; bet1[[1]] = matrix(0,ncol=p,nrow=ns)
  bet2 = vector(mode= "list",length=(reps+burn)) ; bet2[[1]] = matrix(0,ncol=p,nrow=ns)
  bet10 = matrix(0,ncol = p,nrow= (reps +burn))
  bet20 = matrix(0,ncol = p,nrow= (reps +burn))
  gam1 = vector(mode= "list",length=(reps+burn)) ; gam1[[1]] = matrix(0,ncol=n_lags1,nrow=ns)
  gam2 = vector(mode= "list",length=(reps+burn)) ; gam2[[1]] = matrix(0,ncol=n_lags2,nrow=ns)
  
  Sig_inv_b1 = vector(mode= "list",length=(reps+burn)) ; Sig_inv_b1[[1]] = 1e-3*diag(p)
  Sig_inv_b2 = vector(mode= "list",length=(reps+burn)) ; Sig_inv_b2[[1]] = 1e-3*diag(p)
  Sig_inv_g1 = vector(mode= "list",length=(reps+burn)) ; Sig_inv_g1[[1]] = 1e-3*diag(n_lags1)
  Sig_inv_g2 = vector(mode= "list",length=(reps+burn)) ; Sig_inv_g2[[1]] = 1e-3*diag(n_lags2)
  
  gam10 = matrix(0,ncol = n_lags1,nrow= (reps +burn))
  gam20 = matrix(0,ncol = n_lags2,nrow= (reps +burn))
  sig21 = matrix(0,ncol = nh,nrow= (reps +burn)) ; sig21[1,] = 1
  sig22 = matrix(0,ncol = nh,nrow= (reps +burn)) ; sig22[1,] = 1
  tau21 = numeric(reps + burn) ; tau21[1] = 1
  tau22 = numeric(reps + burn) ; tau22[1] = 1
  V1 = matrix(0,ncol = ns,nrow= (reps +burn))
  V2 = matrix(0,ncol = ns,nrow= (reps +burn))
  a12 = numeric(reps + burn) 
  a11 = rep(1,reps+burn)
  
  m1 = mean(Z1_all,na.rm = TRUE)
  m2 = mean(Z2_all,na.rm = TRUE)
  
  for(i in 1:n_hold){
    
    s_idx = hold_dat$loc_ind[i]
    t_idx = hold_dat$time_ind[i]
    o_idx = hold_dat$obs_ind[i]
      
    imp1 = m1
    Z1_loc[[s_idx]][t_idx] = imp1
    Z1_time[[t_idx]][s_idx] = imp1
    Z1_all[o_idx] = imp1
    
    imp2 = m2
    Z2_loc[[s_idx]][t_idx] = imp2
    Z2_time[[t_idx]][s_idx] = imp2
    Z2_all[o_idx] = imp2
    
    n_l = length(lag_times[[i]])
    if(n_l > 0 ){
      for(j in 1:n_l){
        t_ind = lag_times[[i]][j]
        o_ind = o_idx + ns * lags1[j] 
        
        L_loc1[[s_idx]][t_ind,j] = imp1
        L_time1[[t_ind]][s_idx,j] = imp1
        L_all1[o_ind,j] = imp1
        
        L_loc2[[s_idx]][t_ind,j]  = imp2
        L_time2[[t_ind]][s_idx,j] = imp2
        L_all2[o_ind,j] = imp2
      }
    }
    
  }
  
  LtL_loc1 = lapply(L_loc1,function(X) t(X) %*% X )
  LtL_loc2 = lapply(L_loc2,function(X) t(X) %*% X )
  
  st = proc.time()
  
  for(i in 2:(reps + burn)){
    
    ############## Likelihood Variance
    
    sig21[i,] = sig21_update(Z1_all, X_all,bet1[[i-1]], L_all1, gam1[[i-1]], a11[i-1],
                             V1[i-1,],dat$Hour,dat$loc_ind,as1,bs1,ns*nt,nh)
    
    sig22[i,] = sig22_update(Z2_all, X_all,bet2[[i-1]], L_all2,gam2[[i-1]], a12[i-1],V1[i-1,],
                            V2[i-1,],dat$Hour,dat$loc_ind,as2,bs2,ns*nt,nh)
    
    ############## Update Beta
    
    bet1[[i]] = bet1_update(Z1_loc,X_loc,L_loc1,Sig_inv_b1[[i-1]],bet10[i-1,], gam1[[i-1]],
                            a11[i-1],V1[i-1,], sig21[i,], p,ns,nt,nh)
    bet2[[i]] = bet2_update(Z2_loc,X_loc,L_loc2,Sig_inv_b2[[i-1]] ,bet20[i-1,], gam2[[i-1]], 
                            a12[i-1], V1[i-1,],V2[i-1,],sig22[i,] , p,ns,nt,nh)
    bet10[i,] = bet10_update(Sig_inv_b1[[i-1]],bet1[[i]],S_b1_inv, Smb1, ns)
    bet20[i,] = bet20_update(Sig_inv_b2[[i-1]],bet2[[i]] ,S_b2_inv, Smb2, ns)
    Sig_inv_b1[[i]] = Sig10b_update( bet1[[i]], bet10[i,], 1e-3 * diag(p) , p+1 , ns)
    Sig_inv_b2[[i]] = Sig20b_update( bet2[[i]], bet20[i,], 1e-3 * diag(p) , p+1 , ns)
  
    ############## Update Gamma
    
    gam1[[i]] = gam1_update(Z1_loc,X_loc,L_loc1,Sig_inv_g1[[i-1]] ,gam10[i-1,],
                            bet1[[i]], a11[i-1], V1[i-1,], sig21[i,], n_lags1,ns,nt,nh)
    gam2[[i]] = gam2_update(Z2_loc,X_loc,L_loc2,Sig_inv_g2[[i-1]] ,gam20[i-1,],
                            bet2[[i]], a12[i-1],  V1[i-1,], V2[i-1,],sig22[i,] , n_lags2,ns,nt,nh)
    gam10[i,] = gam10_update(Sig_inv_g1[[i-1]],gam1[[i]],S_g1_inv, Smg1, ns)
    gam20[i,] = gam20_update(Sig_inv_g2[[i-1]],gam2[[i]],S_g2_inv, Smg2, ns)
    
    Sig_inv_g1[[i]] = Sig10g_update( gam1[[i]], gam10[i,], 1e-3 * diag(n_lags1) , n_lags1 + 1 , ns)
    Sig_inv_g2[[i]] = Sig20g_update( gam2[[i]], gam20[i,], 1e-3 * diag(n_lags2) ,  n_lags2 + 1 , ns)
    
    ############## Update V
    
    V1[i,] = V1_update(Z1_time,Z2_time, X_time, L_time1, L_time2, a_11 = a11[i-1], a_12 = a12[i-1], 
                       V2[i-1,], sig21 = sig21[i,], sig22 = sig22[i,], tau21 = tau21[i-1],bet1[[i]],
                       bet2[[i]], gam1[[i]],gam2[[i]],Q, nt, ns,nh)
    V1[i,] = scale(V1[i,],scale=FALSE)
    
    V2[i,] = V2_update(Z2_time, X_time, L_time2,  a_12 = a12[i-1],V1[i,], sig22 = sig22[i,], 
                       tau22 = tau22[i-1],bet2[[i]],gam2[[i]],Q, nt, ns,nh)
    V2[i,] = scale(V2[i,],scale=FALSE)
    
    tau21[i] = tau21_update(Q, V1[i,], a_t1 = 1, b_t1 = 1, ns)
    
    tau22[i] = tau22_update(Q, V2[i,], a_t2 = 1, b_t2 = 1, ns)
    
    #  a12[i] = a12[i-1]
    a12[i] = a12_update( Z2_all, X_all,  L_all2,bet2[[i]], gam2[[i]], V1[i,], V2[i,],
                         sig22[i,],dat$Hour,dat$loc_ind, m = 0, s2 = 1,  ns, nt,ns*nt,nh)
    
    ############## Impute missing data
    
    for(j in 1:n_hold){
   #   if( (j == 1) | (s_idx != hold_dat$loc_ind[j]) ){
        s_idx = hold_dat$loc_ind[j]
        X_temp = X_loc[[s_idx]]
        Z1_temp = Z1_loc[[s_idx]]
        Z2_temp = Z2_loc[[s_idx]]
        L1_temp = L_loc1[[s_idx]]
        L2_temp = L_loc2[[s_idx]]
        bet1_temp = bet1[[i]][s_idx,]
        gam1_temp = gam1[[i]][s_idx,]
        bet2_temp = bet2[[i]][s_idx,]
        gam2_temp = gam2[[i]][s_idx,]
  #    }
      
      n_lagsj = length(lag_times[[j]])
      t_idx = hold_dat$time_ind[j]
      o_idx = hold_dat$obs_ind[j]
      h_idx = hold_dat$Hour[j]
      lag_hour = hour_function(h_idx + lags1)

      imp1 = max(0,Z1_update_alt(t_idx,Z1_temp,X_temp, L1_temp,bet1_temp,gam1_temp,1,
                           V1[s_idx],sig21[i,],n_lagsj, lags1,h_idx,lag_hour))
      imp2 = Z2_update_alt(t_idx,Z2_temp,X_temp, L2_temp,bet2_temp,gam2_temp,a12[i],
                      V1[s_idx],V2[s_idx],sig22[i,],n_lagsj, lags2,h_idx,lag_hour)
      
      Z1_loc[[s_idx]][t_idx] = imp1
      Z1_time[[t_idx]][s_idx] = imp1
      Z1_all[o_idx] = imp1
    
      Z2_loc[[s_idx]][t_idx] = imp2
      Z2_time[[t_idx]][s_idx] = imp2
      Z2_all[o_idx] = imp2
      
      n_l = length(lag_times[[j]])
      if(n_l > 0 ){
        for(k in 1:n_l){
          t_ind = lag_times[[j]][k]
          o_ind = o_idx + ns * lags1[j] 
          
          L_loc1[[s_idx]][t_ind,k] = imp1
          L_time1[[t_ind]][s_idx,k] = imp1
          L_all1[o_ind,k] = imp1
          
          L_loc2[[s_idx]][t_ind,k]  = imp2
          L_time2[[t_ind]][s_idx,k] = imp2
          L_all1[o_ind,k] = imp2
          
        }
      }
      
      preds1[i,j] = imp1
      preds2[i,j] = imp2
      
    }
    
    LtL_loc1 = lapply(L_loc1,function(X) t(X) %*% X )
    LtL_loc2 = lapply(L_loc2,function(X) t(X) %*% X )
    
    time_its <- (proc.time() - st)[3] / (i)
    time_used <- round((proc.time() - st)[3]/(60),digits=4)
    time_left <- round(time_its * (reps +burn- i )/(60),digits=4)
    
    cat("\r", i, " of ", reps + burn,"||| Time left: ",floor(time_left/60),
        " hours",time_left%%60," minutes") 
    flush.console()
    
  }
  
  return(list(sig21 = sig21[-(1:burn),], sig22 = sig22[-(1:burn),],bet1 = bet1[-(1:burn)],
              bet2 = bet2[-(1:burn)],bet10 = bet10[-(1:burn),],bet20 = bet20[-(1:burn),],
              Sig_inv_b1 = Sig_inv_b1[-(1:burn)],Sig_inv_b2 = Sig_inv_b2[-(1:burn)],
              gam1 = gam1[-(1:burn)],gam2 = gam2[-(1:burn)],gam10 = gam10[-(1:burn),],
              gam20 = gam20[-(1:burn),], Sig_inv_g1 = Sig_inv_g1[-(1:burn)],
              Sig_inv_g2 = Sig_inv_g2[-(1:burn)],V1 = V1[-(1:burn),],V2 = V2[-(1:burn),],
              tau21 = tau21[-(1:burn)],tau22 = tau22[-(1:burn)],a12 = a12[-(1:burn)],
              preds1 = preds1[-(1:burn),] , preds2 = preds2[-(1:burn),],hold_dat = hold_dat))
  
}