################################################
################ coverage Function   ###########
################################################

in_function = function(x,val,alp=0.1){
  CI = quantile(x,c(alp/2,1 - alp/2))
  1*(val <= CI[2] & val >= CI[1])
}

is.stationary = function(coef,lags){
  n_poly = max(lags) + 1
  coef_vec = numeric(n_poly)
  coef_vec[1] = 1
  coef_vec[lags + 1] = -coef
  all(Mod(polyroot(coef_vec)) > 1)
}

box_tran <- function(dat,lam1,off =0){
  if( lam1 != 0){             ## 1use power transformation in lam != 0
    ((dat + off)^lam1 - 1) / lam1
  }else{                     ## use log(data) if lam =0
    log(dat + off)
  }
}

inv_box_tran <- function(dat,lam1,off = 0){
  if( lam1 != 0){             ## 1use power transformation in lam != 0
    (dat*lam1 + 1 )^(1/lam1) - off
  }else{                     ## use log(data) if lam =0
    exp(dat) - off
  }
}


pred_function = function(k,dat,X_loc,Z1_loc,Z2_loc,L_loc1,L_loc2,bet1,gam1,bet2,gam2,
                         lag_times1,lag_times2,V1,V2,sig21,sig22,a12,lam = 0){
  
  s_idx = dat$loc_ind[k]
  X_temp = X_loc[[s_idx]]
  Z1_temp = Z1_loc[[s_idx]]
  Z2_temp = Z2_loc[[s_idx]]
  L1_temp = L_loc1[[s_idx]]
  L2_temp = L_loc2[[s_idx]]
  bet1_temp = bet1[s_idx,]
  gam1_temp = gam1[s_idx,]
  bet2_temp = bet2[s_idx,]
  gam2_temp = gam2[s_idx,]
  n_lags1j = length(lag_times1[[k]])
  n_lags2j = length(lag_times2[[k]])
  t_idx = dat$time_ind[k]
  pred1 = max(0,Z1_update_alt(t_idx,Z1_temp,X_temp, L1_temp,bet1_temp,gam1_temp,1,
                              V1[s_idx],sig21,n_lags1j, lags1))
  
  pred2 = Z2_update_alt(t_idx,Z2_temp,X_temp, L2_temp,bet2_temp,gam2_temp,a12,
                        V1[s_idx],V2[s_idx],sig22,n_lags2j, lags2)
  
  return( c( pred1^2 , inv_box_tran(pred2,lam) ) )
}

pred_function_fut = function(k,dat,X_loc,Z1_loc,Z2_loc,L_loc1,L_loc2,bet1,gam1,bet2,gam2,
                             lag_times1,lag_times2,V1,V2,sig21,sig22,a12,lam = 0){
  
  s_idx = dat$loc_ind[k]
  X_temp = X_loc[[s_idx]]
  Z1_temp = Z1_loc[[s_idx]]
  Z2_temp = Z2_loc[[s_idx]]
  L1_temp = L_loc1[[s_idx]]
  L2_temp = L_loc2[[s_idx]]
  bet1_temp = bet1[s_idx,]
  gam1_temp = gam1[s_idx,]
  bet2_temp = bet2[s_idx,]
  gam2_temp = gam2[s_idx,]
  n_lags1j = length(lag_times1[[k]])
  n_lags2j = length(lag_times2[[k]])
  t_idx = dat$time_ind[k]
  
  pred1 = max(0,Z1_update_alt_fut(t_idx,Z1_temp,X_temp, L1_temp,bet1_temp,gam1_temp,1,
                                  V1[s_idx],sig21,n_lags1j, lags1))
  
  pred2 = Z2_update_alt_fut(t_idx,Z2_temp,X_temp, L2_temp,bet2_temp,gam2_temp,a12,
                            V1[s_idx],V2[s_idx],sig22,n_lags2j, lags2)
  
  return( c( pred1^2 , inv_box_tran(pred2,lam) ) )
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

MCMC_mexico = function(dat,dat_dec ,lags1,lags2 ,reps ,burn ,ns , nt,core_use = 1,thin = 1){
  
  lam = 0
  lag_times1 <- lag_times2 <- vector(mode = "list",length = (ns*nt))
  
  for(i in 1:(ns*nt)){
    lag_poss1 = dat$time_ind[i] + lags1
    idx = which(lag_poss1 <= nt)
    lag_times1[[i]] = lag_poss1[idx]
  }
  
  for(i in 1:(ns*nt)){
    lag_poss2 = dat$time_ind[i] + lags2
    idx = which(lag_poss2 <= nt)
    lag_times2[[i]] = lag_poss2[idx]
  }
  
  
  Z1_loc = lapply(1:ns,function(x){ sqrt(dat$O3[dat$loc_ind == x])  } )
  Z1_time = lapply(1:nt,function(x){ sqrt(dat$O3[dat$time_ind == x]) } )
  Z1_all =  sqrt(dat$O3)
  Z1_loc_old = lapply(1:ns,function(x){ sqrt(dat_dec$O3[dat_dec$loc_ind == x])  } )
  Z1_old =  sqrt(dat_dec$O3)
  
  Z2_loc = lapply(1:ns,function(x){ box_tran(dat$PM10[dat$loc_ind == x],lam) } )
  Z2_time = lapply(1:nt,function(x){ box_tran(dat$PM10[dat$time_ind == x],lam) } )
  Z2_all =  box_tran(dat$PM10,lam)
  Z2_loc_old = lapply(1:ns,function(x){ box_tran(dat_dec$PM10[dat_dec$loc_ind == x],lam) } )
  Z2_old =  box_tran(dat_dec$PM10,lam)
  
  
  X_all = as.matrix(cbind(1,scale(dat[,c("RH","TMP")],scale=FALSE)))
  X_all = X_all[c(1:(ns*24),1:(nt*ns - ns*24)),]
  X_loc = lapply(1:ns,function(x){ X_all[dat$loc_ind == x,] } )
  XtX_loc = lapply(X_loc,function(X) t(X) %*% X )
  X_time = lapply(1:nt,function(x){ X_all[dat$time_ind == x,] } )
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
  LtL_loc1 = lapply(L_loc1,function(X) t(X) %*% X )
  L_time1 = lapply(1:nt,function(x){ L_all1[dat$time_ind == x,] } )
  
  L_loc2 = lapply(1:ns,function(x){ L_all2[dat$loc_ind == x,] } )
  LtL_loc2 = lapply(L_loc2,function(X) t(X) %*% X )
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
  sig21 = numeric(reps + burn) ; sig21[1] = 1
  sig22 = numeric(reps + burn) ; sig22[1] = 1
  tau21 = numeric(reps + burn) ; tau21[1] = 1
  tau22 = numeric(reps + burn) ; tau22[1] = 1
  V1 = matrix(0,ncol = ns,nrow= (reps +burn))
  V2 = matrix(0,ncol = ns,nrow= (reps +burn))
  a12 = numeric(reps + burn) 
  a11 = rep(1,reps+burn)
  
  for(i in 2:(reps + burn)){
    
    ############## Likelihood Variance
    
    sig21[i] = sig21_update(Z1_loc, X_loc,bet1[[i-1]], L_loc1, gam1[[i-1]], a11[i-1],V1[i-1,],as1,bs1,ns,nt)
    
    sig22[i] = sig22_update(Z2_loc, X_loc,bet2[[i-1]], L_loc2,gam2[[i-1]], a12[i-1],V1[i-1,],V2[i-1,],as2,bs2,ns,nt)
    
    ############## Update Beta
    
    bet1[[i]] = bet1_update(Z1_loc,X_loc,L_loc1,Sig_inv_b1[[i-1]],bet10[i-1,], gam1[[i-1]],
                            a11[i-1],V1[i-1,], sig21[i], XtX_loc , p,ns,nt)
    bet2[[i]] = bet2_update(Z2_loc,X_loc,L_loc2,Sig_inv_b2[[i-1]] ,bet20[i-1,], gam2[[i-1]], 
                            a12[i-1], V1[i-1,],V2[i-1,],sig22[i],XtX_loc , p,ns,nt)
    bet10[i,] = bet10_update(Sig_inv_b1[[i-1]],bet1[[i]],S_b1_inv, Smb1, ns)
    bet20[i,] = bet20_update(Sig_inv_b2[[i-1]],bet2[[i]] ,S_b2_inv, Smb2, ns)
    Sig_inv_b1[[i]] = Sig10b_update( bet1[[i]], bet10[i,], 1e-3 * diag(p) , p+1 , ns)
    Sig_inv_b2[[i]] = Sig20b_update( bet2[[i]], bet20[i,], 1e-3 * diag(p) , p+1 , ns)
    
    ############## Update Gamma
    
    gam1[[i]] = gam1_update(Z1_loc,X_loc,L_loc1,Sig_inv_g1[[i-1]] ,gam10[i-1,],
                            bet1[[i]], a11[i-1], V1[i-1,], sig21[i],LtL_loc1 , n_lags1,ns,nt)
    gam2[[i]] = gam2_update(Z2_loc,X_loc,L_loc2,Sig_inv_g2[[i-1]] ,gam20[i-1,],
                            bet2[[i]], a12[i-1],  V1[i-1,], V2[i-1,],sig22[i],LtL_loc2 , n_lags2,ns,nt)
    gam10[i,] = gam10_update(Sig_inv_g1[[i-1]],gam1[[i]],S_g1_inv, Smg1, ns)
    gam20[i,] = gam20_update(Sig_inv_g2[[i-1]],gam2[[i]],S_g2_inv, Smg2, ns)
    
    Sig_inv_g1[[i]] = Sig10g_update( gam1[[i]], gam10[i,], 1e-3 * diag(n_lags1) , n_lags1 + 1 , ns)
    Sig_inv_g2[[i]] = Sig20g_update( gam2[[i]], gam20[i,], 1e-3 * diag(n_lags2) ,  n_lags2 + 1 , ns)
    
    ############## Update V
    
    V1[i,] = V1_update(Z1_time,Z2_time, X_time, L_time1, L_time2, a_11 = a11[i-1], a_12 = a12[i-1], 
                       V2[i-1,], sig21 = sig21[i], sig22 = sig22[i], tau21 = tau21[i-1],bet1[[i]],
                       bet2[[i]], gam1[[i]],gam2[[i]],Q, nt, ns)
    V1[i,] = scale(V1[i,],scale=FALSE)
    
    V2[i,] = V2_update(Z2_time, X_time, L_time2,  a_12 = a12[i-1],V1[i,],
                       sig22 = sig22[i], tau22 = tau22[i-1],bet2[[i]],gam2[[i]],Q, nt, ns)
    V2[i,] = scale(V2[i,],scale=FALSE)
    
    tau21[i] = tau21_update(Q, V1[i,], a_t1 = 1, b_t1 = 1, ns)
    
    tau22[i] = tau22_update(Q, V2[i,], a_t2 = 1, b_t2 = 1, ns)
    
    #  a12[i] = a12[i-1]
    a12[i] = a12_update( Z2_loc, X_loc,  L_loc2,bet2[[i]], gam2[[i]], V1[i,], V2[i,],
                         sig22[i], m = 0, s2 = 1,  ns, nt)
    
    cat("\r", i, " of ", reps + burn) 
    flush.console()
    
  }
  
  idx_post = -(1:burn)
  sig21 = sig21[idx_post]; sig22 = sig22[idx_post]; bet1 = bet1[idx_post]; bet2 = bet2[idx_post]
  bet10 = bet10[idx_post,] ; bet20 = bet20[idx_post,];Sig_inv_b1 = Sig_inv_b1[idx_post]
  Sig_inv_b2 = Sig_inv_b2[idx_post] ; gam1 = gam1[idx_post]; gam2 = gam2[idx_post]
  gam10 = gam10[idx_post,]; gam20 = gam20[idx_post,]; Sig_inv_g1 = Sig_inv_g1[idx_post]
  Sig_inv_g2 = Sig_inv_g2[idx_post]; V1 = V1[idx_post,]; V2 = V2[idx_post,] 
  tau21 = tau21[idx_post]; tau22 = tau22[idx_post]; a12 = a12[idx_post]
  
  if(thin > 1){
    idx_keep = which(1:reps %% thin == 1)
    sig21 = sig21[idx_keep]; sig22 = sig22[idx_keep]; bet1 = bet1[idx_keep]; bet2 = bet2[idx_keep]
    bet10 = bet10[idx_keep,] ; bet20 = bet20[idx_keep,];Sig_inv_b1 = Sig_inv_b1[idx_keep]
    Sig_inv_b2 = Sig_inv_b2[idx_keep] ; gam1 = gam1[idx_keep]; gam2 = gam2[idx_keep]
    gam10 = gam10[idx_keep,]; gam20 = gam20[idx_keep,]; Sig_inv_g1 = Sig_inv_g1[idx_keep]
    Sig_inv_g2 = Sig_inv_g2[idx_keep]; V1 = V1[idx_keep,]; V2 = V2[idx_keep,] 
    tau21 = tau21[idx_keep]; tau22 = tau22[idx_keep]; a12 = a12[idx_keep]
  }
  
  idx_pred = which(dat$Hour %in% c(10,15,20))
  
  preds = mclapply( 1:(reps/thin) , function(i){
    out = sapply(idx_pred,function(xx){
      pred_function_fut(k = xx,dat = dat,X_loc = X_loc,Z1_loc = Z1_loc,Z2_loc = Z2_loc,
                        L_loc1 = L_loc1,L_loc2 = L_loc2,bet1 = bet1[[i]],gam1 = gam1[[i]],
                        bet2 = bet2[[i]],gam2 = gam2[[i]],lag_times1 = lag_times1,
                        lag_times2 = lag_times2,V1 = V1[i,],V2= V2[i,],sig21=sig21[i],
                        sig22=sig22[i],a12 = a12[i]) 
      
    })
    cat("\r", i, " of ", reps/thin)
    flush.console()
    
    return(out)
  },mc.cores = core_use)
  
  return(list(sig21 = sig21, sig22 = sig22,bet1 = bet1,bet2 = bet2,bet10 = bet10,bet20 = bet20,
              Sig_inv_b1 = Sig_inv_b1,Sig_inv_b2 = Sig_inv_b2,gam1 = gam1,gam2 = gam2,gam10 = gam10,
              gam20 = gam20, Sig_inv_g1 = Sig_inv_g1,Sig_inv_g2 = Sig_inv_g2,V1 = V1,V2 = V2,
              tau21 = tau21,tau22 = tau22,a12 = a12 ,preds = preds,pred_times = dat$time_ind[idx_pred],
              pred_locs = dat$loc_ind[idx_pred]))
  
}

################################################################################################
#################################### Phase stuff ###############################################
################################################################################################

phase_status = function(ex_O3,ex_PM10){
  
  O3_status = max(ex_O3)
  
  nr = length(ex_O3)
  out = numeric(nr)
  
  if(O3_status == 0){
    
    if( sum(ex_PM10 == 2) >= 2 ){
      out[1:nr] = 2
    } else if( sum(ex_PM10 == 1) >= 2 ){
      out = ifelse(ex_PM10 > 1,2,1)
    } else{
      out = ex_PM10
    }
    
  } else if(O3_status == 1){
    
    if( sum(ex_PM10 == 2) >= 2 ){
      out[1:nr] = 2
    } else{
      out = ifelse(ex_PM10 > 1,2,1)
    }
    
  } else{
    out[1:nr] = 2
  }
  return(out)
}

Phase_MC = function(dat,preds,pred_time,pred_loc ,ns,L1_o= 154,L2_o = 204,L1_p = 214,L2_p = 354){
  
  t_un = unique(pred_time)
  nt = length(t_un)
  stat_regs = dat$Region[1:ns]
  reg = sort(unique(stat_regs))
  nr = length(reg)
  phase = matrix(0,nrow = nt,ncol = nr)
  exceed_O3 = matrix(0,nrow = nt,ncol = nr)
  exceed_PM10 = matrix(0,nrow = nt,ncol = nr)
  
  Lo_law1 = 70 ; Lo_law2 = 95 ; Lp_law = 75
  
  O3_now = sapply(t_un,function(x){preds[1,pred_time == x]})
  
  for(i in 4:nt){
    t_i = t_un[i]
    O3_8 =  apply(cbind(sapply((t_i-7):(t_i-1),function(x){dat$O3[dat$time_ind == x]}),
                        preds[1,pred_time == t_i]),1,mean)
    PM10_24 = apply(cbind(sapply((t_i-23):(t_i-1),function(x){dat$PM10[dat$time_ind == x]}),
                          preds[2,pred_time == t_i]),1,mean)
    
    ex_O3 = tapply(O3_now[,i],stat_regs, function(x){ifelse(any(x > L2_o),2,ifelse(any(x > L1_o),1,0 ))})
    ex_PM10 = tapply(PM10_24,stat_regs, function(x){ifelse(any(x > L2_p),2,ifelse(any(x > L1_p),1,0 ))})
    
    phase[i,] = phase_status(ex_O3,ex_PM10)
    
    ex_O3a = tapply(O3_now[,i],stat_regs, function(x){ifelse(any(x > Lo_law2),1,0)})
    ex_O3b = tapply(O3_8,stat_regs, function(x){ifelse(any(x > Lo_law1),1,0)})
    exceed_PM10[i,] = tapply(PM10_24,stat_regs, function(x){ifelse(any(x > Lp_law),1,0 )})
    exceed_O3[i,] = apply(rbind(ex_O3a,ex_O3b),2,max)
  }
  
  return(list(phase = phase, exceed_O3 = exceed_O3 , exceed_PM10 = exceed_PM10))
  
}

# phase_suspend = function(phase,reps, nt){
#    
#   out = phase
#   for(i in 1:reps){
#     phase_temp = phase[,,i]
#     for(j in 2:nt){
#       phase_temp[,j] = ifelse(all(phase_temp[,j] == 0) )
#     }
#     out[,,i] = phase_temp
#   }
#   
#   return(out)
# }

prob_phase = function(phase,nr = nr,reps = reps){
  
  nt = dim(phase)[2]
  reg_day = array(0,c(nr,nt/3,3))
  reg_hour = array(0,c(nr,nt,3))
  
  tot_day = matrix(0,nrow = nt/3,3)
  tot_hour =matrix(0,nrow = nt,3)
  
  for(i in 1:nt){
    reg_hour[,i,] = t(apply(phase[,i,],1,function(x){table(factor(x, levels = 0:2))})) / reps
    temp = apply(phase[,i,],2,max)
    tot_hour[i,] = table(factor(temp, levels = 0:2)) / reps 
  }
  
  for(i in 1:(nt/3)){
    temp = apply(phase[,((i-1)*3 + 1):(i*3),],c(1,3),max)
    reg_day[,i,] = t(apply(temp,1,function(x){table(factor(x, levels = 0:2))}) / reps)
    tot_day[i,] = table(factor(apply(temp,2,max), levels = 0:2)) / reps
  }
  return(list(reg_day = reg_day,reg_hour = reg_hour,tot_day = tot_day,tot_hour = tot_hour))
}


table_funcs = function(phase,nr = nr,reps = reps){
  
  nt = dim(phase)[2]
  reg_hours = array(0,c(reps,nr,3))
  tot_hours = matrix(0,nrow = reps,ncol = 3)
  
  reg_days = array(0,c(reps,nr,3))
  tot_days = matrix(0,nrow = reps,ncol = 3)
  
  for(i in 1:reps){
    reg_hours[i,,] = t(apply(phase[,,i],1,function(x){table(factor(x, levels = 0:2))})) 
    temp = apply(phase[,,i],2,max)
    tot_hours[i,] = table(factor(temp, levels = 0:2))
  }
  
  for(i in 1:reps){
    this = phase[,,i]
    phase_day = matrix(0,ncol = nr,nrow = nt/3)
    for(j in 1:(nt/3)){
      phase_day[j,] = apply(this[,((j-1)*3 + 1):(j*3)],1,max)
    }
    #t(apply(phase_day,2,function(x){table(factor(x, levels = 0:2))}))
    reg_days[i,,] = t(apply(phase_day,2,function(x){table(factor(x, levels = 0:2))}))
    tot_days[i,] = table(factor(apply(phase_day,1,max), levels = 0:2))
  }
  return(list(reg_days = reg_days,reg_hours = reg_hours,tot_days = tot_days,tot_hours = tot_hours))
}


prob_exceed = function(phase_O3,phase_PM10,nr = nr,reps = reps){
  
  nt = dim(phase_O3)[2]
  reg_O3_day = matrix(0,ncol = nr, nrow = nt/3 )
  reg_PM10_day = matrix(0,ncol = nr, nrow = nt/3 )
  reg_O3_hour = matrix(0,ncol = nr, nrow = nt )
  reg_PM10_hour = matrix(0,ncol = nr, nrow = nt )
  
  tot_O3_day = numeric(nt/3)
  tot_PM10_day = numeric(nt/3)
  tot_O3_hour = numeric(nt)
  tot_PM10_hour = numeric(nt)
  
  for(i in 1:nt){
    reg_O3_hour[i,] = (t(apply(phase_O3[,i,],1,function(x){table(factor(x, levels = 0:1))})) / reps)[,2]
    reg_PM10_hour[i,] = (t(apply(phase_PM10[,i,],1,function(x){table(factor(x, levels = 0:1))})) / reps)[,2]
    temp = apply(phase_O3[,i,],2,max)
    tot_O3_hour[i] = mean(temp)
    temp = apply(phase_PM10[,i,],2,max)
    tot_PM10_hour[i] = mean(temp)
  }
  
  for(i in 1:(nt/3)){
    temp = apply(phase_O3[,((i-1)*3 + 1):(i*3),],c(1,3),max)
    reg_O3_day[i,] = (apply(temp,1,function(x){table(factor(x, levels = 0:1))}) / reps)[2,]
    tot_O3_day[i] = mean(apply(temp,2,max))
    
    temp = apply(phase_PM10[,((i-1)*3 + 1):(i*3),],c(1,3),max)
    reg_PM10_day[i,] = (apply(temp,1,function(x){table(factor(x, levels = 0:1))}) / reps)[2,]
    tot_PM10_day[i] = mean(apply(temp,2,max))
    
  }
  return(list(reg_O3_day = reg_O3_day,reg_O3_hour = reg_O3_hour,
              tot_O3_day = tot_O3_day,tot_O3_hour = tot_O3_hour,
              reg_PM10_day = reg_PM10_day,reg_PM10_hour = reg_PM10_hour,
              tot_PM10_day = tot_PM10_day,tot_PM10_hour = tot_PM10_hour))
  
}
