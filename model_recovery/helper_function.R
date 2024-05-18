sim_m1_RL_sliced2_halfStudentT_reparam_confEs_norm <-function(NS,params,group='All',seq_data){
  NT = dim(seq_data$noxInputTrans)[2]
  noxInputTrans = seq_data$noxInputTrans
  TrialTypeAllWMissing = seq_data$TrialTypeAllWMissing
  IndexMissingAll = seq_data$IndexMissingAll
  RatingConf = seq_data$RatingConf
  PredConf = seq_data$PredConf
  confData = matrix(0,NS,NT)
  
  PercIndexArray = seq_data$PercIndexArray
  PredIndexArray= seq_data$PredIndexArray
  
  simulatedData = matrix(NA,NS,NT)
  painRating = matrix(NA,NS,dim(PercIndexArray)[2])
  painPred = matrix(NA,NS,dim(PredIndexArray)[2])
  
  alpha=params$alpha
  gam=params$gam
  eta=params$eta
  E0=params$E0
  cs=params$cs
  
  
  for (s in 1:NS){
    confData[s,PercIndexArray[s,]] = RatingConf[s,]
    confData[s,PredIndexArray[s,]] = PredConf[s,]
    confData[s,IndexMissingAll[s,]] = NA
    
    E = E0[s]
    tmp_store=rep(0,2)
    for (t in 1:NT){
      P = (1-gam[s])*noxInputTrans[s,t] + gam[s]*E;
      PE = P - E;
      E = E + alpha[s]*PE;
      
      tmp_store[1] = P;
      tmp_store[2] = E;
      
      if (!(t %in% IndexMissingAll[s,])){
        simulatedData[s,t]=rnorm(1,tmp_store[TrialTypeAllWMissing[s,t]],eta[s]*exp((1-confData[s,t])/cs[s]))
      }
    }
    
    painRating[s,] = simulatedData[s,PercIndexArray[s,]]
    painPred[s,] = simulatedData[s,PredIndexArray[s,]]
  }
  
  colnames(simulatedData)=NULL
  colnames(simulatedData)=NULL
  colnames(painPred)=NULL
  return(list(simulatedData,painRating,painPred))
}

sim_m2_RL_sliced2_halfStudentT_reparam_confEs_np_norm <-function(NS,params,group='All',seq_data){
  NT = dim(seq_data$noxInputTrans)[2]
  noxInputTrans = seq_data$noxInputTrans
  TrialTypeAllWMissing = seq_data$TrialTypeAllWMissing
  IndexMissingAll = seq_data$IndexMissingAll
  RatingConf = seq_data$RatingConf
  PredConf = seq_data$PredConf
  confData = matrix(0,NS,NT)
  
  PercIndexArray = seq_data$PercIndexArray
  PredIndexArray= seq_data$PredIndexArray
  
  simulatedData = matrix(NA,NS,NT)
  painRating = matrix(NA,NS,dim(PercIndexArray)[2])
  painPred = matrix(NA,NS,dim(PredIndexArray)[2])
  
  alpha=params$alpha
  eta=params$eta
  E0=params$E0
  cs=params$cs
  
  for (s in 1:NS){
    confData[s,PercIndexArray[s,]] = RatingConf[s,]
    confData[s,PredIndexArray[s,]] = PredConf[s,]
    confData[s,IndexMissingAll[s,]] = NA
    
    E = E0[s]
    tmp_store=rep(0,2)
    for (t in 1:NT){
      P = noxInputTrans[s,t]
      PE = P - E;
      E = E + alpha[s]*PE;
      
      tmp_store[1] = P;
      tmp_store[2] = E;
      
      if (!(t %in% IndexMissingAll[s,])){
        simulatedData[s,t]=rnorm(1,tmp_store[TrialTypeAllWMissing[s,t]],eta[s]*exp((1-confData[s,t])/cs[s]))
      }
    }
    
    painRating[s,] = simulatedData[s,PercIndexArray[s,]]
    painPred[s,] = simulatedData[s,PredIndexArray[s,]]
  }
  
  colnames(simulatedData)=NULL
  colnames(simulatedData)=NULL
  colnames(painPred)=NULL
  return(list(simulatedData,painRating,painPred))
}

sim_m3_KF_sliced2_halfStudentT_reparam_confEs_norm <-function(NS,params,group='All',seq_data){
  NT = dim(seq_data$noxInputTrans)[2]
  noxInputTrans = seq_data$noxInputTrans
  TrialTypeAllWMissing = seq_data$TrialTypeAllWMissing
  IndexMissingAll = seq_data$IndexMissingAll
  RatingConf = seq_data$RatingConf
  PredConf = seq_data$PredConf
  confData = matrix(0,NS,NT)
  
  PercIndexArray = seq_data$PercIndexArray
  PredIndexArray= seq_data$PredIndexArray
  
  simulatedData = matrix(NA,NS,NT)
  painRating = matrix(NA,NS,dim(PercIndexArray)[2])
  painPred = matrix(NA,NS,dim(PredIndexArray)[2])
  
  eps = params$eps
  psi = params$psi
  eta = params$eta
  xi = params$xi
  E0 = params$E0
  w0 = params$w0
  cs = params$cs
  
  eps_sq = eps^2
  psi_sq = psi^2
  eta_sq = eta^2
  w0_sq = w0^2
  
  for (s in 1:NS){
    confData[s,PercIndexArray[s,]] = RatingConf[s,]
    confData[s,PredIndexArray[s,]] = PredConf[s,]
    confData[s,IndexMissingAll[s,]] = NA
    
    E = E0[s]
    w = w0_sq[s]
    tmp_store=rep(0,2)
    
    for (t in 1:NT){
      gam = eps_sq[s]/(eps_sq[s]+psi_sq[s]+w);
      P = (1-gam)*noxInputTrans[s,t] + gam*E;
      alpha = w/(psi_sq[s]+w);
      PE = P - E;
      E = E + alpha*PE;
      w = w * (eps_sq[s]+psi_sq[s])/(eps_sq[s]+psi_sq[s]+w) + eta_sq[s];
      
      tmp_store[1] = P;
      tmp_store[2] = E;
      
      if (!(t %in% IndexMissingAll[s,])){
        simulatedData[s,t]=rnorm(1,tmp_store[TrialTypeAllWMissing[s,t]],xi[s]*exp((1-confData[s,t])/cs[s]))
      }
    }
    
    painRating[s,] = simulatedData[s,PercIndexArray[s,]]
    painPred[s,] = simulatedData[s,PredIndexArray[s,]]
  }
  
  colnames(simulatedData)=NULL
  colnames(simulatedData)=NULL
  colnames(painPred)=NULL
  return(list(simulatedData,painRating,painPred))
  
}

sim_m4_KF_sliced2_halfStudentT_reparam_confEs_np_norm <-function(NS,params,group='All',seq_data){
  NT = dim(seq_data$noxInputTrans)[2]
  noxInputTrans = seq_data$noxInputTrans
  TrialTypeAllWMissing = seq_data$TrialTypeAllWMissing
  IndexMissingAll = seq_data$IndexMissingAll
  RatingConf = seq_data$RatingConf
  PredConf = seq_data$PredConf
  confData = matrix(0,NS,NT)
  
  PercIndexArray = seq_data$PercIndexArray
  PredIndexArray= seq_data$PredIndexArray
  
  simulatedData = matrix(NA,NS,NT)
  painRating = matrix(NA,NS,dim(PercIndexArray)[2])
  painPred = matrix(NA,NS,dim(PredIndexArray)[2])
  
  psi = params$psi
  eta = params$eta
  xi = params$xi
  E0 = params$E0
  w0 = params$w0
  cs = params$cs
  
  psi_sq = psi^2
  eta_sq = eta^2
  w0_sq = w0^2
  
  for (s in 1:NS){
    confData[s,PercIndexArray[s,]] = RatingConf[s,]
    confData[s,PredIndexArray[s,]] = PredConf[s,]
    confData[s,IndexMissingAll[s,]] = NA
    
    E = E0[s]
    w = w0_sq[s]
    tmp_store=rep(0,2)
    
    for (t in 1:NT){
      P = noxInputTrans[s,t]
      alpha = w/(psi_sq[s]+w);
      PE = P - E;
      E = E + alpha*PE;
      w = w * (psi_sq[s])/(psi_sq[s]+w) + eta_sq[s];
      
      tmp_store[1] = P;
      tmp_store[2] = E;
      
      if (!(t %in% IndexMissingAll[s,])){
        simulatedData[s,t]=rnorm(1,tmp_store[TrialTypeAllWMissing[s,t]],xi[s]*exp((1-confData[s,t])/cs[s]))
      }
    }
    
    painRating[s,] = simulatedData[s,PercIndexArray[s,]]
    painPred[s,] = simulatedData[s,PredIndexArray[s,]]
  }
  
  colnames(simulatedData)=NULL
  colnames(simulatedData)=NULL
  colnames(painPred)=NULL
  return(list(simulatedData,painRating,painPred))
  
}

sim_m5_C_sliced2_halfStudentT_reparam_confEs_norm<-function(NS,params,group='All',seq_data){
  NT = dim(seq_data$noxInputTrans)[2]
  noxInputTrans = seq_data$noxInputTrans
  TrialTypeAllWMissing = seq_data$TrialTypeAllWMissing
  IndexMissingAll = seq_data$IndexMissingAll
  RatingConf = seq_data$RatingConf
  PredConf = seq_data$PredConf
  confData = matrix(0,NS,NT)
  
  PercIndexArray = seq_data$PercIndexArray
  PredIndexArray= seq_data$PredIndexArray
  
  simulatedData = matrix(NA,NS,NT)
  painRating = matrix(NA,NS,dim(PercIndexArray)[2])
  painPred = matrix(NA,NS,dim(PredIndexArray)[2])
  
  eta = params$eta
  C = params$C
  cs = params$cs
  
  for (s in 1:NS){
    confData[s,PercIndexArray[s,]] = RatingConf[s,]
    confData[s,PredIndexArray[s,]] = PredConf[s,]
    confData[s,IndexMissingAll[s,]] = NA
    
    tmp_store=rep(0,2)
    
    for (t in 1:NT){
      tmp_store[1] = C[s];
      tmp_store[2] = C[s];
      
      if (!(t %in% IndexMissingAll[s,])){
        simulatedData[s,t]=rnorm(1,tmp_store[TrialTypeAllWMissing[s,t]],eta[s]*exp((1-confData[s,t])/cs[s]))
      }
    }
    painRating[s,] = simulatedData[s,PercIndexArray[s,]]
    painPred[s,] = simulatedData[s,PredIndexArray[s,]]
  }
  
  colnames(simulatedData)=NULL
  colnames(simulatedData)=NULL
  colnames(painPred)=NULL
  return(list(simulatedData,painRating,painPred))
  
}
  
  
generate_seqs_and_conf <- function(NS,base_dir,results_dir){
  results_dir = paste(base_dir,results_dir,sep='')
  print(results_dir)
  plotMe = FALSE
  hr = 10/1.25
  wr = 14/2
  window = 10
  lag_spec=5
  lift_perc = 0.15
  lift_pred = -0.2
  sd_perc = 0.05
  sd_pred = 0.075
  
  ###---- Process transformation coefficients ----
  a_coeffs = read.table(paste(results_dir,'/a_coeffs.csv',sep=''),sep=',',header = TRUE)[,1:2]
  b_coeffs = read.table(paste(results_dir,'/b_coeffs.csv',sep=''),sep=',',header = TRUE)[,1:2]
  colnames(a_coeffs)[2]='a_coeff'
  colnames(b_coeffs)[2]='b_coeff'
  coeffs_df = cbind(a_coeffs,b_coeffs['b_coeff'])
  coeffs_df$source = 'data'
  coeffs_df = coeffs_df[,c(1,4,2,3)]
  
  long_coeffs_df = melt(coeffs_df,id.vars = c("PID",'source'))
  
  ##---- Fit distribution to the coefficients ----
  sink('nul')
  fit_ln_a <- fitdist(coeffs_df$a_coeff, "norm")
  fit_ln_b <- fitdist(coeffs_df$b_coeff, "norm")
  sink()
  mu_a=fit_ln_a$estimate["mean"]
  sd_a=fit_ln_a$estimate["sd"]
  mu_b=fit_ln_b$estimate["mean"]
  sd_b=fit_ln_b$estimate["sd"]
  
  Nsamples = 10000
  a_coeff_draws = data.frame(rnorm(Nsamples,mu_a,sd_a))
  b_coeff_draws = data.frame(rnorm(Nsamples,mu_b,sd_b))
  
  colnames(a_coeff_draws) = 'value'
  a_coeff_draws$PID = 901:(900+dim(a_coeff_draws)[1])
  a_coeff_draws$source='draw'
  a_coeff_draws$variable='a_coeff'
  a_coeff_draws = a_coeff_draws[,c(2:4,1)]
  
  colnames(b_coeff_draws) = 'value'
  b_coeff_draws$PID = 901:(900+dim(b_coeff_draws)[1])
  b_coeff_draws$source='draw'
  b_coeff_draws$variable='b_coeff'
  b_coeff_draws = b_coeff_draws[,c(2:4,1)]
  
  oor_val_bool = (a_coeff_draws$value*0+b_coeff_draws$value>=0 & a_coeff_draws$value*0+b_coeff_draws$value<=100) & (a_coeff_draws$value*12+b_coeff_draws$value>=0 & a_coeff_draws$value*12+b_coeff_draws$value<=100)
  
  a_coeff_draws=a_coeff_draws[oor_val_bool,]
  b_coeff_draws=b_coeff_draws[oor_val_bool,]
  
  coeff_draws_df=rbind(a_coeff_draws,b_coeff_draws)
  if (plotMe){
    coeffs_all_df = rbind(long_coeffs_df,coeff_draws_df)  
  }
  
  unique_PID_draw = sort(unique(coeff_draws_df$PID))
  
  model_data_Intens=data.frame(readRDS(paste(results_dir,'/data_for_stan_lin_dropout','.rds',sep=''))$IntesSeq)
  # model_data=readRDS(paste(results_dir,'/data_for_stan_lin_dropout','.rds',sep=''))
  colnames(model_data_Intens) = paste('t',1:dim(model_data_Intens)[2],sep='')
  
  
  pairs_data = cbind(sample(1:27,NS*1.5,replace=TRUE),sample(unique_PID_draw,NS*1.5,replace=TRUE))
  if (anyDuplicated(pairs_data)){
    pairs_data=pairs_data[-c(anyDuplicated(pairs_data)),]
  }
  pairs_data=pairs_data[1:NS,]
  
  sim_data_pretrans = model_data_Intens[pairs_data[,1],]
  rownames(sim_data_pretrans) = paste('s',(1:dim(sim_data_pretrans)[1]),sep='')
  
  sim_coeffs = coeff_draws_df[coeff_draws_df$PID %in% pairs_data[,2],]
  a_sim_coeffs = sim_coeffs[sim_coeffs$variable=='a_coeff',]$value
  b_sim_coeffs =sim_coeffs[sim_coeffs$variable=='b_coeff',]$value
  
  sim_data_trans=a_sim_coeffs*sim_data_pretrans+b_sim_coeffs
  rownames(sim_data_trans) = paste('s',(1:dim(sim_data_trans)[1]),sep='')
  
  colnames(sim_data_trans) = 1:length(colnames(sim_data_trans))
  rownames(sim_data_trans) = 1:length(rownames(sim_data_trans))
  sim_data_trans_array = as.matrix(sim_data_trans)
  sim_data_trans <- reshape2::melt(t(sim_data_trans))
  colnames(sim_data_trans)[1:2]=c('trial','sub')
  
  noxInputTrans = sim_data_trans_array ############## (1) ##############
  
  ##--- Plot sequneces ----
  if (plotMe){
    pl = vector(mode = "list", length = NS)
    for (s in 1:NS){
      data = sim_data_trans[sim_data_trans$sub==s,]
      pl[[s]] = ggplot(data,aes(x=1:320,y=value))+geom_line(size=0.5)+
        ylim(0,100)+
        xlab('Trial')
    }
    eps_ar_plot = grid.arrange(grobs=pl,ncol=4,guide_legend(nrow=2, byrow=TRUE))
    ggplot2::ggsave(paste(results_dir,"/transformed_simulted",'',".png",sep=''), plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
    
  }
  
  ###----- Confidence -----
  confData = data.frame()
  
  PercIndexArrayBig = matrix(0,NS,160)
  PredIndexArrayBig = matrix(0,NS,156)
  RatingConfBig = matrix(0,NS,160)
  PredConfBig = matrix(0,NS,156)
  IndexMissingAllBig = matrix(0,NS,4)
  TrialTypeAllWMissingBig = matrix(0,NS,320)
  
  Tn=320
  for (subno in 1:length(unique(sim_data_trans$sub))){
    type_seq = 1*(runif(Tn,0,1)<0.5)+1
    if (sum(type_seq==1)>Tn/2){
      ind_too_much = which(type_seq==1)
      ind_to_change = ind_too_much[sample.int(n = length(ind_too_much),size = sum(type_seq==1) - Tn/2)]
      type_seq[ind_to_change]=2
    }else if(sum(type_seq==2)>Tn/2){
      ind_too_much = which(type_seq==2)
      ind_to_change = ind_too_much[sample.int(n = length(ind_too_much),size = sum(type_seq==2) - Tn/2)]
      type_seq[ind_to_change]=1
    }
    
    TrialTypeAllWMissingBig[subno,] = type_seq
    tmp_series_perc = sim_data_trans[sim_data_trans$sub==subno,]$value[type_seq==1]
    tmp_series_pred = sim_data_trans[sim_data_trans$sub==subno,]$value[type_seq==2]
    
    conf_store_perc=vector(length = length(tmp_series_perc))
    conf_store_pred=vector(length = length(tmp_series_pred))
    for (t in 1:length(conf_store_perc)){
      if (t<window){
        # conf_store_perc[t] = runif(1,0,1)
        conf_store_perc[t] = as.numeric(acf(tmp_series_perc[t:(window+t-1)],lag=lag_spec,pl=FALSE)[[1]])[2]+rnorm(1,0,sd_perc)+lift_perc
        if (is.na(conf_store_perc[t])  & t>1){
          conf_store_perc[t] = conf_store_perc[t-1]
        }else if(is.na(conf_store_perc[t])  & t==1){
          conf_store_perc[t] = runif(1,0,0.9)
        }
        # conf_store_pred[t] = runif(1,0,0.85)
        conf_store_pred[t] = as.numeric(acf(tmp_series_pred[t:(window+t-1)],lag=lag_spec,pl=FALSE)[[1]])[2]+rnorm(1,0,sd_pred)+lift_pred
        if (is.na(conf_store_pred[t]) & t>1){
          conf_store_pred[t] = conf_store_pred[t-1]
        }else if(is.na(conf_store_pred[t])  & t==1){
          conf_store_pred[t] = runif(1,0,0.7)
        }
      }else{
        conf_store_perc[t] = as.numeric(acf(tmp_series_perc[(t-window):t],lag=lag_spec,pl=FALSE)[[1]])[2]+rnorm(1,0,sd_perc)+lift_perc
        conf_store_pred[t] = as.numeric(acf(tmp_series_pred[(t-window):t],lag=lag_spec,pl=FALSE)[[1]])[2]+rnorm(1,0,sd_pred)+lift_pred
      }
    }
    conf_store_perc = 0.5+conf_store_perc/2
    conf_store_pred = 0.5+conf_store_pred/2
    # plot(1:160,conf_store_perc,type='l',ylim=c(0,1),col='red')
    # lines(1:160,conf_store_pred,type='l',ylim=c(0,1),col='green')
    options(warn=2)
    ma <- function(x, n = sample.int(round(window*0.75),1,replace=TRUE)){stats::filter(x, rep(1 / n, n), sides = 2)}
    
    conf_store_perc_smooth = as.numeric(ma(conf_store_perc))
    conf_store_pred_smooth = as.numeric(ma(conf_store_pred))
    
    na_ind = which(is.na(conf_store_perc_smooth))
    na_sum = sum(is.na(conf_store_perc_smooth))
    
    if((na_sum)){
      if (length(na_ind)>2){
        chg_pt = which(diff(na_ind)>1)+1
        left_ind = na_ind[1:(chg_pt-1)]
        right_ind = na_ind[chg_pt:length(na_ind)]
        conf_store_perc_smooth[left_ind]=conf_store_perc_smooth[left_ind+length(left_ind)]+rnorm(length(left_ind),0,0.05)
        conf_store_perc_smooth[right_ind]=conf_store_perc_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }else if(length(na_ind)==2){
        # chg_pt = which(diff(na_ind)>1)+1
        left_ind = na_ind[1]
        right_ind = na_ind[2]
        conf_store_perc_smooth[left_ind]=conf_store_perc_smooth[left_ind+length(left_ind)]+rnorm(length(left_ind),0,0.05)
        conf_store_perc_smooth[right_ind]=conf_store_perc_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }else{
        right_ind = na_ind
        conf_store_perc_smooth[right_ind]=conf_store_perc_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }
      
    }
    na_ind = which(is.na(conf_store_pred_smooth))
    na_sum = sum(is.na(conf_store_pred_smooth))
    if((na_sum)){
      if (length(na_ind)>2){
        chg_pt = which(diff(na_ind)>1)+1
        left_ind = na_ind[1:(chg_pt-1)]
        right_ind = na_ind[chg_pt:length(na_ind)]
        conf_store_pred_smooth[left_ind]=conf_store_pred_smooth[left_ind+length(left_ind)]+rnorm(length(left_ind),0,0.05)
        conf_store_pred_smooth[right_ind]=conf_store_pred_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }else if(length(na_ind)==2){
        # chg_pt = which(diff(na_ind)>1)+1
        left_ind = na_ind[1]
        right_ind = na_ind[2]
        conf_store_pred_smooth[left_ind]=conf_store_pred_smooth[left_ind+length(left_ind)]+rnorm(length(left_ind),0,0.05)
        conf_store_pred_smooth[right_ind]=conf_store_pred_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }
      else{
        right_ind = na_ind
        conf_store_pred_smooth[right_ind]=conf_store_pred_smooth[right_ind-length(right_ind)]+rnorm(length(right_ind),0,0.05)
      }
      
    }
    
    
    for (c in 1:4){
      ind_preds = which(type_seq[(1+80*(c-1)):(80*c)]==2)
      ind_preds[length(ind_preds)]
      type_seq[(1+80*(c-1)):(80*c)][ind_preds[length(ind_preds)]]=NA
    }
    conf_store_pred_smooth[which(is.na(type_seq[type_seq==2 | is.na(type_seq)]))]=NA
    
    # IndexMissingAll = which(is.na(conf_store_pred_smooth))
    IndexMissingAll = which(is.na(type_seq))
    IndexMissingAllBig[subno,] = IndexMissingAll
    PercIndexArray = which(type_seq==1)
    PredIndexArray = which(type_seq==2)
    
    PercIndexArrayBig[subno,] = PercIndexArray
    PredIndexArrayBig[subno,] = PredIndexArray
    
    conf_store_pred_smooth = conf_store_pred_smooth[!is.na(conf_store_pred_smooth)]
    
    PercConf = conf_store_perc_smooth
    PredConf = conf_store_pred_smooth
    
    # plot(PercIndexArray,PercConf,type='l',ylim=c(0,1),col='red')
    # lines(PredIndexArray,PredConf,type='l',ylim=c(0,1),col='green')
    
    
    data_tmp = cbind(rep(subno,length(PercConf)),PercConf,PercIndexArray,rep('perception',length(PercConf)))
    data_tmp = rbind(data_tmp,cbind(rep(subno,length(PredConf)),PredConf,PredIndexArray,rep('prediction',length(PredConf))))
    data_tmp = data.frame(data_tmp)
    colnames(data_tmp) = c('sub','value','trial','type')
    data_tmp$sub = as.numeric(data_tmp$sub)
    data_tmp$value = as.numeric(data_tmp$value)
    data_tmp$trial = as.numeric(data_tmp$trial)
    
    confData = rbind(confData,data_tmp)
    
    RatingConfBig[subno,] = PercConf
    PredConfBig[subno,] = PredConf
    options(warn=0)
  }
  if (plotMe){
    pl = vector(mode = "list", length = NS)
    for (s in 1:NS){
      data = confData[confData$sub==s,]
      pl[[s]] = ggplot(data,aes(x=trial,y=value,col=type),show.legend = FALSE)+geom_line(size=0.75,show.legend = FALSE)+
        ylim(0,1)+
        xlab('Trial')+
        labs(col='')
      # plot(1:320,sim_data_trans[s,],type='l')
    }
    eps_ar_plot = grid.arrange(grobs=pl,ncol=4,guide_legend(nrow=2, byrow=TRUE))
    ggplot2::ggsave(paste(results_dir,"/conf_simulated",'',".png",sep=''), plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
  }
  
  return (list(noxInputTrans=noxInputTrans,RatingConf=RatingConfBig,PredConf=PredConfBig,TrialTypeAllWMissing=TrialTypeAllWMissingBig,PercIndexArray=PercIndexArrayBig,PredIndexArray=PredIndexArrayBig,IndexMissingAll=IndexMissingAllBig))
}

draw_params <- function(modelName, group, NS=30, 
                        range_mult, range_bd,
                        par_names,base_dir='',extra_dir='/',run_instance=''){
  
  csv_fname = dir(paste(base_dir,run_instance,'/',extra_dir,modelName,'/',sep=''),full.names = T,pattern = paste('group_pars.*',group,'.*.csv',sep=''))
  gr_pars = data.frame(read.csv(csv_fname,sep=',',header = TRUE))
  n_pars = dim(gr_pars)[1]/2
  drawnPars = data.frame()
  for (p in 1:length(range_mult)){
    if (range_mult[p]){
      tmp_draws = pnorm(gr_pars[p,]$value+ gr_pars[p+n_pars,]$value*rnorm(NS))*range_mult[p]
    }else{
      tmp_draws = gr_pars[p,]$value+ gr_pars[p+n_pars,]$value*rnorm(NS*100)
      tmp_draws = tmp_draws[tmp_draws>=range_bd[[p]][1] & tmp_draws<=range_bd[[p]][2]]
      tmp_draws= tmp_draws[sample.int(length(tmp_draws),NS,replace = FALSE)]
    }
    drawnPars = rbind(drawnPars,tmp_draws)
  }
  drawnPars = data.frame(t(drawnPars))
  colnames(drawnPars) = par_names
  rownames(drawnPars) = NULL
  
  return(drawnPars)
  
}