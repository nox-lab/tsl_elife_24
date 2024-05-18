###---- Load libraries and set directories
seed = sample(1:99999,size=1)
set.seed(seed)
library(ggplot2)
# library(dplyr)
library(fitdistrplus)
library("reshape2")
library(gridExtra)
# library(tseries)
# library(astsa)
library(rstan)
library(loo)
library(data.table)


if (Sys.info()["sysname"]=="Darwin"){
  library(shinystan)
  # setwd('/Users/onyskj/Library/CloudStorage/OneDrive-UniversityofCambridge/CamNoxLab - OneDrive/tsl_paper/pain_control_git/paper_submission/experiment1')
  setwd('/Users/onyskj/Library/CloudStorage/OneDrive-UniversityofCambridge/CamNoxLab - OneDrive/tsl_paper/elife_revision/')
  source('model_recovery/helper_function.R')
  
  m_from_sh <- commandArgs(trailingOnly = TRUE)[1]
  m_with_sh <- commandArgs(trailingOnly = TRUE)[2]
  list_names <- commandArgs(trailingOnly = TRUE)[3] #group names
  NS <- as.integer(commandArgs(trailingOnly = TRUE)[4])
  specS_fname <- commandArgs(trailingOnly = TRUE)[5]
  s_iter <- as.integer(commandArgs(trailingOnly = TRUE)[6])
  
  if (is.na(m_from_sh)){
    m_from_sh='m1'
    m_with_sh='m1'
    list_names = 'All'
    NS = 30
    specS_fname='mp_rec_spec_test'
    s_iter=1;
  }
  rstan_options(auto_write = TRUE)
}else{
  # setwd('/rds/project/rds-3IOyKgCQu4I/tsl_paper')
  setwd('/rds/project/rds-3IOyKgCQu4I/tsl_elife_rev_jao57')
  source('model_recovery/helper_function.R')
  
  m_from_sh <- commandArgs(trailingOnly = TRUE)[1]
  m_with_sh <- commandArgs(trailingOnly = TRUE)[2]
  list_names <- commandArgs(trailingOnly = TRUE)[3] #group names
  NS <- as.integer(commandArgs(trailingOnly = TRUE)[4])
  specS_fname <- commandArgs(trailingOnly = TRUE)[5]
  s_iter <- as.integer(commandArgs(trailingOnly = TRUE)[6])
  rstan_options(auto_write = FALSE)
}
source('stan_utility.R')
save_outputs=TRUE

rel_dir = 'model_recovery/'
date.time.append <-dget("date_time_append.R")
script_stamp = date.time.append(paste("_",as.character(round(runif(1)*10e2)),sep=''))

specsM_from=transpose(read.table(file=dir(paste(rel_dir,'hpc/specs/models',sep=''),full.names = T,pattern = paste(m_from_sh,'.*',sep='')),header=FALSE,sep=';',quote = ""))
colnames(specsM_from) <- as.character(specsM_from[1, ])
specsM_from <- specsM_from[-1,]
rownames(specsM_from)<-NULL

specsM_with=transpose(read.table(file=dir(paste(rel_dir,'hpc/specs/models',sep=''),full.names = T,pattern = paste(m_with_sh,'.*',sep='')),header=FALSE,sep=';',quote = ""))
colnames(specsM_with) <- as.character(specsM_with[1, ])
specsM_with <- specsM_with[-1,]
rownames(specsM_with)<-NULL

specsS = transpose(read.table(file=dir(paste(rel_dir,'hpc/specs',sep=''),full.names = T,pattern = paste(specS_fname,'.txt',sep='')),header=FALSE,sep=';',quote = ""))
colnames(specsS) <- as.character(specsS[1, ])
specsS <- specsS[-1,]
rownames(specsS)<-NULL

stan_cores = as.integer(specsS$stanCores)
iters = as.integer(specsS$stanIters)
warmups = as.integer(specsS$stanWarmup)
adapt_delta_value = as.double(specsS$stanAdaptDelta)
max_treedepth_value = as.double(specsS$stanMTD)

shortname_from = gsub('.txt','',dir(paste(rel_dir,'hpc/specs/models',sep=''),full.names = F,pattern = paste(m_from_sh,'.*',sep='')))
shortname_with = gsub('.txt','',dir(paste(rel_dir,'hpc/specs/models',sep=''),full.names = F,pattern = paste(m_with_sh,'.*',sep='')))
model_pars_list_from=eval(parse(text=specsM_from$stanParamList))
model_pars_list_with=eval(parse(text=specsM_with$stanParamList))
range_mult_from = eval(parse(text=specsM_from$range_mult))
range_bd_from = eval(parse(text=specsM_from$range_bd))

options(mc.cores = stan_cores)

specfname = gsub('\\.','-',paste('iter',iters,'delta',adapt_delta_value,'mtd',max_treedepth_value,sep='_'));
model_name_with = paste(rel_dir,'models/',shortname_with,'.stan',sep='')

if (m_from_sh==m_with_sh){
  outputdir_prec = paste(rel_dir,'results/param_rec/',m_from_sh,sep='')
  if(!file.exists(outputdir_prec)){dir.create(outputdir_prec,showWarnings = TRUE,recursive = TRUE)}
  outputdir_prec_extra = paste(rel_dir,'results/param_rec_extra/',m_from_sh,'_extra/','s',toString(s_iter),sep='')
  if(!file.exists(outputdir_prec_extra)){dir.create(outputdir_prec_extra,showWarnings = TRUE,recursive = TRUE)}
}
outputdir_mrec = paste(rel_dir,'results/model_rec/',m_from_sh,'/',m_with_sh,sep='')
if(!file.exists(outputdir_mrec)){dir.create(outputdir_mrec,showWarnings = TRUE,recursive = TRUE)}

# Draw parameters for simulation -------------
drawnPars = draw_params(modelName = shortname_from,group='All', NS=NS,range_mult = range_mult_from,
                        range_bd = range_bd_from,par_names = model_pars_list_from,
                        base_dir='model_recovery/prefit_output/')

print(paste('Drawn parameters from ', shortname_from,sep=''))


# Simulate data -------------
seq_data <- NULL
attempt<-1
while( is.null(seq_data) & attempt<10) {
  options(warn=0)
  print(paste('Attempt: ',attempt,sep=''))
  try(
    seq_data <- generate_seqs_and_conf(NS=NS,base_dir = rel_dir,results_dir='data')
  )
  attempt <- attempt + 1
} 
print('Done simulating sequences and confidence ratings')
# seq_data = generate_seqs_and_conf(NS=NS,base_dir = rel_dir,results_dir='data')

sim_fc = eval(parse(text=paste('sim_',shortname_from,sep='')))
simulatedData = sim_fc(NS,drawnPars,list_names,seq_data)
TN = dim(simulatedData[[1]])[2]
simulatedDataFull = list(N=NS, Tn=TN,N_obs_VAS=dim(seq_data$RatingConf)[2],N_obs_CP=dim(seq_data$PredConf)[2],
                         noxInputTrans=seq_data$noxInputTrans,painRating=simulatedData[[2]],painPred=simulatedData[[3]],
                         RatingConf=seq_data$RatingConf,PredConf=seq_data$PredConf,TrialType=matrix(0,NS,TN),
                         TrialTypeAllWMissing=seq_data$TrialTypeAllWMissing,PainValsAll=matrix(0,NS,TN),
                         PercIndexArrayLong=matrix(0,NS,TN),PredIndexArrayLong=matrix(0,NS,TN),
                         PercIndexArray=seq_data$PercIndexArray,PredIndexArray=seq_data$PredIndexArray,
                         IndexMissingAll=seq_data$IndexMissingAll)
print(paste('Simulated data from ', shortname_from,sep=''))

# Setup and run stan fit -------------
fitted_params = vector(mode = "list", length = length(model_pars_list_with))
fitted_params_summary = vector(mode = "list", length = length(model_pars_list_with))
par_corrs = array(0,length(model_pars_list_with))

print(paste('Started fitting data for',list_names,'with ',shortname_with))
print(paste('Stan specs are: ',specfname,sep=''))

HBA_ind_fit = stan(file=model_name_with,data=simulatedDataFull,cores=stan_cores,verbose=FALSE,save_warmup=FALSE,pars=c('lp_'),include=FALSE,
                   iter=iters,warmup=warmups, control = list(adapt_delta=adapt_delta_value,max_treedepth=max_treedepth_value))

warnings()
print(paste('Done fitting data for',list_names,'with',shortname_with))

# Model recovery -------------
mrec_colname = c('data_from', 'fit_with','looic','waic','group','s_iter','seed','Nsubj','treedepth','energy','div')
options(warn=-1)
sink('nul')
mr_looic = loo(extract_log_lik(HBA_ind_fit))$estimates['looic',1]
mr_waic = waic(extract_log_lik(HBA_ind_fit))$estimates['waic',1]
mrec_df = data.frame(cbind(shortname_from,shortname_with,mr_looic,mr_waic,list_names,s_iter,seed,NS,check_treedepth(HBA_ind_fit),check_energy(HBA_ind_fit),check_div(HBA_ind_fit)))
options(warn=0)
sink()
names(mrec_df) = mrec_colname
try(if (save_outputs){
  mrec_subfname = paste('from_',m_from_sh, '_with_', m_with_sh,'_','s',s_iter,'_g',list_names,'_N',NS,'_',specfname,script_stamp,'.csv',sep='')
  write.csv(mrec_df,file=paste(outputdir_mrec,"/mrec",'_',mrec_subfname,sep=''),row.names=FALSE)
  print(paste('saved model recovery info for ',list_names,'from ',shortname_from,' with ',shortname_with,' s_iter: ',s_iter))
})

# Parameter recovery run -------------
if (m_from_sh==m_with_sh){
  # Prepare individual paramter and group summary tables -------------
  for (j in 1:length(model_pars_list_with)){
    fitted_params[[j]]<-summary(HBA_ind_fit,pars=model_pars_list_with[j])$summary[,'mean']
  }
  
  for (j in 1:length(model_pars_list_with)){
    fitted_params_summary[[j]]=c(mean(fitted_params[[j]]),sd(fitted_params[[j]]))
  }
  names(fitted_params_summary)<-model_pars_list_with
  
  fitted_params_summary_df = data.frame(fitted_params_summary)
  names(fitted_params_summary_df)<-model_pars_list_with
  
  fitted_params_df<-data.frame(fitted_params)
  names(fitted_params_df)<-model_pars_list_with
  rownames(fitted_params_df)<-NULL
  
  # Save p recovery extra to CSV -------------
  try(if (save_outputs){
    prec_extra_subfname = paste(m_with_sh,'_','s',s_iter,'_g',list_names,'_N',NS,'_',specfname,script_stamp,'.csv',sep='')
    write.csv(drawnPars,file=paste(outputdir_prec_extra,"/sim_params_",prec_extra_subfname,sep=''),row.names=FALSE)
    write.csv(fitted_params_summary_df,file=paste(outputdir_prec_extra,"/HBA_summary_",prec_extra_subfname,sep=''),row.names=FALSE)
    write.csv(fitted_params_df,file=paste(outputdir_prec_extra,"/HBA_ind_fit",'_',prec_extra_subfname,sep=''),row.names=FALSE)
    print(paste('saved (drawn) parameters and parameter summary for',list_names,'with',shortname_with,' s_iter: ',s_iter))
  })
  
  # Correlation for parameter recovery + save file-------------
  for (j in 1:length(model_pars_list_with)){
    par_corrs[j] = cor(drawnPars[,j],fitted_params_df[,j])
  }
  # print(paste('Correlation for '))
  par_corrs = transpose(data.frame(par_corrs))
  names(par_corrs)=model_pars_list_with
  
  prec_df = data.frame(cbind(shortname_with,list_names,s_iter,seed,NS,par_corrs))
  prec_colname = c('model','group','s_iter','seed','Nsubj',model_pars_list_with)
  names(prec_df)=prec_colname
  
  try(if (save_outputs){
    prec_subfname = paste(m_with_sh,'_','s',s_iter,'_g',list_names,'_N',NS,'_',specfname,script_stamp,'.csv',sep='')
    write.csv(prec_df,file=paste(outputdir_prec,"/prec",'_',prec_subfname,sep=''),row.names=FALSE)
    print(paste('saved parameter recovery info for ',list_names,'with',shortname_with,' s_iter: ',s_iter))
  })
  
}

warnings()