if (Sys.info()["sysname"]=="Darwin"){
  # setwd('/Users/onyskj/Library/CloudStorage/OneDrive-UniversityofCambridge/CamNoxLab - OneDrive/tsl_paper/pain_control_git/paper_submission/experiment1')
  setwd('/Users/onyskj/Library/CloudStorage/OneDrive-UniversityofCambridge/CamNoxLab - OneDrive/tsl_paper/elife_revision/')
}else{
  setwd('/rds/project/rds-3IOyKgCQu4I/tsl_paper')
}
# source('stan_utility.R')
# library(rstan)
# library(shinystan)
# library(loo)
# library(hBayesDM)
library(data.table)
library(plyr)
library(latex2exp)
library("ggsci")
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyverse)
rel_dir = 'model_recovery/'
results_dir = 'model_recovery/results_revision/'
results_dir2 = 'model_recovery/results_pre_revision/'
# results_dir2 = 'model_recovery/results_continued/'
save_analysis='model_recovery/'
# results_dir = 'model_recovery/results/'
results_subdir = ''
specific_model ='.*'
letters_list = toupper(letters)

# param recovevery permutations setings
runPerm = FALSE
savePerm = FALSE
perms = 500

# Get mrec DFs and find bad runs based on divergences and energy measure -----
csv_files_mrec_1 <- list.files(path = paste(results_dir,'model_rec/',sep=''), recursive = TRUE, pattern = "\\.csv$", full.names = TRUE)
csv_files_mrec_2 <- list.files(path = paste(results_dir2,'model_rec/',sep=''), recursive = TRUE, pattern = "\\.csv$", full.names = TRUE)
#combined two sep runs
csv_files_mrec<- c(csv_files_mrec_1,csv_files_mrec_2)
combined_mrec <- data.frame()
combined_exclude <-data.frame()
for (csv_file in csv_files_mrec) {
  data <- read.csv(csv_file, header = TRUE)
  if (!is_logical((data$energy))){
    data$energy = TRUE
    combined_exclude <-bind_rows(combined_exclude,data)
  }else{
    if (data$div>1 || data$energy){
      combined_exclude <-bind_rows(combined_exclude,data)
    }else{
      combined_mrec <-bind_rows(combined_mrec,data)
    }
  }
}
combined_exclude$s_iter = as.integer(combined_exclude$s_iter)
combined_exclude_prec = combined_exclude[combined_exclude$data_from==combined_exclude$fit_with,]
combined_mrec$s_iter = as.integer(combined_mrec$s_iter)
# combined_mrec %>% count(data_from,fit_with)

# Get Paramter recovery DFs and remove bad runs (combine two sep runs)------------- 
p_rec_m_dirs_1 = dir(paste(results_dir,results_subdir,'param_rec',sep=''),full.names = T,pattern = paste('m',specific_model,sep=''))
p_rec_m_dirs_2 = dir(paste(results_dir2,results_subdir,'param_rec',sep=''),full.names = T,pattern = paste('m',specific_model,sep=''))
p_rec_m_dirs_list = list(p_rec_m_dirs_1,p_rec_m_dirs_2)
p_rec_df_list = vector(mode = "list", length = length(p_rec_m_dirs_1))
for (m in 1:length(p_rec_df_list)){
    p_rec_df_list[[m]] = data.frame()
}
p_rec_df_long = data.frame()
for (p_rec_m_dirs in p_rec_m_dirs_list){
  for (m in 1:length(p_rec_m_dirs)){
    # print(m)
    p_rec_files = dir(p_rec_m_dirs[m],full.names = T,pattern = 'prec_.*')
    if (length(p_rec_files)>0){
      for (f in 1:length(p_rec_files)){
  
        p_rec_temp = (read.table(file = p_rec_files[f], header = TRUE,sep=','))
        s_iter_extract = p_rec_temp$s_iter
        p_rec_temp = p_rec_temp[1,c(1,2,6:dim(p_rec_temp)[2])]
        
        pars_subset = p_rec_temp[3:dim(p_rec_temp)[2]]
        
        # include only good ones
        tmp_row <- data.frame(cbind(rep(p_rec_temp$model,length(pars_subset)),rep(p_rec_temp$group,length(pars_subset)),colnames(pars_subset),as.numeric(pars_subset),rep(s_iter_extract,length(pars_subset))))[1,]
        colnames(tmp_row)<-c('model','group','parameter','value','s_iter')
        tmp_row$s_iter = as.integer(tmp_row$s_iter)
  
        # is_good = all(!((tmp_row$model == combined_exclude$fit_with) & (tmp_row$model == combined_exclude$data_from) & (tmp_row$s_iter == combined_exclude$s_iter)))
        is_good = all(!((tmp_row$model == combined_exclude_prec$fit_with) & (tmp_row$s_iter == combined_exclude_prec$s_iter)))
        if (is_good){
          p_rec_df_long = rbind(p_rec_df_long,cbind(rep(p_rec_temp$model,length(pars_subset)),rep(p_rec_temp$group,length(pars_subset)),
                                                    colnames(pars_subset),as.numeric(pars_subset),rep(s_iter_extract,length(pars_subset))))
           
          p_rec_df_list[[m]] = rbind(p_rec_df_list[[m]],p_rec_temp)
        }
        else{
          print(tmp_row)
          print(is_good)
        }
        
      }
    }
  }
}
colnames(p_rec_df_long)<-c('model','group','parameter','value','s_iter')
p_rec_df_long$s_iter = as.integer(p_rec_df_long$s_iter)
# p_rec_df_long %>% count(model,parameter)

#get model names long and short
model_names_long = unique(p_rec_df_long$model)
model_names_short = gsub('_norm','',gsub('sliced2_halfStudentT_reparam_confEs_','',model_names_long))
model_names_short = gsub('RL','eRL',model_names_short); model_names_short = gsub('KF','eKF',model_names_short)
model_names_short = gsub('eRL_np','RL',model_names_short); model_names_short = gsub('eKF_np','KF',model_names_short)
model_names_short = gsub('_C','_Random',model_names_short)
model_names_short = gsub('m.*_','',model_names_short)


## Calculate mean and sd correlation for each parameter (recovery) ---------
groups = c('All')
model_names = unique(p_rec_df_long$model)
group_names = unique(p_rec_df_long$group)
p_rec_final_list_All = vector(mode = "list", length = length(p_rec_m_dirs_1))
for (m in 1:length(p_rec_m_dirs)){
  p_rec_final_list_All[[m]] = data.frame()
  g=1
  if (groups[g] %in% group_names & nrow(p_rec_df_list[[m]])>1){
  
    p_rec_model_df = p_rec_df_list[[m]][(p_rec_df_list[[m]]$group==groups[g]),]
    mean_prec = round(sapply(p_rec_model_df[,3:dim(p_rec_model_df)[2]],mean),3)
    sd_prec = round(sapply(p_rec_model_df[,3:dim(p_rec_model_df)[2]],sd),3)
    meanSd_prec=paste(as.vector(mean_prec),' (',as.vector(sd_prec),')',sep='')
    p_rec_final_list_All[[m]] = t(data.frame(meanSd_prec))
    colnames(p_rec_final_list_All[[m]])=names(mean_prec)
  }
}

# Save param recovery tables for each model -------
prec_save_dir = paste(save_analysis,'analysed/prec/',sep='')
if(!file.exists(prec_save_dir)){dir.create(prec_save_dir,showWarnings = TRUE,recursive = TRUE)}
for (m in 1:length(p_rec_m_dirs)){
  if (groups[1] %in% group_names & nrow(p_rec_final_list_All[[m]])>0){
    rownames(p_rec_final_list_All[[m]])='mean (sd)'
    # write.csv(p_rec_final_list_All[[m]],paste(prec_save_dir,'p_rec_All_m',m,'.csv',sep=''))
  }
}


# DF with individual parameter vaulues (simulation vs fit) for each model -------
# p_rec_extra_m_dirs = dir(paste(results_dir,results_subdir,'param_rec_extra',sep=''),full.names = T,pattern = paste('m',specific_model,sep=''))
csv_files_1 <- list.files(path = paste(results_dir,'param_rec_extra/',sep=''), recursive = TRUE, pattern = "(.*sim.*|.*fit.*)\\.csv$", full.names = TRUE)
csv_files_2 <- list.files(path = paste(results_dir2,'param_rec_extra/',sep=''), recursive = TRUE, pattern = "(.*sim.*|.*fit.*)\\.csv$", full.names = TRUE)
csv_files <- c(csv_files_1,csv_files_2)
combined_data <- data.frame()
for (csv_file in csv_files) {
  file_parts <- str_split(csv_file, "/")[[1]]
  model_name <- str_split(file_parts[5],"_")[[1]][1]
  s_iter <- as.numeric(sub('[A-z]','',file_parts[6]))
  length(s_iter)
  csv_name<- file_parts[7]
  source_name<- if (grepl("fit",csv_name)) 'fit' else if ((grepl("sim",csv_name))) "sim" else "summary"
  if (source_name !='summary'){
    data <- read.csv(csv_file, header = TRUE)
    data_long <- data %>%
      rownames_to_column(var = "sub") %>%
      pivot_longer(cols = -sub, names_to = "parameter", values_to = "value")%>%
      mutate(source = source_name, s_iter = s_iter)
  }
  data_long['model']=model_name
  combined_data <- bind_rows(combined_data, data_long)
}
# get model_names
model_names = unique(combined_data$model)
combined_data$model <-mapvalues(combined_data$model,from=c('m1','m2','m3','m4','m5'), to=model_names)
combined_data$s_iter = as.integer(combined_data$s_iter)
# exclude bad runs
if (dim(combined_exclude_prec)[1]>0){
  combined_data = combined_data[!(combined_data$model==combined_exclude_prec$fit_with & combined_data$s_iter==combined_exclude_prec$s_iter),]
}
#remove an outliner
combined_data = subset(combined_data, !(s_iter==41 & model=='m3' & parameter=='xi'))
combined_data = subset(combined_data, !(s_iter==70 & model=='m3' & parameter=='xi'))
combined_data = subset(combined_data, !(s_iter==106 & model=='m1' & parameter=='eta'))

# Plot parameter recovyer scatter plots -------
# Create a paramter name to latex format map
par_ind_latex_names = list(m1=c("$\\alpha","$C","$E^0","$\\xi","$\\gamma"),
                           m2=c("$\\alpha","$C","$E^0","$\\xi"),
                           m3=c("C","$E^0","$\\epsilon","$v","$s","$w^0","$\\xi"),
                           m4=c("C","$E^0","$v","$s","$w^0","$\\xi"),
                           m5=c("$R","$C","$\\xi"))

parameter_names <- list(m1=c("alpha", "cs", "E0", "eta", "gam"),
                        m2=c("alpha", "cs", "E0", "eta"),
                        m3=c("cs", "E0", "eps", "eta", "psi", "w0", "xi"),
                        m4=c("cs", "E0", "eta", "psi", "w0", "xi"),
                        m5=c("C", "cs", "eta"))
latex_df <- bind_rows(lapply(names(par_ind_latex_names), function(model) {
  data.frame(model = model, latex_name = unlist(par_ind_latex_names[model]), stringsAsFactors = FALSE)
}))
param_df <- bind_rows(lapply(names(parameter_names), function(model) {
  data.frame(model = model, param_name = unlist(parameter_names[model]), stringsAsFactors = FALSE)
}))
param_name_map <- merge(latex_df, param_df, by = "model")
param_name_map <- param_name_map[, c("model", "param_name", "latex_name")]

# m1 - alpha(lr), cs (conf), E0 (init pain), eta (beh resp noise), gamma (perceptiual weigh)
# m2 - same without gammaj
# m3 - cs (conf), E0 (init pain), eps (subj noise), eta (v; volat.) psi(s; stoch.), w0(init uncertainty), xi (beh resp noise)
# m4 same without epsiolon
# m5 - C - R (exp pain int); cs - C (conf)  ;eta-xi (resp noise)
combined_list <- list()
for (model_name in names(parameter_names)) {
  param_name_map <- data.frame(model = model_name,
                              param_name = parameter_names[[model_name]],
                              latex_name = par_ind_latex_names[[model_name]])
  combined_list[[model_name]] <- param_name_map
}
param_name_map <- do.call(rbind, combined_list)
row.names(param_name_map) <- NULL


wr=23*0.75; hr=5.5*0.75
font_size = 30
font_size_leg=28
for (m in 1:length(model_names)){
  param_names = unique(combined_data[combined_data$model==model_names[m],'parameter'])
  pl = vector(mode = "list", length = length(model_names))
  for (p in 1:length(param_names)){
    df_to_plot = combined_data[combined_data$model==model_names[m] & combined_data$parameter==param_names[p],]
    df_to_plot$id = paste(df_to_plot$sub,df_to_plot$s_iter,sep='_')
    df_to_plot = subset(df_to_plot,select=-c(sub,s_iter))
    df_to_plot = reshape(df_to_plot,idvar=c("id","parameter","model"),timevar = "source",direction="wide")
    p_rec_final_list_All[[m]][param_names[p]]
    
    par_name_latex = param_name_map[(param_name_map$model == substr(model_names[m],1,2)) & (param_name_map$param_name == param_names[p]),'latex_name']
    corr_mean_std_val = p_rec_final_list_All[[m]][colnames(p_rec_final_list_All[[m]])==param_names[p]]
    p_rec_plot_title = unname(TeX((paste(par_name_latex, '$: $', corr_mean_std_val,sep=''))))
    ylab_title = ''
    if (p==1){
      ylab_title = "Value - Fit"
    }
    
    axmax= max(df_to_plot[c('value.fit','value.sim')])
    axmin=min(df_to_plot[c('value.fit','value.sim')])
    
    pl[[p]]<-ggplot(data=df_to_plot, aes(x=value.sim,y=value.fit))+
      geom_point()+
      ylab(ylab_title)+
      xlab("Value - Simulation")+
      labs(subtitle = letters_list[p],title=p_rec_plot_title)+
      coord_cartesian(ylim = c(axmin, axmax),xlim=c(axmin,axmax))+
      theme(axis.text.y=element_text(size=font_size))+
      theme(axis.text.x=element_text(size=font_size,angle=0))+
      theme(plot.subtitle = element_text(size = font_size,face='bold',hjust=-0.15, vjust=0.45))+
      theme(text=element_text(size=font_size))+
      theme(text = element_text(family = "Arial"))+
      theme(plot.title=element_text(hjust = 0.5,vjust=-2.75))+
      scale_fill_discrete(labels=unname(TeX(par_ind_latex_names[[m]])))+
      scale_color_discrete(labels=unname(TeX(par_ind_latex_names[[m]])))
  }
  eps_ar_plot=grid.arrange(grobs=pl,nrow=1,guide_legend(nrow=1, byrow=TRUE), top = textGrob(model_names_short[m], gp=gpar(fontsize=font_size*1.2),vjust=1))
  # ggplot2::ggsave(paste(save_analysis,'analysed/scatter_plots_model_',m,'.png',sep=''),plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
  # ggplot2::ggsave(paste(save_analysis,'analysed/scatter_plots_model_',m,'.eps',sep=''),device = cairo_ps, plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
}




# Model recovery -------------
# Calculate confusion matrix
model_names = unique(combined_mrec$data_from)
sIters = unique(combined_mrec$s_iter)
count_mrec_All = array(0,dim=c(length(model_names),length(model_names)))
s_count_All = rep(0,length(model_names))
for (m_from in 1:length(model_names)){
  g=1
  if (groups[g] %in% group_names){
    for (s in 1:length(sIters)){
      sub_df = combined_mrec[combined_mrec$group==groups[g] & combined_mrec$data_from==model_names[m_from] & combined_mrec$s_iter==sIters[s],]
      if (nrow(sub_df)==5){
        s_count_All[m_from] = s_count_All[m_from]+1 
        min_looic_df = sub_df[sub_df$looic==min(sub_df$looic),]
        which(model_names==min_looic_df$data_from)
        which(model_names==min_looic_df$fit_with)
        count_mrec_All[which(model_names==min_looic_df$data_from),which(model_names==min_looic_df$fit_with)]=count_mrec_All[which(model_names==min_looic_df$data_from),which(model_names==min_looic_df$fit_with)]+1
      }else {
        print(nrow(sub_df))
      }
    }
  }
}

if (groups[1] %in% group_names){
  count_mrec_All = round(t(apply(count_mrec_All,1,function(x) x/sum(x) )),3)
  count_mrec_All = data.frame(count_mrec_All)
  names(count_mrec_All)=model_names_short
  rownames(count_mrec_All)=model_names_short
  # names(count_mrec_All)=model_names
  # rownames(count_mrec_All)=model_names
  # write.csv(count_mrec_All,paste(save_analysis,'analysed/m_rec_All.csv',sep=''))
}



# Plot parameter recovery error (permuted) ------------- 
# it_pts = c(5,10,15,20,25)
# it_pts = c(seq(2,25,2),25)
# it_pts=c(seq(2,50,1))
# sIters = sort(unique(p_rec_df_long$s_iter))
# combined_exclude$s_iter = as.integer(combined_exclude$s_iter)
# p_rec_df_long=p_rec_df_long[!((p_rec_df_long$model %in% combined_exclude$fit_with) & (p_rec_df_long$model %in% combined_exclude$data_from) & (p_rec_df_long$s_iter %in% combined_exclude$s_iter)),]

p_rec_df_long$s_iter = as.integer(p_rec_df_long$s_iter)
groups = c('All')
model_names = unique(p_rec_df_long$model)
group_names = unique(p_rec_df_long$group)
# sIters = sort(unique(p_rec_df_long$s_iter))
errors_list = vector(mode = "list", length = length(group_names))
names(errors_list) = group_names
if (runPerm){
  for (g in 1:length(group_names)){
    errors_list[[g]] = vector(mode = "list", length = length(model_names))
    names(errors_list[[g]]) =model_names 
    for (m in 1:length(model_names)){
      errors_list[[g]][[m]] = vector(mode = "list", length = 3)
      names(errors_list[[g]][[m]]) = c('mean','sd')
      temp_perm_df = data.frame()
      temp_perm_df_pererr = data.frame()
      temp_perm_df_all = data.frame()
      tmp_m_df = p_rec_df_long[p_rec_df_long$model==model_names[m],]
      m_its = sort(unique(tmp_m_df$s_iter))
      for (s in 2:length(m_its)){
        temp_df = data.frame()
        for (p in 1:perms){
          seed = sample(1:99999,size=1)
          # print(seed)
          set.seed(seed)
          # order_s=sample.int(n_its, n_its)
          
          
          # sample s number of simulations from the available list of simulations for that model
          order_s = sample(m_its,s)
  
          temp_sim = p_rec_df_long[p_rec_df_long$model==model_names[m] & p_rec_df_long$group==group_names[g] & p_rec_df_long$s_iter %in% order_s,]
          # temp_sim = p_rec_df_long[p_rec_df_long$model==model_names[m] & p_rec_df_long$group==group_names[g] & p_rec_df_long$s_iter %in% order_s[1:(it_pts[s])],]
          sd_prec=ddply(temp_sim,~parameter,summarise,std=sd(value))
          temp_df=rbind(temp_df,sd_prec)
          
        }
        temp_perm_df=rbind(temp_perm_df,cbind(ddply(temp_df,~parameter,summarise,mean=round(mean(std),3)),rep(s,length(unique(temp_df$parameter)))))
        temp_perm_df_pererr=rbind(temp_perm_df_pererr,
                                  cbind(ddply(temp_df,~parameter,summarise,sd=round(sd(std),3)),rep(s,length(unique(temp_df$parameter))))
                                  )
      }
      temp_perm_df_all = rbind(temp_perm_df_all,cbind(temp_perm_df[,1:2],temp_perm_df_pererr[,2:3]))
      names(temp_perm_df)=c('parameter','mean','pt')
      names(temp_perm_df_pererr)=c('parameter','sd','pt')
      names(temp_perm_df_all)=c('parameter','mean','sd','pt')
      errors_list[[g]][[m]][[1]] = temp_perm_df
      errors_list[[g]][[m]][[2]] = temp_perm_df_pererr
      errors_list[[g]][[m]][[3]] = temp_perm_df_all
    }
  }

  joint_pars_m_g_df = data.frame()
  for (g in 1:length(group_names)){
    for (m in 1:length(model_names)){
      tmp_df = errors_list[[g]][[m]][[3]]
      tmp_df$group=group_names[g]
      tmp_df$model=model_names_short[m]
      joint_pars_m_g_df = rbind(joint_pars_m_g_df,tmp_df)
    }
  }
} else{
  joint_pars_m_g_df <- read.csv(paste(save_analysis,'analysed/pr_error_nperm_',perms,'.csv',sep=''), header = TRUE)
}
if (savePerm && runPerm){
  write.csv(joint_pars_m_g_df,paste(save_analysis,'analysed/pr_error_nperm_',perms,'.csv',sep=''), row.names=FALSE)
}



plot_types=c('SEM','SD')
perms_types=c(perms,perms/perms)
for (plt in 1:length(plot_types)){
  pl = vector(mode = "list", length = length(model_names))
  wr=23*0.75; hr=5.5*0.75
  letters_list = toupper(letters)
  font_size = 30
  font_size_leg=28
  # xlimcalc = max(joint_pars_m_g_df$pt)+(ceiling(max(joint_pars_m_g_df$pt)/10)*10-max(joint_pars_m_g_df$pt))
  # xlimcalc = round(max(joint_pars_m_g_df$pt)*1.05)
  xlimcalc = 110
  for (m in 1:length(model_names)){
    pl[[m]] <- ggplot(data=joint_pars_m_g_df[joint_pars_m_g_df$model==model_names_short[m],],aes(x=pt,y=mean,ymin=mean-sd/sqrt(perms_types[plt]),ymax=mean+sd/sqrt(perms_types[plt]),fill=parameter,col=parameter))+
      geom_line(size=2)+
      geom_ribbon(alpha=0.35,col=NA)+
      ylab(if(m==1){'Mean error'}else{''})+
      xlab('N of simulations')+
      labs(fill='',col='')+
      theme(legend.position = 'bottom')+
      # theme(legend.position = c(-0.1,-0.5))+
      # guides(colour = guide_legend(nrow = 1,label.hjust=-2.5))+
      # guides(fill = guide_legend(nrow = 1,label.hjust=-2.5))+
      guides(colour = guide_legend(nrow = 1,label.hjust=-2))+
      guides(fill = guide_legend(nrow = 1))+
      labs(subtitle = letters_list[m],title=gsub("m.*_","",model_names_short[m]))+
      coord_cartesian(ylim = c(0, 0.3),xlim=c(0,xlimcalc))+
      # coord_cartesian(ylim = c(0, 0.3))+
      theme(axis.text.y=element_text(size=font_size))+
      theme(axis.text.x=element_text(size=font_size*0.7,angle=0))+
      theme(plot.subtitle = element_text(size = font_size,face='bold',hjust=-0.05))+
      theme(text=element_text(size=font_size),
            legend.text = element_text(size = font_size_leg),
            legend.title = element_text(size = font_size_leg))+
      theme(text = element_text(family = "Arial"))+
      theme(plot.title=element_text(hjust = 0.5,vjust=-2.75))+
      # scale_color_d3(labels=unname(TeX(par_ind_latex_names[[m]])))
      scale_fill_manual(labels=unname(TeX(par_ind_latex_names[[m]])), values = c('#1E88E5','#D81B60','#FFC107',"#A9D0C9", '#A78036','#46C3D9','#6A1A8A'))+
      scale_color_manual(labels=unname(TeX(par_ind_latex_names[[m]])), values = c('#1E88E5','#D81B60','#FFC107',"#A9D0C9", '#A78036','#46C3D9','#6A1A8A'))+
      # scale_linetype_manual(labels=unname(TeX(par_ind_latex_names[[m]])), values=c(1,2,3,4,5,6,1))+
      scale_x_continuous(breaks=unique(c(1,floor(seq(floor(sqrt(xlimcalc)),xlimcalc,sqrt(xlimcalc))))))+
      if(m>1){scale_y_continuous(breaks=NULL)}else{NULL}
    # scale_fill_d3(labels=unname(TeX(par_ind_latex_names[[m]])))
    
  }
  eps_ar_plot=grid.arrange(grobs=pl,nrow=1,guide_legend(nrow=1, byrow=TRUE))
  ggplot2::ggsave(paste(save_analysis,'analysed/error_plot_prec_',plot_types[plt],'.png',sep=''),plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
  ggplot2::ggsave(paste(save_analysis,'analysed/error_plot_prec_',plot_types[plt],'.eps',sep=''),device=cairo_ps,plot=eps_ar_plot,dpi=300,width=614*wr,height=376*hr,units="px")
}
