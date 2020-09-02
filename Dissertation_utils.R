
library(data.table)
library(ggplot2)
library(evd)
library(doMC)
library(doRNG)
library(EnvStats)

.NUM_CORES=40
registerDoMC(cores=.NUM_CORES)
#registerDoRNG()

get_a_d_stat<-function(Y,CDF=pgev, ...){
  Y=Y[order(Y)]
  F_i=CDF(Y,...)
  m=length(Y)
  i=1:m
  inner_sum=((2*i)-1)*(log(F_i)+log(1-rev(F_i)))
  return(-m-sum(inner_sum)/m)
}

inrange<-function(x,lower,upper){
  (x<=upper)&(x>=lower)
}

get_m_a_d_u_stat<-function(Y,CDF=pgev, ...){
  Y=Y[order(Y)]
  F_i=CDF(Y,...)
  m=length(Y)
  i=1:m
  inner_sum=(2-((2*i)-1)/m)*log(1-F_i)
  return(m/2-2*sum(F_i)-sum(inner_sum))
}

get_m_a_d_l_stat<-function(Y,CDF=pgev, ...){
  Y=Y[order(Y)]
  F_i=CDF(Y,...)
  m=length(Y)
  i=1:m
  inner_sum=(((2*i)-1)/m)*log(F_i)
  t1=-((3*m)/2)
  t2=2*sum(F_i)
  return(t1+t2-sum(inner_sum))
}


get_gev_return_value_and_std.err<-function(p,fit){
  loc=fit$estimate[['loc']]
  scale=fit$estimate[['scale']]
  shape=fit$estimate[['shape']]
  y=-log(1-p)
  
  return=unname(qgev(p=p,loc=loc,
       scale=scale,
       shape=shape,
       lower.tail=F))
  
  delta_z<-unname(c(1,
             -(1/shape)*(1-y^(-shape)),
             (scale/shape^2)*(1-y^(-shape))-(scale/shape)*(y^(-shape))*log(y)))
  
  se= sqrt(c(t(delta_z))%*%unname(fit$var.cov)%*%delta_z)
  return(c(estimate=return,se=se))
}

get_gpd_return_value_and_std.err<-function(p,fit){
  thresh=fit$threshold
  scale=fit$estimate[['scale']]
  shape=fit$estimate[['shape']]
  p_u=fit$pat
  n=length(fit$data)
  var_p_u=p_u*(1-p_u)/n
  
  w1=matrix(data = c(var_p_u,0,0), nrow = 1, ncol = 3)
  w2=cbind(matrix(0,nrow=2,ncol=1),fit$var.cov)
  w=rbind(w1,w2)
  
  if(abs(shape)>0.000001){
    return= thresh + (scale/shape) * (((p_u/p)^shape)-1)
  }else{
    return = thresh + scale*log(p_u/p) 
  }

  
  delta_z<-unname(c(scale*(1/p)^shape*p_u^(shape-1),
                    (((p_u/p)^shape)-1)/shape,
                    -(scale*((p_u/p)^shape-1)/(shape^2))+(scale*(p_u/p)^shape*log(p_u/p))/shape
  ))
  
  se= sqrt(c(t(delta_z))%*%unname(w)%*%delta_z)
  return(c(estimate=unname(return),se=se))
}

return_estimated_stats<-function(sample,
                                 CDF_function,
                                 fit_function,
                                 return_value_function = get_gev_return_value_and_std.err,
                                 identifier,
                                 loc=loc,
                                 shape=shape,
                                 scale=scale,
                                 group_size=group_size,
                                 replicates=replicates){
  
  tryCatch({fit<-fit_function(sample,std.err = TRUE, corr = T)},
           error = function(e) {
            print('NO_ERRS')
            fit<-fit_function(sample,std.err = FALSE, corr = FALSE)
           }
           )
  if(!exists('fit')){
    print('NO FIT')
    return(NULL)
  }
  if(fit$estimate[['scale']]<=0){
    print("FAIL SCALE")
    return(NULL)
  }else{
  estimated_10_return=return_value_function(1/10,
                                               loc=fit$estimate['loc'],
                                               scale=fit$estimate['scale'],
                                               shape=fit$estimate['shape'])
  estimated_100_return=return_value_function(1/100,
                                             loc=fit$estimate['loc'],
                                             scale=fit$estimate['scale'],
                                             shape=fit$estimate['shape'])
  estimated_1000_return=return_value_function(1/1000,
                                              loc=fit$estimate['loc'],
                                              scale=fit$estimate['scale'],
                                              shape=fit$estimate['shape']) 
    
  return(data.table(identifier=identifier,
                    loc=loc,
                    shape=shape,
                    scale=scale,
                    group_size=group_size,
                    replicates=replicates,
                    deviance=fit$deviance,
               estimated_loc=fit$estimate[['loc']],
               estimated_scale=fit$estimate[['scale']],
               estimated_shape=fit$estimate[['shape']],
               estimated_10_return=estimated_10_return[['estimate']],
               estimated_100_return=estimated_100_return[['estimate']],
               estimated_1000_return=estimated_1000_return[['estimate']],
               n_over_10_return=sum(sample>estimated_10_return[['estimate']]),
               n_over_100_return=sum(sample>estimated_100_return[['estimate']]),
               n_over_1000_return=sum(sample>estimated_1000_return[['estimate']]),
               loc_err=fit$std.err[['loc']],
               scale_err=fit$std.err[['scale']],
               shape_err=fit$std.err[['shape']],
               return_10_err=estimated_10_return[['se']],
               return_100_err=estimated_100_return[['se']],
               return_1000_err=estimated_1000_return[['se']],
               AD=get_a_d_stat(sample,
                               CDF_function,
                               loc=fit$estimate[['loc']],
                               scale=fit$estimate[['scale']],
                               shape=fit$estimate[['shape']]),
               UMAD = get_m_a_d_u_stat(sample,
                                       CDF_function,
                                       loc=fit$estimate[['loc']],
                                       scale=fit$estimate[['scale']],
                                       shape=fit$estimate[['shape']]),
               fit=list(fit)))
  }
}


parallel_stat_acquisition<-function(data,
                                    loc=loc,
                                    shape=shape,
                                    scale=scale,
                                    group_size=group_size,
                                    replicates=replicates,
                    CDF_function=pgev,
                    fit_function=fgev,
                    return_value_function = get_gev_return_value_and_std.err){


  
  setkey(data, "group")
  
  finalRowOrderMatters = FALSE # FALSE can be faster
  return(foreach(x=unique(data[["group"]]), .combine="rbind", .inorder=finalRowOrderMatters) %dopar% 
    data[.(x) ,return_estimated_stats(sample,
                                      CDF_function=CDF_function,
                                      fit_function=fit_function,.(x),
                                      return_value_function=return_value_function,
                                      loc=loc,
                                      shape=shape,
                                      scale=scale,
                                      group_size=group_size,
                                      replicates=replicates)])
  
  
}

get_ad_test_bootstrap<-function(group_size,
                                number_of_groups,
                                loc,
                                shape,
                                scale,
                                CDF_function=pgev,
                                fit_function=fgev,
                                generation_function=rgev,
                                return_value_function = get_gev_return_value){
  print(paste(Sys.time(),': Generating data'))
  data<-data.table(sample=generation_function(number_of_groups*group_size,
                               loc=loc,
                               scale=scale,
                               shape=shape),
                   group=rep(1:number_of_groups,group_size))
  print(paste(Sys.time(),': getting stats'))
  parallel_stat_acquisition(data,
                                    loc=loc,
                                    shape=shape,
                                    scale=scale,
                                    group_size=group_size,
                                    replicates=number_of_groups,
                            CDF_function,
                            fit_function,
                            return_value_function)
}


get_bootstrap_groups<-function(vector_to_bootstrap,
                               number_of_groups,
                               loc,
                               shape,
                               scale,
                               CDF_function=pgev,
                               fit_function=fgev,
                               return_value_function = get_gev_return_value_and_std.err
){
  print(paste(Sys.time(),': Generating data'))
  data<-lapply(1:number_of_groups,function(x){data.table(sample=sample(vector_to_bootstrap,replace = T),
                                                   group=x)})
  data<-rbindlist(data)
  print(paste(Sys.time(),': getting stats'))
  parallel_stat_acquisition(data,
                            loc=loc,
                            shape=shape,
                            scale=scale,
                            group_size=length(data),
                            replicates=number_of_groups,
                            CDF_function,
                            fit_function,
                            return_value_function = return_value_function)
}

generate_grouped_data<-function(number_of_groups,
                                group_size,
                                generation_function=rnorm,
                                num_cores=1,
                                ...){
  n=number_of_groups*group_size
  if(num_cores==1){
    rand=generation_function(n,...)    
  }else{
    rand<-mclapply(1:num_cores,function(x){
      generation_function(n/num_cores,...)
    },mc.cores=num_cores)
    rand<-unlist(rand)
  }

  data.table(rand=rand,group=rep(1:number_of_groups,group_size))
}

return_specified_stats<-function(data,
                                 CDF_function,
                                 return_value_function = qgev,
                                 loc=loc,
                                 shape=shape,
                                 scale=scale,
                                 group_size=group_size,
                                 replicates=replicates){
  
  estimated_100_return=return_value_function(1/100,
                                             loc,
                                             scale,
                                             shape,
                                             lower.tail=FALSE)
  estimated_1000_return=return_value_function(1/1000,
                                              loc,
                                              scale,
                                              shape,
                                              lower.tail=FALSE) 
    
    return(data[,list(identifier=identifier,
                      loc=loc,
                      shape=shape,
                      scale=scale,
                      group_size=group_size,
                      replicates=replicates,
                      estimated_100_return=estimated_100_return,
                      estimated_1000_return=estimated_1000_return,
                      n_over_100_return=sum(rand>estimated_100_return),
                      n_over_1000_return=sum(rand>estimated_1000_return),
                      AD=get_a_d_stat(rand,
                                      CDF_function,
                                      loc=loc,
                                      scale=scale,
                                      shape=shape),
                      UMAD = get_m_a_d_u_stat(rand,
                                              CDF_function,
                                              loc=loc,
                                              scale=scale,
                                              shape=shape)
                      ),
    by=list(identifier=group)])
}

gev_prof_ret_mu<-function(return_level, p, scale, shape){
  if(abs(shape)<=0.0000001){
    return_level + scale*log(-log(1-p))
  }else{
    return_level + (scale/shape)*(1-(-log(1-p))^(-shape))
  }
  
}

gev_prof_lik_return<-function(par, Y, return_level, p){
  scale=par[['scale']]
  shape=par[['shape']]
  if(scale<=0){
    -.Machine$integer.max
  }
  loc=gev_prof_ret_mu(return_level, p, scale, shape)
  if (is.na(loc)|is.infinite(loc)){
    #Do not allow invalid parameters
    -.Machine$integer.max
  }
  m=length(Y)
  alpha = -m*log(scale)
  S=(Y-loc)/scale
  if(any((1+(shape*S))<=0)){
    #Do not allow invalid parameters
    -.Machine$integer.max
  }
  if(abs(shape)<=0.0000001){
    -(alpha - sum(S) - sum(exp(-S)))
  }else{
    inner_sum_1<-log(1+(shape*S))
    inner_sum_2<-(1+(shape*S))^(-1/shape)
    -(alpha - (1+(1/shape))*sum(inner_sum_1) - sum(inner_sum_2))
  }
}


get_gev_profile_return_log_likelihood<-function(return_levels, Y, p, start_scale, start_shape){
  results<-numeric(length(return_levels))
  par=c(scale=start_scale,
    shape=start_shape)
  for (i in 1:length(return_levels)){
    fit<-optim(par=par,
                   fn=gev_prof_lik_return,
                   method="BFGS",
                   Y=Y,
                   return_level=return_levels[i],
                   p=p)
    par=fit$par
    results[i]<-fit$value
  }
  results
}

get_earthquake_data<-function(folder_path){
  files<-list.files(folder_path,pattern="*.catalog",full.names = T)
  all_data<-vector(mode='list',length=length(files))
  for(i in 1:length(files)){
    data<-fread(files[[i]],skip=6)
    setnames(data,
             c("#YYY/MM/DD","HH:mm:SS.ss"),
             c("DATE","TIME"))
    names(data)<-tolower(names(data))
    all_data[[i]]<-data
  }
  all_data<-rbindlist(all_data)
  all_data[,timestamp:=as.POSIXct(paste(date,time),tz='PDT')]
  all_data
}

get_globaltemp_data<-function(folder_path){
  all_data<-fread(file.path(folder_path,'GlobalTemperatures.csv'))
  all_data[,timestamp:=as.POSIXct(dt)]
  all_data
}

get_profile_confint<-function(logliks,parameter,max_loglik,p=0.95,plot=FALSE){
  threshold=max_loglik-qchisq(p=p,df=1)/2
  if(max(logliks)<threshold){
    print('No log likelihoods are high enough to be plausible. Check model assumptions and parameter span.')
    return(c(lower=NA,upper=NA))
  }
  max_lik_par=parameter[logliks==max(logliks)]
  #Get smooth spline for predicting intervals. Probably a more efficient way to do this, but sod it.
  spline_1<-splinefun(x=logliks[parameter>=max_lik_par],
                      y=parameter[parameter>=max_lik_par])
  spline_2<-splinefun(x=logliks[parameter<=max_lik_par],
                      y=parameter[parameter<=max_lik_par])
  
  if(plot==T){
    print(ggplot(data=data.table(x=parameter,y=logliks),aes(x=x,y=y))+geom_point()+
      geom_hline(yintercept = threshold,col='red')+
      geom_vline(xintercept = spline_1(threshold),col='red')+
      geom_vline(xintercept = spline_2(threshold),col='red'))
  }
  c(lower=spline_2(threshold),upper=spline_1(threshold))
}


get_all_profile_confints<-function(fit,p=0.1){
  confint<-confint(
        profile(fit))
  return_val_check<-get_gev_return_value_and_std.err(p,fit)
  
  seq_of_interest<-seq(return_val_check['estimate']-2.8*return_val_check['se'],
                       return_val_check['estimate']+2.8*return_val_check['se'],
                       length.out=30)
  
  
  thing<--get_gev_profile_return_log_likelihood(seq_of_interest, fit$data, 
                                                p, 1, 0)
  return_confints<-get_profile_confint(thing,seq_of_interest,-fit$deviance/2,plot=F)
  names(return_confints)<-c('return_lower','return_upper')
  rbind(confint,t(return_confints))
  
}

get_profile_confints<-function(fits){
  
  fitlist<-vector('list',.NUM_CORES)
  for(i in 1:.NUM_CORES){
    block_length<-ceiling(length(fits)/.NUM_CORES)
    fitlist[[i]]<-fits[seq(block_length*(i-1)+1,block_length*i)]
  }
  profiles=mclapply(fitlist,function(x){
    lapply(x,function(y){
      get_all_profile_confints(y)}
    )},mc.cores = .NUM_CORES)
  profiles<-rbindlist(unlist(profiles,recursive = FALSE))
}


get_pot_fits<-function(data){
  #Requires data of type [identifier,group,sample] where sample is all data points
  data[,list(identifier=identifier[1],
             fit=list(fpot(sample,threshold=quantile(sample,0.95)))),by=group]
  
}

get_pot_fits_no_std<-function(data){
  #Requires data of type [identifier,group,sample] where sample is all data points
  data[,list(identifier=identifier[1],
             fit=list(fpot(sample,threshold=min(sample)-0.000000001,
                           std.err=F))),by=group]
  
}

get_gev_fits<-function(data){
  #Requires data of type [identifier,group,sample] where sample is all data points
  data[,list(identifier=identifier[1],
             fit=list(fgev(sample))),by=group]
  
}

get_both_fits<-function(data){
  #Expects data with block (for blocking), group (for fit), identifier and sample
  grouped_maxima<-data[,list(identifier=identifier[1],
                             sample=max(sample)),
                       by=list(group,block)]
  
  gev_fits<-grouped_maxima[,list(identifier=identifier[1],
                                 fit=list(fgev(sample))),by=group]
  pot_fits<-data[,list(identifier=identifier[1],
                       fit=list(fpot(sample,threshold=quantile(sample,0.95)))),by=group]
  return<-merge(gev_fits,pot_fits,by=c('identifier','group'),suffixes=c('gev','gpd'))
  return(return)
}

get_both_fits_shape_0<-function(data){
  #Expects data with block (for blocking), group (for fit), identifier and sample
  grouped_maxima<-data[,list(identifier=identifier[1],
                             sample=max(sample)),
                       by=list(group,block)]
  
  gev_fits<-grouped_maxima[,list(identifier=identifier[1],
                                 fit=list(fgev(sample,shape=0))),by=group]
  pot_fits<-data[,list(identifier=identifier[1],
                       fit=list(fpot(sample,threshold=quantile(sample,0.95),shape=0))),by=group]
  return<-merge(gev_fits,pot_fits,by=c('identifier','group'),suffixes=c('gev','gpd'))
  return(return)
}

fit_safety=function(fit_function,...){
  tryCatch({return(fit_function(...))},
           error = function(e) {
            return(NULL)
           }
  )

}


get_both_fits_and_profiles<-function(data){
  #Expects data with block (for blocking), group (for fit), identifier and sample
  grouped_maxima<-data[,list(identifier=identifier[1],
                             sample=max(sample)),
                       by=list(group,block)]
  
  gev_fits<-grouped_maxima[,list(identifier=identifier[1],
                         fit=list(fit_safety(fgev,sample)),
                         fit_10=list(fit_safety(fgev,sample,prob=0.1)),
                         fit_100=list(fit_safety(fgev,sample,prob=0.01)),
                         fit_1000=list(fit_safety(fgev,sample,prob=0.001))),by=group]
  gev_fits[,`:=`(nullfit=is.null(fit[[1]]),
                 nullfit10=is.null(fit_10[[1]]),
                 nullfit100=is.null(fit_100[[1]]),
                 nullfit1000=is.null(fit_1000[[1]])),by=group]
  
  pot_fits<-data[,list(identifier=identifier[1],
                       fit=list(fit_safety(fpot,x=sample,threshold=quantile(sample,0.95))),
                       fit_10=list(fit_safety(fpot,x=sample,threshold=quantile(sample,0.95),
                                              mper=10,npp=100)),
                       fit_100=list(fit_safety(fpot,x=sample,threshold=quantile(sample,0.95),
                                               mper=100,npp=100)),
                       fit_1000=list(fit_safety(fpot,x=sample,threshold=quantile(sample,0.95),
                                                mper=1000,npp=100))),by=group]

  
  pot_fits[,`:=`(nullfit=is.null(fit[[1]]),
                 nullfit10=is.null(fit_10[[1]]),
                 nullfit100=is.null(fit_100[[1]]),
                 nullfit1000=is.null(fit_1000[[1]])),by=group]
  return<-merge(gev_fits,pot_fits,by=c('identifier','group'),suffixes=c('gev','gpd'))
  return[nullfitgev==F,gev_profile:=list(list(profile(fitgev[[1]]))),by=group]
  return[nullfitgpd==F,gpd_profile:=list(list(profile(fitgpd[[1]]))),by=group]
  return[nullfit10gev==F,gev_profile_10:=list(list(profile(fit_10gev[[1]],which='quantile'))),by=group]
  return[nullfit10gpd==F,gpd_profile_10:=list(list(profile(fit_10gpd[[1]],which='rlevel'))),by=group]
  return[nullfit100gev==F,gev_profile_100:=list(list(profile(fit_100gev[[1]],which='quantile'))),by=group]
  return[nullfit100gpd==F,gpd_profile_100:=list(list(profile(fit_100gpd[[1]],which='rlevel'))),by=group]
  return[nullfit1000gev==F,gev_profile_1000:=list(list(profile(fit_1000gev[[1]],which='quantile'))),by=group]
  return[nullfit1000gpd==F,gpd_profile_1000:=list(list(profile(fit_1000gpd[[1]],which='rlevel'))),by=group]
  return(return)
}

get_parallel_fits<-function(data,fits_function,num_cores=1){
  # Expects blocked or subset data ready to fit of format group, identifier, sample
  print(paste('running with ',num_cores,' cores.'))
  data[,par_group:=ceiling((group/max(group))*num_cores)]
  setkey(data,'par_group')
  return(foreach(x=1:num_cores, .combine="rbind", .inorder=FALSE) %dopar% 
           fits_function(data[.(x)]))
}

assign_stats_from_fits<-function(data,include_profiles=T){
  data[,`:=`(
   gpd_confint=list(confint(fitgpd[[1]])),
   gev_confint=list(confint(fitgev[[1]]))

   ),by=group]
  
  
  if(include_profiles==T){
    data[,`:=`(
      gev_profile_confints=list(confint(gev_profile[[1]])),
      gpd_profile_confints=list(confint(gpd_profile[[1]])),
      gev_10_return=list(get_gev_return_value_and_std.err(1/10,
                                                          fitgev[[1]])),
      gev_100_return=list(get_gev_return_value_and_std.err(1/100,
                                                           fitgev[[1]])),
      gev_1000_return=list(get_gev_return_value_and_std.err(1/1000,
                                                            fitgev[[1]])),
      gpd_10_return=list(get_gpd_return_value_and_std.err(1/10,
                                                          fitgpd[[1]])),
      gpd_100_return=list(get_gpd_return_value_and_std.err(1/100,
                                                           fitgpd[[1]])),
      gpd_1000_return=list(get_gpd_return_value_and_std.err(1/1000,
                                                            fitgpd[[1]])),
      gpd_10000_return=list(get_gpd_return_value_and_std.err(1/10000,
                                                             fitgpd[[1]])),
      gpd_100000_return=list(get_gpd_return_value_and_std.err(1/100000,
                                                              fitgpd[[1]]))
    ),by=group]
    data[nullfit10gev==F,`:=`(gev_10_confint=list(confint(fit_10gev[[1]])),
                              gev_profile_10_confints=list(confint(gev_profile_10[[1]]))),by=group]
    data[nullfit100gev==F,`:=`(gev_100_confint=list(confint(fit_100gev[[1]])),
                               gev_profile_100_confints=list(confint(gev_profile_100[[1]]))),by=group]
    data[nullfit10gpd==F,`:=`(gpd_10_confint=list(confint(fit_10gpd[[1]])),
                              gpd_profile_10_confints=list(confint(gpd_profile_10[[1]]))),by=group]
    data[nullfit100gpd==F,`:=`(gpd_100_confint=list(confint(fit_100gpd[[1]])),
                               gpd_profile_100_confints=list(confint(gpd_profile_100[[1]]))),by=group]
  }
  
  data
}

get_runforward_data<-function(data,splits=11){
  #Expects data in format identifier, sample, block. Group will be created.
  #Work out thresholds for data selection: Exxlude thresholds at ends. 
  #Even groups are post-break, odd pre-break.
  breaks<-ceiling(seq(data[,min(block)],data[,max(block)],length.out = splits))[2:(splits-1)]
  results<-vector(mode='list',length(breaks))
  for(i in 1:length(breaks)){
    results[[i]]<-data[block<=breaks[[i]],list(breakpoint=breaks[[i]],
                             identifier,
                             sample,
                             block,
                             group=i)]
  }
  results<-rbindlist(results)
  results
}

get_estimate_agreements<-function(data){
  data[,list(gev_shape_in_gpd_normal=inrange(fitgev[[1]]$estimate['shape'],
                                                   gpd_confint[[1]]['shape','2.5 %'],
                                                   gpd_confint[[1]]['shape','97.5 %']),
                   gpd_shape_in_gev_normal=inrange(fitgpd[[1]]$estimate['shape'],
                                                   gev_confint[[1]]['shape','2.5 %'],
                                                   gev_confint[[1]]['shape','97.5 %']),
                   gev_100_in_gpd_normal=inrange(gev_100_return[[1]][['estimate']],
                                                 gpd_100_confint[[1]]['rlevel','2.5 %'],
                                                 gpd_100_confint[[1]]['rlevel','97.5 %']),
                   gpd_100_in_gev_normal=inrange(gpd_10000_return[[1]][['estimate']],
                                                 gev_100_confint[[1]]['quantile','2.5 %'],
                                                 gev_100_confint[[1]]['quantile','97.5 %']),
                   
                   gev_shape_in_gpd_profile=inrange(fitgev[[1]]$estimate['shape'],
                                                    gpd_profile_confints[[1]]['shape','lower'],
                                                    gpd_profile_confints[[1]]['shape','upper']),
                   gpd_shape_in_gev_profile=inrange(fitgpd[[1]]$estimate['shape'],
                                                    gev_profile_confints[[1]]['shape','lower'],
                                                    gev_profile_confints[[1]]['shape','upper']),
                   gev_100_in_gpd_profile=inrange(gev_100_return[[1]][['estimate']],
                                                  gpd_profile_100_confints[[1]]['rlevel','lower'],
                                                  gpd_profile_100_confints[[1]]['rlevel','upper']),
                   gpd_100_in_gev_profile=inrange(gpd_10000_return[[1]][['estimate']],
                                                  gev_profile_100_confints[[1]]['quantile','lower'],
                                                  gev_profile_100_confints[[1]]['quantile','upper'])
  ),by=group][,list(gev_shape_in_gpd_normal=mean(gev_shape_in_gpd_normal),
                    gpd_shape_in_gev_normal=mean(gpd_shape_in_gev_normal),
                    gev_100_in_gpd_normal=mean(gev_100_in_gpd_normal,na.rm=T),
                    gpd_100_in_gev_normal=mean(gpd_100_in_gev_normal,na.rm=T),
                    gev_shape_in_gpd_profile=mean(gev_shape_in_gpd_profile),
                    gpd_shape_in_gev_profile=mean(gpd_shape_in_gev_profile),
                    gev_100_in_gpd_profile=mean(gev_100_in_gpd_profile,na.rm=T),
                    gpd_100_in_gev_profile=mean(gpd_100_in_gev_profile,na.rm=T))]
}

get_differences_table<-function(data,ident,upper=1,lower=1){
  differences=data[,list(identifier=ident,I=.I,
                         shape = (upper*gev_confint[[1]]['shape','97.5 %'])-(lower*gev_confint[[1]]['shape','2.5 %']),
                         shape_profile=(upper*gev_profile_confints[[1]]['shape','upper'])-(lower*gev_profile_confints[[1]]['shape','lower']),
                         loc = (upper*gev_confint[[1]]['loc','97.5 %'])-(lower*gev_confint[[1]]['loc','2.5 %']),
                         loc_profile=(upper*gev_profile_confints[[1]]['loc','upper'])-(lower*gev_profile_confints[[1]]['loc','lower']),
                         scale = (upper*gev_confint[[1]]['scale','97.5 %'])-(lower*gev_confint[[1]]['scale','2.5 %']),
                         scale_profile=(upper*gev_profile_confints[[1]]['scale','upper'])-(lower*gev_profile_confints[[1]]['scale','lower']),
                         return_10 = (upper*gev_10_confint[[1]]['quantile','97.5 %'])-(lower*gev_10_confint[[1]]['quantile','2.5 %']),
                         return_10_profile=(upper*gev_profile_10_confints[[1]]['quantile','upper'])-(lower*gev_profile_10_confints[[1]]['quantile','lower']),
                         return_100 = (upper*gev_100_confint[[1]]['quantile','97.5 %'])-(lower*gev_100_confint[[1]]['quantile','2.5 %']),
                         return_100_profile=(upper*gev_profile_100_confints[[1]]['quantile','upper'])-(lower*gev_profile_100_confints[[1]]['quantile','lower']),
                         gpd_shape = (upper*gpd_confint[[1]]['shape','97.5 %'])-(lower*gpd_confint[[1]]['shape','2.5 %']),
                         gpd_shape_profile=(upper*gpd_profile_confints[[1]]['shape','upper'])-(lower*gpd_profile_confints[[1]]['shape','lower']),
                         gpd_scale = (upper*gpd_confint[[1]]['scale','97.5 %'])-(lower*gpd_confint[[1]]['scale','2.5 %']),
                         gpd_scale_profile=(upper*gpd_profile_confints[[1]]['scale','upper'])-(lower*gpd_profile_confints[[1]]['scale','lower']),
                         gpd_return_10 = (upper*gpd_10_confint[[1]]['rlevel','97.5 %'])-(lower*gpd_10_confint[[1]]['rlevel','2.5 %']),
                         gpd_return_10_profile=(upper*gpd_profile_10_confints[[1]]['rlevel','upper'])-(lower*gpd_profile_10_confints[[1]]['rlevel','lower']),
                         gpd_return_100 = (upper*gpd_100_confint[[1]]['rlevel','97.5 %'])-(lower*gpd_100_confint[[1]]['rlevel','2.5 %']),
                         gpd_return_100_profile=(upper*gpd_profile_100_confints[[1]]['rlevel','upper'])-(lower*gpd_profile_100_confints[[1]]['rlevel','lower'])
  ),by=group]
  differences
}

offset_lines_and_data<-function(differences,column_name,display_name,offset,prefix){
  ggplot(data=differences[identifier=='w=0.5',
                                   list(x=get(column_name),
                                        y=get(paste0(column_name,'_profile')),
                                        Distribution=identifier)],
         aes(x=x,y=y,col=Distribution))+
       geom_point(size=0.01)+
       geom_point(data=differences[identifier=='w=1',
                                   list(x=get(column_name),
                                        y=get(paste0(column_name,'_profile'))+offset,
                                        Distribution=identifier)],
                  size=0.01)+
      geom_point(data=differences[identifier=='w=0',
                                  list(x=get(column_name),
                                       y=get(paste0(column_name,'_profile'))-offset,
                                       Distribution=identifier)],
                 size=0.01)+
       geom_abline(slope=1,intercept=0,linetype='dotted')+
       geom_abline(slope=1,intercept=offset,linetype='dotted')+
       geom_abline(slope=1,intercept=-offset,linetype='dotted')+
    theme_classic()+theme(legend.position='none')+
  xlab(paste0(prefix," normal approximation 95% CI of ", display_name))+
  scale_y_continuous(name = paste0(prefix,' profile Likelihood 95% CI\nof ', display_name),breaks=NULL)
}

confidence_range_comparison_plots<-function(rand_stats, broken_stats,
                                            very_broken_stats, upper=1, lower=1,
                                            prefix="Width of"){
  differences<-rbind(get_differences_table(rand_stats,'w=0',upper,lower),
  get_differences_table(broken_stats,'w=0.5',upper,lower),
  get_differences_table(very_broken_stats,'w=1',upper,lower))
  
  return(list(
    gev_shape=offset_lines_and_data(differences,'shape', 'GEV shape parameter', 0.05, prefix),
    gev_loc=offset_lines_and_data(differences,'loc', 'GEV location parameter', 0.05, prefix),
    gev_scale=offset_lines_and_data(differences,'scale', 'GEV scale parameter', 0.05, prefix),
    gev_return_10=offset_lines_and_data(differences,'return_10', 'GEV 10 block return level', 0.2, prefix),
    gev_return_100=offset_lines_and_data(differences,'return_100', 'GEV 100 block return level', 2, prefix),
    gpd_shape=offset_lines_and_data(differences,'gpd_shape', 'GPD shape parameter', 0.05, prefix),
    gpd_scale=offset_lines_and_data(differences,'gpd_scale', 'GPD 10 scale parameter', 0.05, prefix),
    gpd_return_10=offset_lines_and_data(differences,'return_10', 'GPD 10 block return level', 0.5, prefix),
    gpd_return_100=offset_lines_and_data(differences,'return_100', 'GPD 100 block return level', 2, prefix)
  ))
}


