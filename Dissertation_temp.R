#Data exploration and plot generation

#Data exploration



#TEMPERATURES
temp<-get_globaltemp_data('~/dissertation/')[!is.na(LandMaxTemperature)&dt>"1900-01-01"]
temp_groups<-temp[!is.na(LandMaxTemperature)&dt>"1900-01-01",
                  list(temperature=max(LandMaxTemperature),
                       .N),
                  by=list(year=year(dt))]


thresh=temp[order(-LandMaxTemperature)][,quantile(LandMaxTemperature,0.95)]

bt=temp[LandMaxTemperature<=thresh,
         list(LandMaxTemperature,Date=as.Date(dt),
              POT = 'Below Threshold')]
abt=temp[LandMaxTemperature>thresh,
         list(LandMaxTemperature,Date=as.Date(dt),
              POT = 'Above Threshold')]
toplot=rbind(abt,bt)
ggplot(toplot,aes(x=Date,y=LandMaxTemperature))+
  geom_point(alpha=0.5,aes(col=POT))+
  geom_line(alpha=0.5)+
  geom_smooth(se=F,method='loess')+
  theme_classic()+
  xlab("Date")+
  ylab("Temperature (Degrees Celsius)")+
  geom_hline(yintercept=thresh,col='red')+
  scale_color_manual(values=c("red", "black"))

ggsave("~/temp_pot.png")

toplot<-rbind(
  temp[,
       list(LandMaxTemperature = max(LandMaxTemperature),
            LandMaxTop=max(LandMaxTemperature+LandMaxTemperatureUncertainty),
            LandMaxBottom=max(LandMaxTemperature-LandMaxTemperatureUncertainty)),
       by=list(Date=as.Date(paste0(year(as.Date(dt)),'-06-01')))]
  )

ggplot(toplot,aes(x=Date,y=LandMaxTemperature))+
  geom_point()+
  geom_errorbar(aes(ymin=LandMaxBottom, ymax=LandMaxTop), width=.2,
                position=position_dodge(.9),alpha=0.5) +
  geom_smooth(se=F,method='loess',alpha=0.5)+
  theme_classic()+
  xlab("Date")+
  ylab("Temperature (Degrees Celsius)")
ggsave("~/temp_groups.png")

temp_data<-temp[,list(identifier='temp_unmodified-100',
                        sample=LandMaxTemperature,
                        block =year(dt),
                        group=1)]
temp_models<-get_both_fits_and_profiles(temp_data)

temp_models<-assign_stats_from_fits(temp_models)
temp_models[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
             by=group]
temp_models[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)-3.345345e-05),
             by=group]
temp_models[,inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err)]
temp_models[,pnorm(shape_diff,0,shape_diff_err)*100]
plot(temp_models[,fitgev[[1]]])
plot(temp_models[,fitgpd[[1]]])

temp_models[,fitgev[[1]]]

temp_models[,list(AD_gev=get_a_d_stat(fitgev[[1]]$data,CDF=pgev,
                                       loc=fitgev[[1]]$estimate[['loc']],
                                       scale=fitgev[[1]]$estimate[['scale']],
                                       shape=fitgev[[1]]$estimate[['shape']]),
                   UMAD_gev=get_m_a_d_u_stat(fitgev[[1]]$data,CDF=pgev,
                                             loc=fitgev[[1]]$estimate[['loc']],
                                             scale=fitgev[[1]]$estimate[['scale']],
                                             shape=fitgev[[1]]$estimate[['shape']]),
                   AD_gpd=get_a_d_stat(fitgpd[[1]]$exceedances,CDF=pgpd,
                                       loc=fitgpd[[1]]$threshold,
                                       scale=fitgpd[[1]]$estimate[['scale']],
                                       shape=fitgpd[[1]]$estimate[['shape']]),
                   UMAD_gpd=get_m_a_d_u_stat(fitgpd[[1]]$exceedances,CDF=pgpd,
                                             loc=fitgpd[[1]]$threshold,
                                             scale=fitgpd[[1]]$estimate[['scale']],
                                             shape=fitgpd[[1]]$estimate[['shape']]))
]
#### AD_STATS ####
set.seed(336)
temp_sim_gev_data<-data.table(sample=rgev(100000*116,
                                           loc=temp_models[,fitgev[[1]]$estimate[['loc']]],
                                           scale=temp_models[,fitgev[[1]]$estimate[['scale']]],
                                           shape=temp_models[,fitgev[[1]]$estimate[['shape']]]),
                               group=rep(1:100000,116),
                               identifier='temp_simulated_gev-116')
temp_sim_gev_fits<-get_parallel_fits(temp_sim_gev_data,get_gev_fits,num_cores=8)
temp_sim_gev_fits[,`:=`(ad_stat=get_a_d_stat(Y=fit[[1]]$data,
                                             CDF=pgev,
                                             loc=fit[[1]]$estimate[['loc']],
                                             scale=fit[[1]]$estimate[['scale']],
                                             shape=fit[[1]]$estimate[['shape']]),
                        umad_stat=get_m_a_d_u_stat(Y=fit[[1]]$data,
                                                   CDF=pgev,
                                                   loc=fit[[1]]$estimate[['loc']],
                                                   scale=fit[[1]]$estimate[['scale']],
                                                   shape=fit[[1]]$estimate[['shape']])),
                  by=group]


temp_sim_gpd_data<-data.table(sample=rgpd(100000*69,
                                           loc=temp_models[,fitgpd[[1]]$threshold],
                                           scale=temp_models[,fitgpd[[1]]$estimate[['scale']]],
                                           shape=temp_models[,fitgpd[[1]]$estimate[['shape']]]),
                               group=rep(1:100000,69),
                               identifier='temp_simulated-gpd-69')
temp_sim_gpd_fits<-get_parallel_fits(temp_sim_gpd_data,get_pot_fits_no_std,num_cores=8)
temp_sim_gpd_fits[,`:=`(ad_stat=get_a_d_stat(Y=fit[[1]]$exceedances,
                                             CDF=pgpd,
                                             loc=fit[[1]]$threshold,
                                             scale=fit[[1]]$estimate[['scale']],
                                             shape=fit[[1]]$estimate[['shape']]),
                        umad_stat=get_m_a_d_u_stat(Y=fit[[1]]$exceedances,
                                                   CDF=pgpd,
                                                   loc=fit[[1]]$threshold,
                                                   scale=fit[[1]]$estimate[['scale']],
                                                   shape=fit[[1]]$estimate[['shape']])),
                  by=group]

temp_sim_gev_fits[,list('0.1'=quantile(ad_stat,0.9),
                         '0.05'=quantile(ad_stat,0.95),
                         '0.01'=quantile(ad_stat,0.99),
                         '0.005'=quantile(ad_stat,0.995))]
temp_sim_gev_fits[,list('0.1'=quantile(umad_stat,0.9),
                         '0.05'=quantile(umad_stat,0.95),
                         '0.01'=quantile(umad_stat,0.99),
                         '0.005'=quantile(umad_stat,0.995))]

temp_sim_gpd_fits[,list('0.1'=quantile(ad_stat,0.9),
                         '0.05'=quantile(ad_stat,0.95),
                         '0.01'=quantile(ad_stat,0.99),
                         '0.005'=quantile(ad_stat,0.995))]
temp_sim_gpd_fits[,list('0.1'=quantile(umad_stat,0.9),
                         '0.05'=quantile(umad_stat,0.95),
                         '0.01'=quantile(umad_stat,0.99),
                         '0.005'=quantile(umad_stat,0.995))]

####Runforward####

temp_runforward_data=get_runforward_data(temp_data,splits = 99)
temp_runforward_data<-temp_runforward_data[group>5]

temp_runforward_models<-get_both_fits(temp_runforward_data[,list(block,group,
                                                                     identifier='temp_runforward',
                                                                     sample=sample)])
temp_runforward_models<-merge(temp_runforward_models,
                              temp_runforward_data[,list(Year=max(block)),by=group],
                              by='group')
temp_runforward_models<-assign_stats_from_fits(temp_runforward_models,include_profiles = F)
temp_runforward_models[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
            by=group]
temp_runforward_models[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)-3.345345e-05),
            by=group]
temp_runforward_models[,inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err)]
ggplot(temp_runforward_models[,list(Year,Quantile=pnorm(shape_diff,0,shape_diff_err)*100)],
       aes(x=Year,y=Quantile))+geom_point()+theme_classic()
ggsave('~/temp_runforward.png')

####Bootstraps####
set.seed(10)
temp_bootstrap_data<-lapply(1:10000,function(x){temp_data[,list(sample=sample(temp_data$sample,replace = T),
                                                                  block,
                                                                  identifier='temp_bootstrap-100',
                                                                  group=x)]})
temp_bootstrap_data<-rbindlist(temp_bootstrap_data)

temp_bootstrap_fits<-get_parallel_fits(temp_bootstrap_data,get_both_fits_no_err,num_cores=8)
temp_bootstrap_cor = temp_bootstrap_fits[,list(a=fitgpd[[1]]$estimate[['shape']],
                                                 b=fitgev[[1]]$estimate[['shape']]),
                                           by=group][,cor(b,a)]

temp_models[,shape_diff_err_known_corr:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)-
                                                2*temp_bootstrap_cor*fitgpd[[1]]$std.err[['shape']]*fitgev[[1]]$std.err[['shape']]),
             by=group]
temp_models[,inrange(0,shape_diff-1.96*shape_diff_err_known_corr,shape_diff+1.96*shape_diff_err_known_corr)]
temp_models[,pnorm(shape_diff,0,shape_diff_err_known_corr)*100]
