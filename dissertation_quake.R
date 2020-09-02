quake<-get_earthquake_data('~/dissertation/SCEC_DC/')
quake<-quake[year(date)>1932&mag>3.25]
quake[,I:=.I]
quake[,yearmonth:=as.Date(paste0(year(date),'-',month(date),'-','01'))]
thresh=quake[order(-mag)][,quantile(mag,0.95)]
bt=quake[mag<=thresh,
         list(mag,I,Date=as.Date(date),
              POT = 'Below Threshold',
              alpha=0.05)]
abt=quake[mag>thresh,
          list(mag,I,Date=as.Date(date),
               POT = 'Above Threshold',
               alpha=1)]
toplot=rbind(abt,bt)

ggplot(toplot,aes(x=I,y=mag))+
  geom_point(alpha=0.07,data=bt,aes(col=POT))+
  geom_point(alpha=0.5,data=abt,aes(col=POT))+
  geom_smooth(method='loess',se=F,span=0.1)+theme_classic()+
  xlab("Index")+
  ylab("Magnitude")+
  geom_hline(yintercept=thresh,col='red')+
  scale_color_manual(values=c("red", "black"))

ggsave("~/earthquake_pot.png")

ggplot(quake[year(date)>1932&mag>3.25,
             list(mag=max(mag)),by=list(group=round(I/100))],aes(x=group,y=mag))+
  geom_point()+theme_classic()+
  xlab("Group")+
  ylab("Magnitude")

ggsave("~/earthquake_groups.png")

quake_data<-quake[year(date)>1932&mag>3.25,list(identifier='quake_unmodified-100',
                                                sample=mag,
                                                block = round(I/100),
                                                group=1)]
quake_models<-get_both_fits_and_profiles(quake_data)

quake_models<-assign_stats_from_fits(quake_models)

quake_models[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
                  by=group]
quake_models[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)-3.345345e-05),
                  by=group]
quake_models[,inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err)]
quake_models[,pnorm(shape_diff,0,shape_diff_err)*100]

par(mfrow=c(2,2))
plot(quake_models[,fitgev[[1]]])
plot(quake_models[,fitgpd[[1]]])
par(mfrow=c(1,1))

quake_models[,list(AD_gev=get_a_d_stat(fitgev[[1]]$data,CDF=pgev,
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
set.seed(607)
quake_sim_gev_data<-data.table(sample=rgev(100000*129,
                                            loc=quake_models[,fitgev[[1]]$estimate[['loc']]],
                                            scale=quake_models[,fitgev[[1]]$estimate[['scale']]],
                                            shape=quake_models[,fitgev[[1]]$estimate[['shape']]]),
                 group=rep(1:100000,129),
                 identifier='quake_simulated-129')
quake_sim_gev_fits<-get_parallel_fits(quake_sim_gev_data,get_gev_fits,num_cores=8)
quake_sim_gev_fits[,`:=`(ad_stat=get_a_d_stat(Y=fit[[1]]$data,
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

quake_sim_gpd_data<-data.table(sample=rgpd(100000*641,
                                           loc=quake_models[,fitgpd[[1]]$threshold],
                                           scale=quake_models[,fitgpd[[1]]$estimate[['scale']]],
                                           shape=quake_models[,fitgpd[[1]]$estimate[['shape']]]),
                               group=rep(1:100000,641),
                               identifier='quake_simulated-gpd-641')
quake_sim_gpd_fits<-get_parallel_fits(quake_sim_gpd_data,get_pot_fits_no_std,num_cores=8)
quake_sim_gpd_fits[,`:=`(ad_stat=get_a_d_stat(Y=fit[[1]]$data,
                                              CDF=pgpd,
                                              loc=fit[[1]]$threshold,
                                              scale=fit[[1]]$estimate[['scale']],
                                              shape=fit[[1]]$estimate[['shape']]),
                         umad_stat=get_m_a_d_u_stat(Y=fit[[1]]$data,
                                                    CDF=pgpd,
                                                    loc=fit[[1]]$threshold,
                                                    scale=fit[[1]]$estimate[['scale']],
                                                    shape=fit[[1]]$estimate[['shape']])),
                   by=group]

quake_sim_gev_fits[,list('0.1'=quantile(ad_stat,0.9),
                         '0.05'=quantile(ad_stat,0.95),
                         '0.01'=quantile(ad_stat,0.99),
                         '0.005'=quantile(ad_stat,0.995))]
quake_sim_gev_fits[,list('0.1'=quantile(umad_stat,0.9),
                         '0.05'=quantile(umad_stat,0.95),
                         '0.01'=quantile(umad_stat,0.99),
                         '0.005'=quantile(umad_stat,0.995))]

quake_sim_gpd_fits[,list('0.1'=quantile(ad_stat,0.9),
                         '0.05'=quantile(ad_stat,0.95),
                         '0.01'=quantile(ad_stat,0.99),
                         '0.005'=quantile(ad_stat,0.995))]
quake_sim_gpd_fits[,list('0.1'=quantile(umad_stat,0.9),
                         '0.05'=quantile(umad_stat,0.95),
                         '0.01'=quantile(umad_stat,0.99),
                         '0.005'=quantile(umad_stat,0.995))]



####runforward####

quake_runforward_data=get_runforward_data(quake_data,splits=64)
quake_runforward_data=quake_runforward_data[group>=9]

quake_runforward_models<-get_both_fits(quake_runforward_data[,list(block,group,
                                                                 identifier=breakpoint,
                                                                 sample=sample)])

quake_runforward_models<-assign_stats_from_fits(quake_runforward_models,include_profiles = F)
quake_runforward_models[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
             by=group]
quake_runforward_models[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)-3.345345e-05),
             by=group]
quake_runforward_models[,mean(inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err))]
quake_runforward_models[,pnorm(shape_diff,0,shape_diff_err)*100]
ggplot(quake_runforward_models[,list(identifier,Quantile=pnorm(shape_diff,0,shape_diff_err)*100)],
       aes(x=identifier,y=Quantile))+geom_point()+theme_classic()
ggsave('~/quake_runforward.png')


