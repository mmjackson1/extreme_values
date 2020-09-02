get_both_ad_of_interest<-function(Y,CDF=pgev, ...){
  data.table(AD=get_a_d_stat(Y,CDF,...),
             UMAD = get_m_a_d_u_stat(Y,CDF,...))
}

get_triple_threat_ads<-function(m,
                                loc=0,
                                scale = 1,
                                shape=0,
                                CDF=pgev){
  sample=rgev(m,
              loc=loc,
              scale = scale,
              shape = shape)
  fit<-fgev(sample,std.err = FALSE)
  return(get_both_ad_of_interest(sample,
                                 CDF=CDF,
                                 loc=fit$estimate[['loc']],
                                 scale=fit$estimate[['scale']],
                                 shape=fit$estimate[['shape']]))
}

bootstrap_ad_umad_test<-function(m=50,
                                 loc = 0,
                                 scale = 1,
                                 shape = 0,
                                 bootstrap_size=100000){
  results<-mclapply(1:bootstrap_size,function(x){get_triple_threat_ads(m,
                                                                       loc=loc,
                                                                       scale = scale,
                                                                       shape=shape)},
                    mc.cores = 32)
  results<-rbindlist(results)
  
  AD_crits<-quantile(results$AD,c(0.9,0.95,0.99,0.995))
  UMAD_crits<-quantile(results$UMAD,c(0.9,0.95,0.99,0.995))
  return(data.table(size=m,
                    loc=loc,
                    scale=scale,
                    shape=shape,
                    bootstrap_size=bootstrap_size,
                    AD.9=AD_crits[[1]],
                    AD.95=AD_crits[[2]],
                    AD.99=AD_crits[[3]],
                    AD.995=AD_crits[[4]],
                    AD.SD=sd(results$AD),
                    UMAD.9=UMAD_crits[[1]],
                    UMAD.95=UMAD_crits[[2]],
                    UMAD.99=UMAD_crits[[3]],
                    UMAD.995=UMAD_crits[[4]],
                    UMAD.SD=sd(results$UMAD)))
  
}

#Designed to be used as a small part of multiple overnight runs.
sizes = c(200,300,400,750)
shapes = c(-0.4,-0.3,-0.05,0.05,0.3,0.4)
scales = c(0.3,0.75,1.5,1,1.25)
locs = c(-2,-1.5,0,1.5,2)

n=length(sizes)*length(shapes)*length(scales)*length(locs)

pars<-data.table(size=sizes,ID=1)
pars<-merge(pars,
            data.table(shape=shapes,ID=1),
            allow.cartesian=T)
pars<-merge(pars,
            data.table(scale=scales,ID=1),
            allow.cartesian=T)
pars<-merge(pars,
            data.table(loc=locs,ID=1),
            allow.cartesian=T)
pars[,ID:=NULL]

#Use if part completed run
if(exists(results_bound)){
  pars=setDT(pars)[!results_bound, on = c('size','scale','loc','shape')]
}


results<-vector('list',length=pars[,.N])

for(i in 1:pars[,.N]){
  print(i/pars[,.N])
  test_stats<-bootstrap_ad_umad_test(m=pars[i,size],
                                     shape=pars[i,shape],
                                     scale=pars[i,scale],
                                     loc=pars[i,loc],
                                     bootstrap_size=100000)
  results[[i]]<-test_stats
}

results_bound<-rbindlist(results[!is.na(results)])

results_bk<-copy(results)

###Begin from here if analysing already-saved results

results_bound_1<-fread('~/bootstrapped_crits_1.csv')
results_bound_2<-fread('~/bootstrapped_crits_2.csv')
results_bound_3<-fread('~/bootstrapped_crits_3.csv')
results_bound_4<-fread('~/bootstrapped_crits_4.csv')


results_composite<-rbind(results_bound_1,
                         results_bound_2,
                         results_bound_3,
                         results_bound_4,
                         fill=T)

results_composite<-unique(results_composite)
results_composite

ggplot(results_composite[shape==0&size>25&bootstrap_size==1e+05],
       aes(x=scale,y=AD.95,
           col=as.factor(size)))+
  geom_point(alpha=0.5)+theme_classic()+
  scale_colour_discrete(name = "Sample size")+
  xlab('Scale parameter')+ylab('95th quantile of AD statistic')
ggsave('~/ADstats_scale.png')
ggplot(results_composite[shape==0&size>25&bootstrap_size==1e+05],
       aes(x=scale,y=UMAD.95,
           col=as.factor(size)))+
  scale_colour_discrete(name = "Sample size")+
  geom_point(alpha=0.5)+theme_classic()+
  xlab('Scale parameter')+ylab('95th quantile of UMAD statistic')
ggsave('~/UMADstats_scale.png')

ggplot(results_composite[shape==0&size>25&bootstrap_size==1e+05],
       aes(x=loc,y=AD.95,
           col=as.factor(size)))+
  geom_point(alpha=0.5)+theme_classic()+
  scale_colour_discrete(name = "Sample size")+
  xlab('Location parameter')+ylab('95th quantile of AD statistic')
ggsave('~/ADstats_loc.png')
ggplot(results_composite[shape==0&size>25&bootstrap_size==1e+05],
       aes(x=loc,y=UMAD.95,
           col=as.factor(size)))+
  geom_point(alpha=0.5)+theme_classic()+
  scale_colour_discrete(name = "Sample size")+
  xlab('Location parameter')+ylab('95th quantile of AD statistic')
ggsave('~/UMADstats_loc.png')

ggplot(results_composite[size>25&shape>-0.3],
       aes(x=shape,y=AD.95,
           col=as.factor(size)))+
  scale_colour_discrete(name = "Sample size")+
  geom_point(alpha=0.5)+theme_classic()+geom_smooth()+
  xlab('Shape parameter')+ylab('95th quantile of AD statistic')
ggsave('~/ADstats_shape.png')
ggplot(results_composite[size>25&shape>-0.3],
       aes(x=shape,y=UMAD.95,
           col=as.factor(size)))+
  scale_colour_discrete(name = "Sample size")+
  geom_point(alpha=0.5)+theme_classic()+geom_smooth()+
  xlab('Shape parameter')+ylab('95th quantile of AD statistic')
ggsave('~/UMADstats_shape.png')

#Fixed size test to validate results vs dagostino and stephens


get_fixed_shape<-function(m,
                          loc=0,
                          scale = 1,
                          shape=0,
                          estshape = 0){
  sample=rgev(m,
              loc=loc,
              scale = scale,
              shape = shape)
  return(get_both_ad_of_interest(sample,
                                 CDF=pgev,
                                 loc=loc,
                                 scale=scale,
                                 shape=shape))
}



#Check to make sure we're matching the appropriate tables

bootstrap_fixed_size_test<-function(m=50,
                                    loc = 0,
                                    scale = 1,
                                    shape = 0,
                                    bootstrap_size=100000){
  results<-mclapply(1:bootstrap_size,function(x){get_fixed_shape(m,
                                                                 loc=loc,
                                                                 scale = scale,
                                                                 shape=shape)},
                    mc.cores = 32)
  results<-rbindlist(results)
  
  AD_crits<-quantile(results$AD,c(0.9,0.95,0.99))
  UMAD_crits<-quantile(results$UMAD,c(0.9,0.95,0.99))
  return(list(AD=AD_crits,UMAD=UMAD_crits))
  
}

sizes_2=c(20,30,40,50,60,80,90,100,125,150,175,200,300,400,500,600,700,800,900,1000)
checkit_size_2<-vector(mode="list",length=length(sizes_2))
for(i in 1:length(sizes_2)){
  print(sizes_2[[i]])
  checkit_size_2[[i]]<-bootstrap_fixed_size_test(m=sizes_2[[i]],
                                                 bootstrap_size = 10000)
}
checkit_size_2
ADs_size_2<-vector(mode="list",length=length(sizes_2))
MADs_size_2<-vector(mode="list",length=length(sizes_2))
for(i in 1:length(sizes_2)){
  ADs_size_2[[i]]<-as.data.table(checkit_size_2[[i]]$AD[[2]])
  MADs_size_2[[i]]<-as.data.table(checkit_size_2[[i]]$UMAD[[2]])
}
ADs_size_2<-rbindlist(ADs_size_2)
MADs_size_2<-rbindlist(MADs_size_2)

plot(sizes_2[sizes_2>=20],ADs_size_2$V1[sizes_2>=20])
plot(sizes_2[sizes_2>=20],MADs_size_2$V1[sizes_2>=20])