#Checking correlations of shapes across multiple parameter sets
load(file='~/dissertation/rand_stats.RData')
rand_corr<-rand_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                     b=fitgev[[1]]$estimate[['shape']]),
               by=group][,cor(b,a)]
rand_stats<-NULL
gc()

#Test with altered group sizes but same threshold selection mechanism and underlying parameters
set.seed(57)
cor_test_data<-generate_grouped_data(2000000,50,rexp,num_cores=40)
setnames(cor_test_data,'group','block')
cor_test_data<-cor_test_data[order(block)]
cor_test_data[,group:=ceiling(block/200)]
cor_test_data[,rand:=rand-log(50)]
setkey(cor_test_data,'block')

cor_test_stats_groupsize<-get_parallel_fits(cor_test_data[,list(block,group,identifier='rand_stats',sample=rand)],
                              get_both_fits_and_profiles,num_cores=20)
save(cor_test_stats_groupsize,file='~/dissertation/cor_test_stats_groupsize.RData')
group_size_cor = cor_test_stats_groupsize[,list(a=fitgpd[[1]]$estimate[['shape']],
                            b=fitgev[[1]]$estimate[['shape']]),
                      by=group][,cor(b,a)]
cor_test_stats_groupsize<-NULL
gc()


#Test with altered loc
set.seed(57)
cor_test_data<-generate_grouped_data(1000000,100,rexp,num_cores=40)
setnames(cor_test_data,'group','block')
cor_test_data<-cor_test_data[order(block)]
cor_test_data[,group:=ceiling(block/100)]
cor_test_data[,rand:=10+rand-log(100)]
setkey(cor_test_data,'block')

cor_test_stats_loc<-get_parallel_fits(cor_test_data[,list(block,group,identifier='rand_stats',sample=rand)],
                                  get_both_fits_and_profiles,num_cores=20)
save(cor_test_stats_loc,file='~/dissertation/cor_test_stats_loc.RData')
cor_test_loc = cor_test_stats_loc[,list(a=fitgpd[[1]]$estimate[['shape']],
                                    b=fitgev[[1]]$estimate[['shape']]),
                              by=group][,cor(b,a)]
cor_test_stats_loc<-NULL
gc()


#Test with altered scale
set.seed(57)
cor_test_data<-generate_grouped_data(1000000,100,rexp,num_cores=40)
setnames(cor_test_data,'group','block')
cor_test_data<-cor_test_data[order(block)]
cor_test_data[,group:=ceiling(block/100)]
cor_test_data[,rand:=(rand-log(100))*50]
setkey(cor_test_data,'block')

cor_test_stats_scale<-get_parallel_fits(cor_test_data[,list(block,group,identifier='rand_stats',sample=rand)],
                                      get_both_fits_and_profiles,num_cores=20)
save(cor_test_stats_scale,file='~/dissertation/cor_test_stats_scale.RData')
cor_test_scale = cor_test_stats_scale[,list(a=fitgpd[[1]]$estimate[['shape']],
                                    b=fitgev[[1]]$estimate[['shape']]),
                              by=group][,cor(b,a)]
cor_test_stats_scale<-NULL
gc()

#Test with altered shape
set.seed(57)
cor_test_data<-generate_grouped_data(1000000,100,rnorm,num_cores=40)
setnames(cor_test_data,'group','block')
cor_test_data<-cor_test_data[order(block)]
cor_test_data[,group:=ceiling(block/100)]
setkey(cor_test_data,'block')

cor_test_stats_shape<-get_parallel_fits(cor_test_data[,list(block,group,identifier='rand_stats',sample=rand)],
                                        get_both_fits_and_profiles,num_cores=20)
save(cor_test_stats_shape,file='~/dissertation/cor_test_stats_shape.RData')
cor_test_shape = cor_test_stats_shape[,list(a=fitgpd[[1]]$estimate[['shape']],
                                            b=fitgev[[1]]$estimate[['shape']]),
                                      by=group][,cor(b,a)]
cor_test_stats_shape<-NULL
gc()

#Test with altered threshold
set.seed(57)
cor_test_data<-generate_grouped_data(1000000,100,rexp,num_cores=40)
setnames(cor_test_data,'group','block')
cor_test_data<-cor_test_data[order(block)]
cor_test_data[,group:=ceiling(block/100)]
cor_test_data[,rand:=rand-log(100)]
setkey(cor_test_data,'block')

cor_test_stats_thresh<-get_parallel_fits(cor_test_data[,list(block,group,identifier='rand_stats',sample=rand)],
                                        get_both_fits_and_profiles,num_cores=20,threshold_quantile=0.9)
save(cor_test_stats_thresh,file='~/dissertation/cor_test_stats_thresh.RData')
cor_test_thresh = cor_test_stats_thresh[,list(a=fitgpd[[1]]$estimate[['shape']],
                                            b=fitgev[[1]]$estimate[['shape']]),
                                      by=group][,cor(b,a)]
cor_test_stats_thresh<-NULL
gc()
