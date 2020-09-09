
options(scipen=999)
range=seq(-6,6,0.01)
dist=rbindlist(list(data.table(x=range,
                               Density=dgev(range,
                                            loc=0,
                                            2,
                                            0.5),
                               Distribution = 'Fréchet (ξ=0.5)'),
                    data.table(x=range,
                               Density=dgev(range,
                                      loc=0,
                                      2,
                                      0),
                               Distribution = 'Gumbel (ξ=0)'),
                    data.table(x=range,
                               Density=dgev(range,
                                      loc=0,
                                      2,
                                      -0.5),
                               Distribution = 'Reverse Weibull (ξ=-0.5)')
))

ggplot(dist,aes(x=x))+geom_line(aes(y=Density,col=Distribution),alpha=0.9)+
  geom_vline(xintercept = -4,col='red',linetype="dashed")+
  geom_vline(xintercept = 4,col='blue',linetype="dashed")+
  xlim(-6,6)+
  theme_classic()+
  scale_y_continuous(name = 'Density',breaks=NULL)
ggsave("~/3_density_dist.png")

####Generate normalisation sequence plot####
set.seed(57)
data=generate_grouped_data(10000,100000,rexp,num_cores=40)
data[,in_group_index:=1:.N,by=group]
data[,block:=ceiling(in_group_index/100)]
data[,w_1:=cos(2*pi*(in_group_index-1)/100)]
data[,w_2:=cos(2*pi*(in_group_index-1)/1000)]
data[,w_constant:=w_1+w_2]
data[,very_broken_rand:=rand+1*w_constant]
setkey(data,'group')
grouped_maximums<-data[,list(M_10=max(rand[in_group_index<=10]),
                             M_100=max(rand[in_group_index<=100]),
                             M_1000=max(rand[in_group_index<=1000]),
                             M_10000=max(rand[in_group_index<=10000]),
                             M_100000=max(rand[in_group_index<=100000])
)
,by=group]


grouped_maximums<-melt(grouped_maximums[,!c("group"),with=F])
grouped_maximums[,int_n:=as.integer(gsub('M_','',variable))]
grouped_maximums[,n:=as.factor(int_n)]


ggplot(grouped_maximums,aes(x=value,col=n,group=n))+
  geom_density()+
  xlab(expression("M"[n]))+
  scale_y_continuous(name = 'Density', limits = c(0,0.5),breaks=NULL)+
  theme_classic()
ggsave("~/exp_plots_1.png")

ggplot(grouped_maximums,aes(x=(value-log(int_n)),col=n,group=n))+
  geom_density()+
  geom_line(data=data.table(x=seq(-3,5,0.01),
                            y=dgev(seq(-3,5,0.01),loc = 0,scale=1,shape=0)),
            aes(x=x,y=y,group='GEV',linetype='GEV'),col='black',linetype='dashed')+
  scale_x_continuous(name=expression("M"[n]^"*"),limits=c(-3,5))+
  scale_y_continuous(name = 'Density', limits = c(0,0.5),breaks=NULL)+
  theme_classic()
ggsave("~/exp_plots_2.png")



broken_grouped_maximums<-rbind(data[in_group_index<100000,
                                    list(position=position[which.max(very_broken_rand)],
                                         max=max(very_broken_rand)-log(sum(.N)),
                                         n=100000
                                    ),by=group],
                               data[in_group_index<10000,
                                    list(position=position[which.max(very_broken_rand)],
                                         max=max(very_broken_rand)-log(sum(.N)),
                                         n=10000
                                    ),by=group],
                               data[in_group_index<1000,
                                    list(position=position[which.max(very_broken_rand)],
                                         max=max(very_broken_rand)-log(sum(.N)),
                                         n=1000
                                    ),by=group],
                               data[in_group_index<100,
                                    list(position=position[which.max(very_broken_rand)],
                                         max=max(very_broken_rand)-log(sum(.N)),
                                         n=100
                                    ),by=group],
                               data[in_group_index<10,
                                    list(position=position[which.max(very_broken_rand)],
                                         max=max(very_broken_rand)-log(sum(.N)),
                                         n=10
                                    ),by=group]
)



ggplot(broken_grouped_maximums[,list(max,n=as.factor(n))],aes(x=max,col=n))+
  geom_density()+
  geom_line(data=data.table(x=seq(-3,10,0.01),
                            y=dgev(seq(-3,10,0.01),loc = 0,scale=1,shape=0)),
            aes(x=x,y=y),col='black',linetype='dashed')+
  geom_line(data=data.table(x=seq(-3,10,0.01),
                            y=dgev(seq(-3,10,0.01),loc = 2,scale=1,shape=0)),
            aes(x=x,y=y),col='black',linetype='dashed')+
  scale_x_continuous(name=expression("M"[n]^"*"),limits=c(-3,10))+
  scale_y_continuous(name = 'Density', limits = c(0,0.5),breaks=NULL)+
  theme_classic()
ggsave("~/e_d_normalised.png")

set.seed(102)
data_2=generate_grouped_data(10000,100000,rnorm,num_cores=40)
data_2[,in_group_index:=1:.N,by=group]
setkey(data_2,'group')
grouped_maximums_2<-data_2[,list(M_10=max(rand[in_group_index<=10]),
                                 M_100=max(rand[in_group_index<=100]),
                                 M_1000=max(rand[in_group_index<=1000]),
                                 M_10000=max(rand[in_group_index<=10000]),
                                 M_100000=max(rand[in_group_index<=100000]))
                           ,by=group]
grouped_maximums_2<-melt(grouped_maximums_2[,!c("group"),with=F])
grouped_maximums_2[variable=='M_10',n:=10]
grouped_maximums_2[variable=='M_100',n:=100]
grouped_maximums_2[variable=='M_1000',n:=1000]
grouped_maximums_2[variable=='M_10000',n:=10000]
grouped_maximums_2[variable=='M_100000',n:=100000]
grouped_maximums_2[,n:=as.factor(n)]

ggplot(grouped_maximums_2,aes(x=value,col=n,group=n))+
  geom_density()+
  xlab(expression("M"[n]))+
  scale_y_continuous(name = 'Density', limits = c(0,1.8),breaks=NULL)+
  theme_classic()
ggsave("~/a_b_plots.png")


set.seed(57)

data_3=generate_grouped_data(1000000,100,rexp,num_cores=40)
setnames(data_3,'group','block')
data_3<-data_3[order(block)]
data_3[,group:=ceiling(block/100)]
data_3[,in_group_index:=(.I-1)%%10000]
data_3[,rand:=rand-log(100)]
setkey(data_3,'block')
threshold<-data_3[,quantile(rand,0.95)]

data_3[,w_1:=cos(2*pi*(in_group_index)/100)]
data_3[,w_2:=cos(2*pi*(in_group_index)/1000)]

data_3[,w_constant:=w_1+w_2]
data_3[,broken_rand:=rand+0.5*w_constant]
broken_threshold<-data_3[,quantile(broken_rand,0.95)]
data_3[,very_broken_rand:=rand+1*w_constant]
very_broken_threshold<-data_3[,quantile(very_broken_rand,0.95)]



block_maxima<-data_3[,list(group=group[1],
                           rand=max(rand),
                           broken_rand=max(broken_rand),
                           very_broken_rand=max(very_broken_rand)),by=block]

##Generate models

rand_model=fgev(block_maxima[1:55000,rand])
rand_pot_model=fpot(data_3[block<=55000,rand],
                    threshold=data_3[block<=55000,quantile(rand,0.95)])
rand_loc=round(rand_model$estimate[['loc']],2)
rand_shape=round(rand_model$estimate[['shape']],2)
rand_scale=round(rand_model$estimate[['scale']],2)

b_model=fgev(block_maxima[1:55000,broken_rand])
b_pot_model=fpot(data_3[block<=55000,broken_rand],
                  threshold=broken_threshold)
b_loc=round(b_model$estimate[['loc']],2)
b_shape=round(b_model$estimate[['shape']],2)
b_scale=round(b_model$estimate[['scale']],2)

vb_model=fgev(block_maxima[1:55000,very_broken_rand])
vb_pot_model=fpot(data_3[block<=55000,very_broken_rand],
                    threshold=very_broken_threshold)
vb_loc=round(vb_model$estimate[['loc']],2)
vb_shape=round(vb_model$estimate[['shape']],2)
vb_scale=round(vb_model$estimate[['scale']],2)

####Initial plots####

to_plot_1<-rbind(data_3[rand>threshold,list(x=rand-threshold,
                                            Distribution='Exp(1)',
                                            threshold=threshold)],
                 data_3[broken_rand>broken_threshold,list(x=broken_rand-broken_threshold,
                                                          Distribution='w=0.5',
                                                          threshold=broken_threshold)],
                 data_3[very_broken_rand>very_broken_threshold,list(x=very_broken_rand-very_broken_threshold,
                                                                    Distribution='w=1',
                                                                    threshold=very_broken_threshold)])

ggplot()+geom_density(data=to_plot_1,
                      aes(x=x,y=..density..,col=Distribution),bounds=c(0,Inf))+
  geom_line(data=data.table(x=seq(0.01,20,0.01),
                            y=dgpd(seq(0.01,20,0.01),loc = 0,scale=1,shape=0)),
            aes(x=x,y=y,group='GPD',linetype='GPD'),col='red',linetype='dashed')+
  scale_y_continuous(name = 'Density',breaks=NULL)+
  theme_classic()+xlab(expression("X"))+
  ggtitle("Densities of theoretical and observed distributions.")


to_plot_2<-rbind(data_3[,list(x=max(rand),
                              Type='Observed',
                              Distribution='w=0'),
                        by=block],
                 data_3[,list(x=max(broken_rand),
                              Type='Observed',
                              Distribution='w=0.5'),
                        by=block],
                 data_3[,list(x=max(very_broken_rand),
                              Type='Observed',
                              Distribution='w=1'),
                        by=block])

to_plot_expected<-rbind(
                 data.table(x=seq(-5,10,0.01),
                            y=dgev(seq(-5,10,0.01),
                                   loc = 0,scale=1,shape=0),
                            Distribution='w=0',
                            Type='Estimated'),
                 data.table(x=seq(-5,10,0.01),
                            y=dgev(seq(-5,10,0.01),
                                   loc = b_loc,scale=b_scale,shape=b_shape),
                            Distribution='w=0.5',
                            Type='Estimated'),
                 data.table(x=seq(-5,10,0.01),
                            y=dgev(seq(-5,10,0.01),
                                   loc = vb_loc,scale=vb_scale,shape=vb_shape),
                            Distribution='w=1',
                            Type='Estimated'))

ggplot()+geom_density(data=to_plot_2,
                      aes(x=x,y=..density..,col=Distribution, linetype=Type),fill='white',alpha=0)+
  theme_classic()+
  geom_line(data=to_plot_expected,
            aes(x=x,y=y,col=Distribution,linetype=Type))+
  xlab(expression("X"))+
  ylab(expression("Density"))+
scale_linetype_manual(values=c("dashed", "solid"))
ggsave('~/e_d_comparison.png')

to_plot_3<-rbind(data_3[block<=20,list(i=.I,X=rand,
                                       Distribution='Standard'),
                        by=block],
                 data_3[block<=20,list(i=.I,X=broken_rand,
                                       Distribution='w=0.5'),
                        by=block],
                 data_3[block<=20,list(i=.I,X=very_broken_rand,
                                       Distribution='w=1'),
                        by=block])

ggplot(to_plot_3,
       aes(x=i,y=X,col=Distribution))+geom_line(alpha=0.3)+
  geom_smooth(se = F,method='loess',span=0.1)+theme_classic()
ggsave('~/In-group-dependence.png')


to_plot_4<-rbind(data_3[block<=20,list(i=.I,X=rand,
                                       Distribution='Standard')][X>quantile(X,0.95)],
                 data_3[block<=20,
                        list(i=.I,X=broken_rand,
                             Distribution='w=0.5')][X>quantile(X,0.95)],
                 data_3[block<=20,
                        list(i=.I,X=very_broken_rand,
                             Distribution='w=1')][X>quantile(X,0.95)])

ggplot(to_plot_4,
       aes(x=i,y=X,col=Distribution))+geom_point(alpha=0.6)+
  theme_classic()+
  ggtitle("Above threshold sampled values  Standard and\nTime Dependent distributions by index (i)")
#100 blocks per group for multi-sample runs

ggplot(block_maxima[block<=100],aes(x=block,y=very_broken_rand))+geom_point()+
  theme_classic()+
  ylab('Maximum value')+
  ggtitle("Sampled values of maxima from the Time Dependent distribution (w=1, block size=100)")
####Parallel Fit getting####
if(FALSE){
  rand_stats<-get_parallel_fits(data_3[,list(block,group,identifier='rand_stats',sample=rand)],
                                get_both_fits_and_profiles,num_cores=25)
  save(rand_stats,file='~/dissertation/rand_stats.RData')
  rand_stats<-NULL
}else{
  load(file='~/dissertation/rand_stats.RData')
}

if(FALSE){
  broken_stats<-get_parallel_fits(data_3[,list(block,group,identifier='broken_stats',sample=broken_rand)],
                                get_both_fits_and_profiles,num_cores=25)
  save(broken_stats,file='~/dissertation/broken_stats.Rdata')
  broken_stats<-NULL
}else{
  load(file='~/dissertation/broken_stats.Rdata')
}

if(FALSE){
  very_broken_stats<-get_parallel_fits(data_3[,list(block,group,identifier='very_broken_stats',
                                                    sample=very_broken_rand)],
                                  get_both_fits_and_profiles,num_cores=25)
  save(very_broken_stats,file='~/dissertation/very_broken_stats.Rdata')
  very_broken_stats<-NULL
}else{
  load(file='~/dissertation/very_broken_stats.Rdata')
}

rand_stats_fixed_shape<-get_parallel_fits(data_3[,list(block,group,identifier='rand_stats',sample=rand)],
                              get_both_fits_shape_0,num_cores=5)
broken_stats_fixed_shape<-get_parallel_fits(data_3[,list(block,group,identifier='rand_stats',sample=broken_rand)],
                                          get_both_fits_shape_0,num_cores=5)
very_broken_stats_fixed_shape<-get_parallel_fits(data_3[,list(block,group,identifier='rand_stats',sample=very_broken_rand)],
                                          get_both_fits_shape_0,num_cores=5)
####AGREEMENTS####
rand_stats<-assign_stats_from_fits(rand_stats)
get_estimate_agreements(rand_stats)
broken_stats<-assign_stats_from_fits(broken_stats)
get_estimate_agreements(broken_stats)
very_broken_stats<-assign_stats_from_fits(very_broken_stats)
get_estimate_agreements(very_broken_stats)

rand_stats[,gev_gpd_scale:=fitgev[[1]]$estimate[['scale']]+
             fitgev[[1]]$estimate[['shape']]*
             (fitgpd[[1]]$threshold-fitgev[[1]]$estimate[['loc']]),by=group]

multiplerror(x,y,sex,sey){
  sqrt(sex/s+sey/y)*(x+y)
}

rand_stats[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
                  by=group]
rand_stats[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)),
                  by=group]
broken_stats[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
                  by=group]
broken_stats[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)),
                  by=group]
very_broken_stats[,shape_diff:=fitgev[[1]]$estimate[['shape']]-fitgpd[[1]]$estimate[['shape']],
           by=group]
very_broken_stats[,shape_diff_err:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+(fitgev[[1]]$std.err[['shape']]^2)),
           by=group]

ggplot(rand_stats,aes(x=shape_diff))+geom_density()+
  geom_line(data=data.table(x=seq(-0.25,0.25,0.01),
                            y=dnorm(seq(-0.25,0.25,0.01),
                                    mean=0,
                                    sd=mean(rand_stats[,shape_diff_err_2]))),
            aes(x=x,y=y),col='blue',alpha=0.5)+
  geom_line(data=data.table(x=seq(-0.25,0.25,0.01),
                            y=dnorm(seq(-0.25,0.25,0.01),
                                    mean=0,
                                    sd=sd(rand_stats[,shape_diff]))),
            aes(x=x,y=y),col='green',alpha=0.5)+
  geom_line(data=data.table(x=seq(-0.25,0.25,0.01),
                            y=dnorm(seq(-0.25,0.25,0.01),
                                    mean=0,
                                    sd=sqrt(
                                      sd(rand_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)^2+
                                      sd(rand_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1)^2-
                                      rand_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                                                       b=fitgev[[1]]$estimate[['shape']]),
                                                 by=group][,2*cov(b,a)]))),
            aes(x=x,y=y),col='red',alpha=0.5)

rand_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err))]
broken_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err))]
very_broken_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err,shape_diff+1.96*shape_diff_err))]

#With sample SD and Cov 
rand_cov = rand_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                           b=fitgev[[1]]$estimate[['shape']]),
                     by=group][,cov(b,a)]
rand_cor = rand_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                            b=fitgev[[1]]$estimate[['shape']]),
                      by=group][,cor(b,a)]

rand_sd=sqrt(
  sd(rand_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)^2+
    sd(rand_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1)^2-
    2*rand_cov)
rand_sd_2=sqrt(
     sd(rand_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)^2+
         sd(rand_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1)^2-
         2*rand_cor*sd(rand_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)*
       sd(rand_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1))


broken_cov = broken_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                                b=fitgev[[1]]$estimate[['shape']]),
                          by=group][,cov(b,a)]
broken_sd=sqrt(
  sd(broken_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)^2+
    sd(broken_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1)^2-
    2*broken_cov)

vb_cov=very_broken_stats[,list(a=fitgpd[[1]]$estimate[['shape']],
                               b=fitgev[[1]]$estimate[['shape']]),
                         by=group][,cov(b,a)]
very_broken_sd=sqrt(
  sd(very_broken_stats[,fitgev[[1]]$estimate[['shape']],by=group]$V1)^2+
    sd(very_broken_stats[,fitgpd[[1]]$estimate[['shape']],by=group]$V1)^2-
    2*vb_cov)

rand_stats[,1-mean(inrange(0,shape_diff-1.96*rand_sd,shape_diff+1.96*rand_sd))]
broken_stats[,1-mean(inrange(0,shape_diff-1.96*rand_sd,shape_diff+1.96*rand_sd))]
very_broken_stats[,1-mean(inrange(0,shape_diff-1.96*rand_sd,shape_diff+1.96*rand_sd))]


ggplot(rbind(rand_stats[,list(Distribution="w=0",shape_diff)],
             broken_stats[,list(Distribution="w=0.5",shape_diff)],
             very_broken_stats[,list(Distribution="w=1",shape_diff)]),
       aes(x=shape_diff,col=Distribution))+geom_density()

rand_stats[,shape_diff_err_2:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+
                                   (fitgev[[1]]$std.err[['shape']]^2)-(2*
             rand_cor*fitgpd[[1]]$std.err[['shape']]*fitgev[[1]]$std.err[['shape']])),
           by=group]
broken_stats[,shape_diff_err_2:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+
                                     (fitgev[[1]]$std.err[['shape']]^2)-2*
               rand_cor*fitgpd[[1]]$std.err[['shape']]*fitgev[[1]]$std.err[['shape']]),
             by=group]
very_broken_stats[,shape_diff_err_2:=sqrt((fitgpd[[1]]$std.err[['shape']]^2)+
                                          (fitgev[[1]]$std.err[['shape']]^2)-2*
                rand_cor*fitgpd[[1]]$std.err[['shape']]*fitgev[[1]]$std.err[['shape']]),
              by=group]

rand_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err_2,shape_diff+1.96*shape_diff_err_2))]
broken_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err_2,shape_diff+1.96*shape_diff_err_2))]
very_broken_stats[,1-mean(inrange(0,shape_diff-1.96*shape_diff_err_2,shape_diff+1.96*shape_diff_err_2))]

corr_check<-merge(rand_stats[,list(I=1,gpdshape=fitgpd[[1]]$std.err[['shape']]),by=group],
             rand_stats[,list(I=1,gevshape=fitgev[[1]]$std.err[['shape']]),by=group],
             by="I",allow.cartesian=T)
corr_check[,shape_diff:=gevshape-gpdshape]
corr_check[,cor(gevshape,gpdshape)]
corr_check[,sd(shape_diff)]

###AD REJECTION_RATES####
get_specified_ad_agreements<-function(data,gev_shape=NULL,
                                      gev_loc=NULL,
                                      gev_scale=NULL,
                                      gpd_shape=NULL,
                                      gpd_scale=NULL){
  data[,list(gev_AD=get_a_d_stat(fitgev[[1]]$data,
                                 pgev,
                                 loc=ifelse(is.null(gev_loc),fitgev[[1]]$estimate[['loc']],gev_loc),
                                 shape=ifelse(is.null(gev_shape),fitgev[[1]]$estimate[['shape']],gev_shape),
                                 scale=ifelse(is.null(gev_scale),fitgev[[1]]$estimate[['scale']],gev_scale)),
             gpd_AD=get_a_d_stat(fitgpd[[1]]$data[fitgpd[[1]]$data>fitgpd[[1]]$threshold],
                                 pgpd,
                                 loc=fitgpd[[1]]$threshold,
                                 shape=ifelse(is.null(gpd_shape),fitgpd[[1]]$estimate[['shape']],gpd_shape),
                                 scale=ifelse(is.null(gpd_scale),fitgpd[[1]]$estimate[['scale']],gpd_scale)),
             gev_UMAD=get_m_a_d_u_stat(fitgev[[1]]$data,
                                       pgev,
                                       loc=ifelse(is.null(gev_loc),fitgev[[1]]$estimate[['loc']],gev_loc),
                                       shape=ifelse(is.null(gev_shape),fitgev[[1]]$estimate[['shape']],gev_shape),
                                       scale=ifelse(is.null(gev_scale),fitgev[[1]]$estimate[['scale']],gev_scale)),
             gpd_UMAD=get_m_a_d_u_stat(fitgpd[[1]]$data[fitgpd[[1]]$data>fitgpd[[1]]$threshold],
                                       pgpd,
                                       loc=fitgpd[[1]]$threshold,
                                       shape=ifelse(is.null(gpd_shape),fitgpd[[1]]$estimate[['shape']],gpd_shape),
                                       scale=ifelse(is.null(gpd_scale),fitgpd[[1]]$estimate[['scale']],gpd_scale))
  ),by=group][,list(gev_AD=mean(gev_AD>1.93),
                   gev_UMAD=mean(gev_UMAD>1),
                   gpd_AD=mean(gpd_AD>1.93),
                   gpd_UMAD=mean(gpd_UMAD>1))]
}
get_specified_ad_agreements(rand_stats)
get_specified_ad_agreements(broken_stats)
get_specified_ad_agreements(very_broken_stats)

get_specified_ad_agreements(rand_stats_fixed_shape,gev_shape=0,gpd_shape=0)
get_specified_ad_agreements(broken_stats_fixed_shape,gev_shape=0,gpd_shape=0)
get_specified_ad_agreements(very_broken_stats_fixed_shape,gev_shape=0,gpd_shape=0)

get_specified_ad_agreements(rand_stats,gev_shape=0,gev_loc=0,gev_scale=1,gpd_shape = 0,gpd_scale=1)
get_specified_ad_agreements(broken_stats,
                            gev_shape=b_model$estimate[['shape']],
                            gev_scale=b_model$estimate[['scale']],
                            gev_loc=b_model$estimate[['loc']],
                            gpd_shape=0,
                            gpd_scale=1)
get_specified_ad_agreements(very_broken_stats,
                            gev_shape=vb_model$estimate[['shape']],
                            gev_scale=vb_model$estimate[['scale']],
                            gev_loc=vb_model$estimate[['loc']],
                            gpd_shape=0,
                            gpd_scale=1)

get_specified_ad_agreements(rand_stats_fixed_shape,gev_shape=0,gev_loc=0,gev_scale=1,gpd_shape = 0,gpd_scale=1)

#TODO rejection rates with known shape. (NEEDS YET ANOTHER FIT ADDING)

####Range comparisons####

range_comparisons=confidence_range_comparison_plots(rand_stats,broken_stats,very_broken_stats)
upper_comparisons=confidence_range_comparison_plots(rand_stats,broken_stats,very_broken_stats,lower=0,prefix = "Upper end of")
lower_comparisons=confidence_range_comparison_plots(rand_stats,broken_stats,very_broken_stats,upper=0,lower=-1,prefix="Lower end of")

for(name in names(range_comparisons)){
  ggsave(filename=paste0('~/range_',name,'.png'),plot = range_comparisons[[name]])
  ggsave(filename=paste0('~/upper_',name,'.png'),plot = upper_comparisons[[name]])
  ggsave(filename=paste0('~/lower_',name,'.png'),plot = lower_comparisons[[name]])
}

####Return level comparisons####
#get actual returns

block_return_levels<-block_maxima[,list(rand_return_level_10=quantile(rand,0.9),
                   rand_return_level_100=quantile(rand,0.99),
                   rand_return_level_1000=quantile(rand,0.999),
                   broken_return_level_10=quantile(broken_rand,0.9),
                   broken_return_level_100=quantile(broken_rand,0.99),
                   broken_return_level_1000=quantile(broken_rand,0.999),
                   very_broken_return_level_10=quantile(very_broken_rand,0.9),
                   very_broken_return_level_100=quantile(very_broken_rand,0.99),
                   very_broken_return_level_1000=quantile(very_broken_rand,0.999))]

observation_return_levels<-data_3[,list(return_level_10=quantile(rand,0.999),
             return_level_100=quantile(rand,0.9999),
             return_level_1000=quantile(rand,0.99999),
             broken_return_level_10=quantile(broken_rand,0.999),
             broken_return_level_100=quantile(broken_rand,0.9999),
             broken_return_level_1000=quantile(broken_rand,0.99999),
             very_broken_return_level_10=quantile(very_broken_rand,0.999),
             very_broken_return_level_100=quantile(very_broken_rand,0.9999),
             very_broken_return_level_1000=quantile(very_broken_rand,0.99999))]

osc_block_return_levels<-block_maxima[,list(rand_return_level_10=quantile(rand,0.9),
                                        rand_return_level_100=quantile(rand,0.99),
                                        rand_return_level_1000=quantile(rand,0.999),
                                        broken_return_level_10=quantile(broken_rand,0.9),
                                        broken_return_level_100=quantile(broken_rand,0.99),
                                        broken_return_level_1000=quantile(broken_rand,0.999),
                                        very_broken_return_level_10=quantile(very_broken_rand,0.9),
                                        very_broken_return_level_100=quantile(very_broken_rand,0.99),
                                        very_broken_return_level_1000=quantile(very_broken_rand,0.999)),
                                  by=list(offset=(block-1)%%10)]

osc_obs_return_levels<-data_3[,list(.N,rand_return_level_10=quantile(rand,0.999),
                                            rand_return_level_100=quantile(rand,0.9999),
                                            rand_return_level_1000=quantile(rand,0.99999),
                                            broken_return_level_10=quantile(broken_rand,0.999),
                                            broken_return_level_100=quantile(broken_rand,0.9999),
                                            broken_return_level_1000=quantile(broken_rand,0.99999),
                                            very_broken_return_level_10=quantile(very_broken_rand,0.999),
                                            very_broken_return_level_100=quantile(very_broken_rand,0.9999),
                                            very_broken_return_level_1000=quantile(very_broken_rand,0.99999)),
                                      by=list(offset=in_group_index%%1000)]

get_return_violations<-function(stats,actual_return){
  stats[,list(I=1,Distribution='w=0',
                          retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                          retgev=fit_10gev[[1]]$estimate[['quantile']],
                          gev_conf=gev_10_confint[[1]][['quantile','97.5 %']],
                          gpd_conf=gpd_10_confint[[1]][['rlevel','97.5 %']],
                          gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                          gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                    by=group][,list(retgpd=mean(retgpd>actual_return,na.rm=T),
                                    retgev=mean(retgev>actual_return),
                                    gev_conf=mean(gev_conf>actual_return),
                                    gpd_conf=mean(gpd_conf>actual_return,na.rm=T),
                                    gev_prof=mean(gev_prof>actual_return),
                                    gpd_prof=mean(gpd_prof>actual_return,na.rm=TRUE))]
  
}

get_return_violations(rand_stats,observation_return_levels$return_level_10)
get_return_violations(broken_stats,observation_return_levels$broken_return_level_10)
get_return_violations(very_broken_stats,observation_return_levels$very_broken_return_level_10)

toplot<-rbind(osc_obs_return_levels[,list(Distribution='w=0', i=offset,
                                          ret_level=rand_return_level_10)],
              osc_obs_return_levels[,list(Distribution='w=0.5', i=offset,
                                          ret_level=broken_return_level_10)],
              osc_obs_return_levels[,list(Distribution='w=1', i=offset,
                                          ret_level=very_broken_return_level_10)])

ggplot(toplot,aes(x=i, y=ret_level, col=Distribution))+
  geom_point(size=0.5)+#geom_smooth(method='loess',span=0.1)+
  geom_hline(yintercept=observation_return_levels$return_level_10,col='red')+
  geom_hline(yintercept=observation_return_levels$broken_return_level_10,col='green')+
  geom_hline(yintercept=observation_return_levels$very_broken_return_level_10,col='blue')+
  ylab(expression("z"["0.001"]))+theme_classic()
ggsave('~/oscillating_returns.png')


ggplot(osc_obs_return_levels,aes(x=offset,y=very_broken_return_level_10))+
  geom_point()+geom_smooth(method='loess',span=0.1)+
  geom_point(aes(x=offset,y=very_broken_return_level_100))+
  geom_hline(yintercept=observation_return_levels$return_level_10,col='green')+
  geom_hline(yintercept=observation_return_levels$very_broken_return_level_10,col='red')

####returns vs models 1####
test<-rbind(merge(rand_stats[,list(I=1,Distribution='w=0',
                                   retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                   retgev=fit_10gev[[1]]$estimate[['quantile']],
                                  gev_conf=gev_10_confint[[1]][['quantile','97.5 %']],
                                  gpd_conf=gpd_10_confint[[1]][['rlevel','97.5 %']],
                                  gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                  gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                             by=group],
              osc_obs_return_levels[,list(I=1,offset,
                                          actual_return=rand_return_level_10)],
              by='I',allow.cartesian=T),
            merge(broken_stats[,list(I=1,Distribution='w=0.5',
                                   retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                   retgev=fit_10gev[[1]]$estimate[['quantile']],
                                   gev_conf=gev_10_confint[[1]][['quantile','97.5 %']],
                                   gpd_conf=gpd_10_confint[[1]][['rlevel','97.5 %']],
                                   gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                   gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                             by=group],
                  osc_obs_return_levels[,list(I=1,offset,
                                              actual_return=broken_return_level_10)],
                  by='I',allow.cartesian=T),
            merge(very_broken_stats[,list(I=1,Distribution='w=1',
                                   retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                   retgev=fit_10gev[[1]]$estimate[['quantile']],
                                   gev_conf=gev_10_confint[[1]][['quantile','97.5 %']],
                                   gpd_conf=gpd_10_confint[[1]][['rlevel','97.5 %']],
                                   gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                   gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                             by=group],
                  osc_obs_return_levels[,list(I=1,offset,
                                              actual_return=very_broken_return_level_10)],
                  by='I',allow.cartesian=T)
)

toplot<-test[,list(gpd_ret_diff=actual_return-retgpd,
                   gev_ret_diff=actual_return-retgev,
                   gpd_ret_diff_top_bound=actual_return-gpd_conf,
                   gev_ret_diff_top_bound=actual_return-gev_conf,
                   gpd_ret_diff_top_prof=actual_return-gpd_prof,
                   gev_ret_diff_top_prof=actual_return-gev_prof,
                   gpd_ret_over_top_bound= actual_return>gpd_conf,
                   gev_ret_over_top_bound= actual_return>gev_conf,
                   gpd_ret_over_top_prof= actual_return>gpd_prof,
                   gev_ret_over_top_prof= actual_return>gev_prof),
             by=list(group,offset,Distribution)][,list(gpd_ret_diff=mean(gpd_ret_diff,na.rm=T),
                                                       gpd_ret_upper=quantile(gpd_ret_diff, 0.975,na.rm=T),
                                                       gev_ret_diff=mean(gev_ret_diff,na.rm=T),
                                                       gev_ret_upper=quantile(gev_ret_diff, 0.975,na.rm=T),
                                                       gpd_ret_diff_top_prof=mean(gpd_ret_diff_top_prof,na.rm=T),
                                                       gpd_ret_diff_top_upper=quantile(gpd_ret_diff_top_prof, 0.975,na.rm=T),
                                                       gev_ret_diff_top_prof=mean(gev_ret_diff_top_prof,na.rm=T),
                                                       gev_ret_diff_top_upper=quantile(gev_ret_diff_top_prof, 0.975,na.rm=T)
             ),
                                                 by=list(offset,Distribution)]

ggplot(toplot,aes(x=offset,y=gpd_ret_diff,col=Distribution))+
  geom_point()+xlab('i')+ylab(expression(bar("Δz"[10])))+theme_classic()
ggsave('~/return_diffs_gpd.png')
ggplot(toplot,aes(x=offset,y=gev_ret_diff,col=Distribution))+
  geom_point()+xlab('i')+ylab(expression(bar("Δz"[10])))+theme_classic()
ggsave('~/return_diffs_gev.png')
ggplot(toplot,aes(x=offset,y=gpd_ret_diff_top_prof,col=Distribution))+
  geom_point()+xlab('i')+ylab(expression(bar("Δz"[10])))+theme_classic()
ggsave('~/return_diffs_top_gpd.png')
ggplot(toplot,aes(x=offset,y=gev_ret_diff_top_prof,col=Distribution))+
  geom_point()+xlab('i')+ylab(expression(bar("Δz"[10])))+theme_classic()
ggsave('~/return_diffs_top_gev.png')


####Crazymerge time####
data_3[,offset:=in_group_index%%1000]
setkey(data_3,offset)
get_return_stats<-function(index){
  print(i)
  rand_to_use<-rand_stats[group==index,list(I=1,
                                   retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                   retgev=fit_10gev[[1]]$estimate[['quantile']],
                                   gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                   gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                            by=group]
  broken_to_use<-broken_stats[group==index,list(I=1,
                                               retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                               retgev=fit_10gev[[1]]$estimate[['quantile']],
                                               gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                               gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                             by=group]
  vb_to_use<-very_broken_stats[group==index,list(I=1,
                                               retgpd=fit_10gpd[[1]]$estimate[['rlevel']],
                                               retgev=fit_10gev[[1]]$estimate[['quantile']],
                                               gev_prof=gev_profile_10_confints[[1]][['quantile','upper']],
                                               gpd_prof=gpd_profile_10_confints[[1]][['rlevel','upper']]),
                             by=group]
  rbind(data_3[,list(retgpd=mean(rand>rand_to_use$retgpd),
                retgev=mean(rand>rand_to_use$retgev),
                gev_prof=mean(rand>rand_to_use$gev_prof),
                gpd_prof=mean(rand>rand_to_use$gpd_prof),
                group=index,
                Distribution='w=0'),
          by=offset],
        data_3[,list(retgpd=mean(broken_rand>broken_to_use$retgpd),
                retgev=mean(broken_rand>broken_to_use$retgev),
                gev_prof=mean(broken_rand>broken_to_use$gev_prof),
                gpd_prof=mean(broken_rand>broken_to_use$gpd_prof),
                group=index,
                Distribution='w=0.5'),
          by=offset],
        data_3[,list(retgpd=mean(very_broken_rand>vb_to_use$retgpd),
                retgev=mean(very_broken_rand>vb_to_use$retgev),
                gev_prof=mean(very_broken_rand>vb_to_use$gev_prof),
                gpd_prof=mean(very_broken_rand>vb_to_use$gpd_prof),
                group=index,
                Distribution='w=1'),
          by=offset])
}

results<-mclapply(1:1000,get_return_stats)
results2<-rbindlist(results)
save(results2,file='~/combined_100_returns_2.Rdata')
ggplot(results2[offset!=-1,list(gpd_prof=mean(gpd_prof)),by=list(offset,Distribution)],
       aes(x=offset,y=gpd_prof,col=Distribution))+geom_point(size=0.1)+
theme_classic()+xlab('i')+ylab(expression("p"["z0.001"]))
ggsave('~/gpd_97.png')
ggplot(results2[offset!=-1,list(gev_prof=mean(gev_prof)),by=list(offset,Distribution)],
       aes(x=offset,y=gev_prof,col=Distribution))+geom_point(size=0.1)+
  theme_classic()+xlab('i')+ylab(expression("p"["z0.001"]))
ggsave('~/gev_97.png')


