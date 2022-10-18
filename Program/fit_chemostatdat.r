
#Script to fit a linear model to predict acetate upatke from mu
#Data extracted from publication Imam et al. (2015) The plant journal:
mu=c(0.082, 0.082, 0.061, 0.061)
nh4up=c(3.035, 0.748, 1.259, 1.071)
po4up=c(0.099, 0.078, 0.047, 0.041)
acup=c(1.364, 2.261, 1.059, 1.234)

saheed_dat=data.frame(mu, nh4up, po4up, acup)
mod1=lm(acup ~mu, saheed_dat)
print(summary(mod1))
print('using a linear model for mu=0.0937 the predicted uptake rate is')
#Experimental growth rates
proteomics_dat=data.frame(Condition=c('control', 'highcell', 'highsalt', 'hightemp', 'noshaking','UVM4', 'Stop1', 'Stop2'),
           mu=c(0.11, 0.02, 0.02, 0.02, 0.04, 0.0937,0.0864,0.0837))
pred=predict(mod1, proteomics_dat , se.fit=TRUE)
proteomics_dat$EX_ac_e=pred$fit+pred$se.fit
proteomics_dat=rbind(proteomics_dat, c('dark', '0.02', '1.6')) #only set light for dark to zero, otherwise take default alvue
proteomics_dat$EX_photonVis_e=c(rep(80, 8), 0)
proteomics_dat$model=c(rep('mixo', 8), 'hetero')
write.table(proteomics_dat, '/Data/QconCAT_David20220124/fitted_acup.tsv', row.names=FALSE, sep='\t')
