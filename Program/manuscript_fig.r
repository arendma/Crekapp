#Plot manuscript figures
library(ggplot2)
source("Program/deps/plotting/plot_nidle_kapp.r")
source("Program/deps/plotting/plotexternalcomp.r")
#library(xlsx)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#A132F0", "#0072B2", "#D55E00", "#CC79A7")
cnames=c("control", "highcell","highsalt", "hightemp", "noshaking", "UVM4", "Stop1", "Stop2", "dark")
clabs=c('Control CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690','Control UVM4', 'SDP OE1 UVM4', 'SDP OE2 UVM4',  'Dark CC1690')
names(clabs)=cnames

subsys_kapp_ndat=plot_nidle_kapp(cnames, clabs, TRUE) 

plotexternalcomp(cnames, subsys_kapp_ndat) 

plotvdis = function(cnames, resdir) {
  #Function to plot the comulative distribution of flux values for the different conditions
  #- char cnames:  character vector with column names of conditions in NIDLE output
  # - char resdir: character vector giving the path to save figures to
  source('Program/deps/utilities/expandcolv3.r', local=TRUE)
  #Import fluxdata
  fluxdat=read.delim("Results/NIDLE/kcat_n_flux.tsv")
  #regenerate color code from kapp pca plots 
  colorval=scales::hue_pal()(length(cnames)+1)[1:length(cnames)]
  names(colorval)=c(cnames)

  #remove blocked reactions
  fluxdat=fluxdat[apply(fluxdat[,cnames]==0, 1, sum)<length(cnames),]
  #add rank
  plot_cumsum = function(fluxdat, resdir, cnames) {
    #function to plot cumulative sum of fluxes for each condition ordered by the reactions in control condition
    fluxdat=fluxdat[order(fluxdat$control),]
    fluxdat$rancont=1:nrow(fluxdat)
    plotdat=fluxdat
    plotdat[,cnames]=apply(plotdat[cnames], 2, cumsum)
    plotdat=expandcol(plotdat, which(colnames(plotdat) %in% cnames))
    p1=ggplot(plotdat, aes(x=rancont, y=x1, color=x2)) + geom_line() + theme_bw() + 
      scale_color_manual(values = colorval, labels=clabs, name="Conditions") +theme(text=element_text(size=25), axis.text.x = element_blank()) + 
      scale_y_log10() + xlab("Reactions") + ylab("Cumulative flux sum") # + xlab(expression(mmol~gDW{-1}h{-1}))
    ggsave(file.path(resdir, "cumfluxdist.pdf"), p1)
  }
  plot_cumsum(fluxdat, resdir, cnames)
 
  plot_fluxdens= function(fluxdat, resdir, cnames) {
    plotdat=expandcol(fluxdat, which(colnames(fluxdat) %in% cnames))
    #remove 0 values
    plotdat=plotdat[plotdat$x1!=0,]
    #check for 0 or negative values
    if (any(plotdat$x1<=0)) {stop("values <=0 found in flux values")}
  plotdat$x1=log10(plotdat$x1)
  #convert to logarithmic 
  p=ggplot(plotdat,aes(x=x1, color=x2)) + geom_density( adjust=1.5, fill=NA) +
    scale_color_discrete(labels=clabs, name= "Conditions") + theme_bw() + 
    theme(text=element_text(size=25)) + xlab("Flux")
  ggsave(file.path(resdir, "fluxdens.pdf"), p)
  }
  plot_fluxdens(fluxdat, resdir, cnames)
}
plotvdis(cnames, "Results/NIDLE")


plot_ecmodcomp=function(cnames) {
  source('Program/deps/utilities/expandcolv3.r', local=TRUE)
  predgrowth=read.csv('Results/eccomp/chemostat_comp.txt')
  plot_predgrowth=expandcol(predgrowth, 2:ncol(predgrowth))
  plot_predgrowth$oldcols=gsub('_.*$', '', plot_predgrowth$oldcols)
  shps=c(8, 16, 16, 16,17, 17)
  names(shps)=c("exp_mu", "FBA","GKOraw", "GKOadp", "NDLraw", "NDLadp")
  
  colnames(plot_predgrowth)=c('Condition', 'mu', 'Datatype')
  plot1=ggplot(plot_predgrowth, aes(x=Condition ,y=mu, shape= Datatype,  colour=Datatype)) + geom_jitter(size=8,alpha=0.7, height=0, width=0.1) +
    theme_bw() + theme(text=element_text(size=25), legend.position="bottom") + 
    scale_colour_manual(values=c('black', 'yellow2', 'seagreen', 'skyblue4', 'seagreen2', 'lightskyblue'), labels=c('Experimental', 'FBA','GECKO adapted' , 'GECKO raw', 'GECKO adapted + NIDLE', 'GECKO raw + NIDLE')) + 
    scale_shape_manual(values=shps, labels=c('Experimental', 'FBA','GECKO adapted' , 'GECKO raw', 'GECKO adapted + NIDLE', 'GECKO raw + NIDLE'))+
    scale_y_log10() + ylab('Growth rate') +  theme(text=element_text(size=25), axis.text.x=element_text(angle=60, hjust=1), legend.position="bottom") + 
    coord_flip()
  ggsave('Results/eccomp/chemostat_comp.pdf', plot1, width=13, height=6, useDingbats=FALSE)
  spearprot=read.delim('Results/eccomp/abunvse/spearman.tsv')
  spearprot$Row=c('Control CC1690', 'Dark CC1690','UVM4', 'SDP OE1 UVM4', 'SDP OE2 UVM4')
  plot_spearprot=expandcol(spearprot, 2:ncol(spearprot))
  colnames(plot_spearprot)=c('Condition', 'sp_cor', 'Model')
  shps2=c(16,16,17,17)
  names(shps2)=c('adpGKO', 'rawGKO', 'adpGKONDL', 'rawGKONDL')
  cols2=c('seagreen',  'skyblue4', 'seagreen2',  'lightskyblue')
  names(cols2)=c('adpGKO', 'rawGKO', 'adpGKONDL', 'rawGKONDL')
  plot2=ggplot(plot_spearprot, aes(x=sp_cor ,y=Condition, color=Model, shape=Model)) + geom_jitter(size=5,alpha=0.7, height=0.1, width=0) +
    theme_bw() + theme(text=element_text(size=25), legend.position="bottom") + 
    scale_color_manual(values=cols2, labels=c("GECKO adapted", "GECKO adapted + NIDLE", "GECKO raw", "GECKO raw + NIDLE")) + 
    scale_shape_manual(values=shps2, labels=c("GECKO adapted", "GECKO adapted + NIDLE", "GECKO raw", "GECKO raw + NIDLE")) +
    xlab('Spearman correlation')
  ggsave('Results/eccomp/abunvse/spearman.pdf', plot2, width=13, height=4, useDingbats=FALSE)
}
plot_ecmodcomp(cnames)

kcat_ana=function() {
  #code to extract information on the coverage of photosynthetic organisms in BRENDA an SABIO-RK database
  #Import database
  bs_db=read.delim("Data/kcats/kcats-full-lineages.tsv")
  #get yeast ecnumbers
  n_cere=length(unique(bs_db$ec_number[grep("cerevisiae", bs_db$lineage)]))
  #get chlamy entries
  n_chlamydomonas=length(unique(bs_db$ec_number[grep("Chlamydomonas", bs_db$lineage)]))
  print(paste("Number of Chlamydomonas entries:", n_chlamydomonas))
  print(paste("Relative amount of yeast entries compared to chlamy: ", signif(n_cere/n_chlamydomonas,2)))
  #get fungi vs plant number of kcats
  n_fungi=length(unique(bs_db$ec_number[grep("Dikarya", bs_db$lineage)]))
  print(paste("Number of fungal EC NUmbers:", n_fungi))
  n_viri=length(unique(bs_db$ec_number[grep("Viridiplantae", bs_db$lineage)]))
  print(paste("Number of Viridiplantae EC Numbers:", n_viri))
  #relative number
  r_viri=length(grep("Viridiplantae", bs_db$lineage))/nrow(bs_db)
  print(paste("Relative amount of viridiplantae entries in sabio and brenda:", r_viri))
  
  r_chlamy=length(grep("Chlorophyta", bs_db$lineage))/nrow(bs_db)
  print(paste("Relative amount of chlorophyta entries in sabio and brenda:", signif(r_chlamy, 2)))
}
