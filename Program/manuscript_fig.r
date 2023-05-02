#Plot manuscript figures
library(ggplot2)
#library(xlsx)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#A132F0", "#0072B2", "#D55E00", "#CC79A7")
cnames=c("control", "highcell","highsalt", "hightemp", "noshaking", "UVM4", "Stop1", "Stop2", "dark")
clabs=c('Control CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690','Control UVM4', 'SDP OE1 UVM4', 'SDP OE2 UVM4',  'Dark CC1690')
names(clabs)=cnames
plot_nidle_kcat= function(cnames, clabs, same_set) {
  #Code to plot the internal comparison of NIDLE kcats
  #Input:
  # - logical same_set: If true than in the scatter plot of nidlekapps and kcats only reactions with 
  #                     kcat values for pro & eukaryots are taken into account
  #OUTPUT:
  # - df subsys_kapp_ndat:  Data frame of raw NIDLE values with kappmax column (maximum over conditions)
  #                         and cleaned up subsystem info from universal irreversible model
  source('Program/deps/utilities/expandcolv3.r', local=TRUE)
  require(ggfortify)
  resdir='Results/NIDLE/'
  if (!(dir.exists(resdir))) {
    dir.create(resdir)
  }
  
  #Import NIDLE kcats 
  kapp_ndat=read.delim(file.path(resdir,'kcat_n.tsv'), na.strings = 'NaN')
  #check if all cnames are contained in the NIDLE data
  if (!all(cnames %in% colnames(kapp_ndat))) {
    stop("Not all cnames are contained as column names in kcat_n.tsv")
  }
  colnames(kapp_ndat)[colnames(kapp_ndat) %in% 'proc_Kcats']='prok_Kcats'
  
  kapp_stat_barplot= function() {
    #function to plot count statistics on the overall number of kapps from NIDLE 
    #INPUT:
    # - char resdir: character vector giving the path to save figures to
    
    ##plot a barplot giving the number of maximum kcats per condition
    max_idx=cnames[apply(kapp_ndat[, colnames(kapp_ndat) %in% cnames],1,which.max)]
    ggsave(file.path(resdir, 'maxpercond.pdf'), ggplot(data.frame(max_idx), aes(max_idx)) + geom_bar() + 
             scale_x_discrete(labels = clabs) + 
             labs(x='Condition', y='Count') + theme_bw() + theme(axis.text.x = element_text(angle=50, hjust=1), text=element_text(size=20), panel.grid.major.x=element_blank()), 
           height=4, width=6)
    write.table(table(cnames[max_idx]), file = file.path(resdir, 'maxpercond.tsv'), sep = "\t", row.names=FALSE)
    ##plot number of conditions in which a kapp could be calculated 
    num_kapp=apply(!(is.na(kapp_ndat[, colnames(kapp_ndat) %in% cnames])) & !(kapp_ndat[, colnames(kapp_ndat) %in% cnames]==0),1, sum)
    #save histogram 
    ggsave(file.path(resdir, 'histkapppercond.pdf'), ggplot(data.frame(num_kapp), aes(num_kapp)) + geom_bar() + 
             scale_x_continuous(breaks=1:length(cnames), labels=1:length(cnames)) +
             labs(x='No of Conditions', y='Count') + theme_bw()+ theme(text=element_text(size=20), panel.grid.major.x=element_blank()), 
           height=4, width=5.5)
    write.table(table(num_kapp), file = file.path(resdir, 'histkappercond.tsv'),sep = "\t", row.names=FALSE)
  }
  kapp_stat_barplot()
  #calculate maximum colum
  kapp_ndat$kappmax=apply(kapp_ndat[, colnames(kapp_ndat) %in% cnames],1,max, na.rm=TRUE)
  
  
  plot_man_kcat_scatter = function(same_set) {
    #TO_DO: remove pseudocount addition since only max column is used in plotting and this is garanteered to be non 0
    #function to plot scatter plot of manual queried (conservative matching) kcats vs NIDLE kapps
    #INPUT:
    # - logical same_set: If true than in the scatter plot of nidlekapps and kcats only reactions with 
    #                     kcat values for pro & eukaryots are taken into account
    warning("ATTENTION for the manual kcat scatter plot currently the minimum nonzero value times 10^-4 is added as pseudocount before log 10 transformation")
    
    kappmat=kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]
    #add pseodocount
    psco=min(kappmat[kappmat!=0 & !(is.na(kappmat))])*0.0001
    kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]+psco
    if (any(kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]<=0, na.rm = TRUE)){stop("Values<=0 in kcat data detected. Aborting log transformation...")}
    #log transfrom
    kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=log10(kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")])
    #plot pca 
    if (same_set) {
      euk_idx=!is.na(kapp_ndat$euk_Kcats)&!is.na(kapp_ndat$prok_Kcats)
      proc_idx=!is.na(kapp_ndat$euk_Kcats)&!is.na(kapp_ndat$prok_Kcats)
    } else{
      euk_idx=!is.na(kapp_ndat$euk_Kcats)
      proc_idx=!is.na(kapp_ndat$prok_Kcats)
    } #finish this for unmatched
    #euk_pearsrho=cor(x=kapp_ndat[euk_idx, 'kappmax'], y=kapp_ndat[euk_idx, 'euk_Kcats'])
    #proc_pearsrho = cor(x=kapp_ndat[proc_idx, 'kappmax'], y=kapp_ndat[proc_idx, 'prok_Kcats'])
    plotdat_euk=data.frame(reaction=kapp_ndat[euk_idx, 'Rxns'],
                           kappmax=kapp_ndat[euk_idx, 'kappmax'], 
                           Kcats=kapp_ndat[euk_idx, 'euk_Kcats'], 
                           taxon=rep('eukaryot', sum(euk_idx)))
    plotdat_proc=data.frame(reaction=kapp_ndat[proc_idx, 'Rxns'],
                            kappmax=kapp_ndat[proc_idx, 'kappmax'], 
                            Kcats=kapp_ndat[proc_idx, 'prok_Kcats'], 
                            taxon=rep('prokaryot', sum(proc_idx)))
    plotdat=rbind(plotdat_euk, plotdat_proc)
    euk_pearsrho=cor(x=plotdat$kappmax[plotdat$taxon %in% "eukaryot"], y=plotdat$Kcats[plotdat$taxon %in% "eukaryot"])
    prok_pearsrho=cor(x=plotdat$kappmax[plotdat$taxon %in% "prokaryot"], y=plotdat$Kcats[plotdat$taxon %in% "prokaryot"])
    #Plot overlay of eukaryot and prokaryot comparison 
    ggsave(file.path(resdir,'kappmaxvxeucprockcatsamdet.pdf'),ggplot(plotdat,aes(x=kappmax, y=Kcats, colour=as.factor(taxon))) +
             annotate(geom='text', x=min(plotdat$kappmax)+3, y=max(plotdat$Kcats)-1, size=5,  label=paste(' euk. \u03c1 = ', round(euk_pearsrho, 2) , '\nprok. \u03c1 = ', round(prok_pearsrho, 2))) +
             scale_color_manual('Taxon', values=cbPalette[c(1,3)],  labels=c('Eukaryots', 'Bacteria'))+geom_point(size=3) + theme_bw() + theme(text = element_text(size=20)) + labs(x= 'log(kapp)', y='literature log(kcat)'), device=grDevices::cairo_pdf, width=6, height = 5)
    viri_pearsrho=cor(x=kapp_ndat$viri_Kcats[!(is.na(kapp_ndat$viri_Kcats))], y=kapp_ndat$kappmax[!(is.na(kapp_ndat$viri_Kcats))])
    ggsave(file.path(resdir, 'loglitvsmaxkapp2det.pdf'),ggplot(kapp_ndat[!is.na(kapp_ndat$viri_Kcats),],aes(x=kappmax, y=viri_Kcats, colour=as.factor(viri_Kcat_matches))) + annotate(geom='text', x=min(kapp_ndat$kappmax[!is.na(kapp_ndat$viri_Kcats)])+6, y=max(kapp_ndat$viri_Kcats[!is.na(kapp_ndat$viri_Kcats)]), size=5,  label=paste('Pears. corr.=', round(viri_pearsrho, 2) )) +
             scale_color_discrete('Matched by', label=c('EC+substr.+taxon', 'EC+substr.', 'EC+taxon', 'EC'))+geom_point(size=3) + theme_bw() + theme(text = element_text(size=20)) + labs(x= 'log(k_app)', y='literature log(k_cat)', title='Comparison with viridiplantae Kcats'), width=6, height = 5)
  }
  plot_man_kcat_scatter(same_set)
  
  plot_count_dat = function (){
    #Plot count statistics on subsystems without analysing actual kapp values, so no log transformation here
    #OUTPUT:
    # - kapp_ndat: a kapp_ndat data frame with subsystem info added but with raw values
    
    #get subsystem/pathway info from global kcat data for irreversible model
    glob_kcats=read.delim("Data/Cre1355/Cre1355_lit_kcats.txt", na.strings='NaN')
    # if (any(glob_kcats[,c("viri_Kcats", "euk_Kcats", "proc_Kcats")]<=0, na.rm = TRUE)) {
    #   stop("Found model reactions with assigned kcat <0 log transformation impossible. Aborting...")
    # }
    # glob_kcats[,c("viri_Kcats", "euk_Kcats", "proc_Kcats")]=log10(glob_kcats[,c("viri_Kcats", "euk_Kcats", "proc_Kcats")])
    glob_kcats=merge(glob_kcats, kapp_ndat[, c(1,21)],by.x="Reaction", by.y="Rxns", all.x=TRUE)
    kapp_ndat=merge(glob_kcats[,1:2], kapp_ndat, by.x="Reaction", by.y="Rxns", all.y = TRUE)
    #plot a stacked barplot with the kcat information available from viridiplantae
    subsyslevels=unique(unlist(strsplit(glob_kcats$SubSystem, split = "\t")))
    plotdat6=lapply(subsyslevels, function(x) {
      no_kcat=sum(!is.na(glob_kcats$viri_Kcats[grepl(x, glob_kcats$SubSystem)]))
      no_kapp=sum(is.na(glob_kcats$viri_Kcats[grepl(x, glob_kcats$SubSystem)]) &
                    !is.na(glob_kcats$kappmax[grepl(x, glob_kcats$SubSystem)]))
      if (any(c(no_kcat, no_kapp)!=0)) {
        return(data.frame(subsys=rep(x,2), type=c("Viridiplantae kcat", "NIDLE kapp"),
                          value=c(no_kcat, no_kapp)))
      } else {
        return(NULL)
      }})
    plotdat6=do.call(rbind, plotdat6)
    plotdat6$subsys=factor(plotdat6$subsys, levels=plotdat6$subsys[which(plotdat6$type %in% "NIDLE kapp")[order(plotdat6$value[plotdat6$type %in% "NIDLE kapp"])]])
    ggsave(file.path(resdir,'determinedkmax.pdf'), ggplot(plotdat6, aes(x=subsys, y=value, fill=type)) + geom_bar(position = "stack", stat="identity") +
             scale_fill_manual(values=cbPalette[c(2,1)], name="Kinetic data source") + ylab('No of Reactions') + xlab("Subsystem") + theme_bw() + 
             theme(axis.text.x = element_blank(), text=element_text(size=20)), width = 12, height = 4)
    ggsave(file.path(resdir,'determinedkmaxlabels.pdf'), ggplot(plotdat6, aes(x=subsys, y=value, fill=type)) + geom_bar(position = "stack", stat="identity") +
             scale_fill_manual(values=cbPalette[c(2,1)], name="Kinetic data source") + ylab('No of Reactions') + xlab("Subsystem") + theme_bw() + 
             theme(axis.text.x = element_text(angle=60, hjust=1), text=element_text(size=20)), width = 14, height = 10)
    write.table(plotdat6, file = file.path(resdir, 'determinedkmax.tsv'),sep = "\t", row.names=FALSE)
    #trim subsystems
    kapp_ndat$SubSystem=gsub("Biosynthesis of ", "",kapp_ndat$SubSystem)
    sparseSubSystem=unlist(lapply(kapp_ndat$SubSystem, function (x) {
      res=strsplit(x, c(' |,|\\t'))[[1]]
      return(paste(res[1:min(2, length(res))], sep='' ,collapse=' '))}))
    kapp_ndat=cbind(sparseSubSystem,kapp_ndat)
    rm(sparseSubSystem)
    kapp_ndat$sparseSubSystem=gsub(' $| and$| /$', '',kapp_ndat$sparseSubSystem)
    kapp_ndat$sparseSubSystem[kapp_ndat$sparseSubSystem %in% names(table(kapp_ndat$sparseSubSystem)[table(kapp_ndat$sparseSubSystem)<5])]="Other"
    plotdat3=data.frame(SubSystem=names(table(kapp_ndat$sparseSubSystem)), No=as.vector(table(kapp_ndat$sparseSubSystem)))
    #Sort subsystems into large metabolic entities
    Macrosystem=rep("Other", length(table(kapp_ndat$sparseSubSystem)))
    Macrosystem[plotdat3$SubSystem %in% c('Alanine', 'Arginine', 'Glutamate metabolism', 'Glycine', 'Histidine metabolism', 'Lysine biosynthesis', 'Methionine metabolism', 'Phenylalanine', 'Valine')]="Amino Acid"
    Macrosystem[plotdat3$SubSystem %in% c('Carbon fixation', 'Glycolysis', 'Propanoate metabolism', 'Pyruvate metabolism')]="Central Carbon"
    Macrosystem[plotdat3$SubSystem %in% c('Butanoate metabolism', 'Fatty acid', 'Glycerolipid metabolism','Glycerophospholipid metabolism', 'steroids', 'unsaturated fatty')]='Lipid'
    Macrosystem[plotdat3$SubSystem %in% c('Carotenoid biosynthesis', 'N-Glycan biosynthesis', 'One carbon', 'Porphyrin', 'Riboflavin metabolism', 'Thiamine metabolism')]='Cofactor and Modification'
    Macrosystem[plotdat3$SubSystem %in% c('Purine metabolism', 'Pyrimidine metabolism')]='Nucleotide'
    Macrosystem[plotdat3$SubSystem %in% c('Transport')]="Transport"
    plotdat3=cbind(plotdat3, Macrosystem)
    #add macorsystem info to kapp_ndat
    Macrosystem=sapply(kapp_ndat$sparseSubSystem, function(x){return(plotdat3$Macrosystem[plotdat3$SubSystem %in% x])})
    kapp_ndat=cbind(Macrosystem, kapp_ndat)
    rm(Macrosystem)
    #final cleanup
    plotdat3$SubSystem=gsub(' metabolism$| biosynthesis$', '', plotdat3$SubSystem)
    plotdat3$SubSystem[plotdat3$SubSystem %in% 'One carbon']="One carbon pool"
    plotdat3$SubSystem[plotdat3$SubSystem %in% 'unsaturated fatty']='Unsaturated fatty acid'
    #plot number of kapps for each subsystem 
    ggsave(file.path(resdir, 'kapp_no_subsys.pdf'), ggplot(plotdat3, aes(x=SubSystem, y=No, fill=Macrosystem)) + geom_bar(stat='identity') +  scale_x_discrete(limits=plotdat3$SubSystem[order(plotdat3$No)]) +
             theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust=1), text = element_text(size=20)) + scale_fill_manual(values=cbPalette) +
             labs(y='No of kapps', x='Metabolic Subsystem'), height=9, width=9)
    write.table(plotdat3, file = file.path(resdir, 'kapp_no_subsys.tsv'),sep = "\t", row.names=FALSE)
    # plot stacked column of macrosystmes
    plotdat4=data.frame(No=as.vector(by(plotdat3$No, plotdat3$Macrosystem, sum)), Macrosystem=names(by(plotdat3$No, plotdat3$Macrosystem, sum)), 
                        x=rep('Macrosystem', length(unique(plotdat3$Macrosystem))))
    plotdat4$Macrosystem=factor(plotdat4$Macrosystem, levels=plotdat4$Macrosystem[order(plotdat4$No)])
    ggsave(file.path(resdir, 'kapp_no_macsys.pdf'), ggplot(plotdat4, aes(x=x, y=No, fill=Macrosystem)) + geom_bar(position='stack', stat='identity') +
             scale_fill_manual(limits=levels(plotdat4$Macrosystem), values=cbPalette) + ylab('No of kapps') +theme_bw() + 
             theme(axis.text.x = element_blank(),  axis.title.x = element_blank(), text=element_text(size=20)), 
           height=3, width=5)
    write.table(plotdat4, file = file.path(resdir, 'kkapp_no_macsys.tsv'),sep = "\t", row.names=FALSE)
    
    # return the kapp data frame with subsystem info added 
    return(kapp_ndat)
  }
  subsys_kapp_ndat=plot_count_dat()
  
  plot_kapp_nsparse = function(kapp_ndat) {
    #plot overviews of nonsparse (!is.na & !=0) kapp and flux data no spseudocount addition before log transformtion
    #INPUT:
    #-kapp_ndat: a kapp_ndat data frame with subsystem information as returned by the plot_count_dat function
    #
    kapp_ndat_nons=kapp_ndat[apply(is.na(kapp_ndat[, colnames(kapp_ndat) %in% cnames]) | kapp_ndat[,colnames(kapp_ndat) %in% cnames]==0, 1, sum)==0,]
    if (any(kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]<=0, na.rm = TRUE)){stop("Values<=0 in kcat data detected. Aborting log transformation...")}
    #log transfrom
    kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=log10(kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")])
    plotdat5=data.frame(m=apply(kapp_ndat_nons[colnames(kapp_ndat) %in% cnames], 1, mean), s=apply(kapp_ndat_nons[colnames(kapp_ndat) %in% cnames], 1, sd), Macrosystem=kapp_ndat_nons$Macrosystem)
    ggsave(file.path(resdir, 'nonsparce_sdvsmean.pdf'),
           ggplot(plotdat5, aes(x=m, y=s, color=Macrosystem)) + geom_point(size=4, alpha=0.7) +
             scale_color_manual(values=cbPalette) + 
             ylab('standard deviation') + xlab('mean') +
             theme_bw() + theme(text=element_text(size=25)),
           width=7, height = 5)
    #plot PCA plot of kapp conditions and kappmax
    pca1 <- prcomp(t(kapp_ndat_nons[,colnames(kapp_ndat) %in% c(cnames, "kappmax")]))
    bp1 <- autoplot(pca1, data=data.frame(sample=factor(c(cnames, "kappmax"), levels =c(cnames, "kappmax"))),colour='sample', size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') +
      scale_color_discrete(labels = c(clabs, 'maxKapp'))
    ggsave(file.path(resdir, 'nonsparce_pca.pdf'), bp1, width = 7, height=5, useDingbats=FALSE)
    write.table(pca1[["x"]][,c("PC1", "PC2")], file = file.path(resdir, 'nonsparce_pca.txt'), sep="\t")
    #plot PCA plot of log transformed fluxes for nonsparse nonzero entries
    nflux=read.delim(file.path(resdir,'kcat_n_flux.tsv'), na.strings = 'NaN')
    nflux_nons=nflux[apply(is.na(nflux[,colnames(nflux) %in% cnames])|nflux[,colnames(nflux) %in% cnames]==0, 1, sum)==0,]
    nflux_nons[,colnames(nflux) %in% cnames]=log10(nflux_nons[,colnames(nflux) %in% cnames])
    
    pca2= prcomp(t(nflux_nons[,colnames(nflux) %in% cnames]))
    bp2=  autoplot(pca2, data=data.frame(sample=factor(cnames, levels = cnames)),colour='sample', size=5) + 
      theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') +
      scale_color_discrete(labels =clabs)
    ggsave(file.path(resdir, 'nonsparseflux_pca.pdf'), bp2, width=7, height=5, useDingbats=FALSE)
    write.table(pca1[["x"]][,c("PC1", "PC2")], file = file.path(resdir, 'nonsparseflux_pca.txt'), sep="\t")
    
    
    
    
    ####### old plotting function that requires sparse data to be log transformed and thus only works with a kapp_ndat that is added a pseudocount
    
    ########
  }
  
  plot_kapp_nsparse(subsys_kapp_ndat)
  
  plot_kapp_subsyssparse = function(subsys_kapp_ndat) {
    #function to plot log values for all enzymes with kapp (including sparse) by subsystem. 
    #INPUT: 
    # - df kapp_ndat: a dataframe of kapp values as produced by NIDLE with unscaled raw kapp values
    #Pseudocount is added
    warning("ATTENTION for the boxplots aof all reactions currently the minimum nonzero value times 10^-4 is added as pseudocount before log 10 transformation")
    kappmat=subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]
    #add pseodocount
    psco=min(kappmat[kappmat!=0 & !(is.na(kappmat))])*0.0001
    subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]+psco
    if (any(subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]<=0, na.rm = TRUE)){stop("Values<=0 in kcat data detected. Aborting log transformation...")}
    #log transfrom
    subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=log10(subsys_kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")])
    #only select data from subsystems with min 3 entries for NIDLE kappsa
    occur=table(subsys_kapp_ndat$sparseSubSystem)
    plotdat1=subsys_kapp_ndat[subsys_kapp_ndat$sparseSubSystem %in% names(occur)[occur>2],c('sparseSubSystem', 'kappmax', 'euk_Kcats', 'prok_Kcats')]
    #only select data from subsystems with min 3 entries from NIDLE eukaryots and procaryots
    occur=table(subsys_kapp_ndat$sparseSubSystem[(!is.na(subsys_kapp_ndat$euk_Kcats))&(!is.na(subsys_kapp_ndat$prok_Kcats))])
    plotdat2=subsys_kapp_ndat[subsys_kapp_ndat$sparseSubSystem %in% names(occur)[occur>2],c('sparseSubSystem', 'kappmax', 'euk_Kcats', 'prok_Kcats')]
    plot_bysubsys= function(plotdat) {
      #plot a boxplot of kcat values by subsystem
      plotdat=expandcol(plotdat, 2:4)#
      colnames(plotdat)=c("Process", "Value", "Source")
      #remove values below -10 as outliers
      plotdat=plotdat[plotdat$Value>-10,]
      plotdat$Source[plotdat$Source %in% "kappmax"]="NIDLE"
      orderssubsys=sort(by(plotdat$Value, as.factor(plotdat$Process), median, na.rm=TRUE))
      return(ggplot(plotdat, aes(x=Process, y=Value, color=Source)) +geom_boxplot() +  scale_x_discrete(limits=names(orderssubsys)) +scale_color_manual(values=cbPalette) + #geom_jitter(shape=16) +
               theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust=1), axis.title.x=element_blank(), text = element_text(size=20))+   labs(y= 'log(Kcat|Kapp)'))
    }
    ggsave(file.path(resdir, 'kcat_kapp_subsys_all.pdf'), plot_bysubsys(plotdat1), height=7, width=14)
    ggsave(file.path(resdir, 'kcat_kapp_subsys.pdf'), plot_bysubsys(plotdat2), height=4, width=13)
  }
  plot_kapp_subsyssparse(subsys_kapp_ndat)
  
  ### heatmap plotting functions called by sparse_hm and nsparse_hm
  abs_hm=function(logkapp_ndat, fname,shrinkfac) {
    #plot the mean centered absolute values shrinking the dynamic range of the color code to 1/shrinkfac
    # - char fname: character vector giving the filename, plot is saved under <fname>_abs.pdf
    #               source dat is saved under <fname>_abs.tsv
    # - num shrinkfac:  an integer giving the factor by which the dynamic area of the color range is shrinked 
    #                   in the heatmap with unscaled mean centere kapp values
    
    library(pheatmap)
    #since we are going to mean center the data remove rows with only one entry
    hm_kapp_ndat=logkapp_ndat[apply(!(is.na(logkapp_ndat[,cnames])), 1, sum)>1,c("Rxns", cnames)]
    #mean center the cnames
    hm_kapp_ndat[,cnames]=hm_kapp_ndat[,cnames]-apply(hm_kapp_ndat[, cnames], 1, mean, na.rm=TRUE)
    #since we include NA values we have incomparable rows in the matrix we want to cluster that leads
    #to errors when pheatmap calls hclust with the output from dist that includes NA
    #therefore manually set NA in the dist output to Inf
    d=dist(hm_kapp_ndat[,cnames])
    d[is.na(d)]=10^50
    #calculate focussed range brake (1/3) of actual values for color code
    mybreaks=(1:101)*(max(hm_kapp_ndat[,cnames], na.rm=TRUE)-min(hm_kapp_ndat[,cnames], na.rm=TRUE))/(101*shrinkfac)+min(hm_kapp_ndat[,cnames], na.rm=TRUE)/shrinkfac
    pdf(file.path(resdir, paste(fname,"_abs.pdf", sep="", collapse="")))
    p=pheatmap(hm_kapp_ndat[,cnames], breaks = mybreaks, clustering_distance_rows = d, clustering_method = "single",
               show_rownames = FALSE, treeheight_row = 0, labels_col = clabs)
    print(p)
    dev.off()
    write.table(hm_kapp_ndat[p$tree_row$order,cnames], file =file.path(resdir, paste(fname,"_abs.tsv", sep="", collapse="")), row.names = FALSE, sep="\t" )
  }
  
  rel_hm=function(logkapp_ndat, fname) {
    #plot the z scaled kapp values clustered by 
    #INPUT:
    # - df logkapp_ndat: data.frame containing log transformedkapp values as produced by NIDLE
    # - char fname: character vector giving the filename, plot is saved under <fname>_zs.pdf
    #               source dat is saved under <fname>_zs.tsv
    library(pheatmap)
    
    #remove all entries with less than 2 values
    hm_kapp_ndat=logkapp_ndat[apply(!(is.na(logkapp_ndat[,cnames])), 1, sum)>2,c("Rxns", cnames)]
    #calculate correlation distances\
    d=amap::Dist(hm_kapp_ndat[,cnames], method ="correlation")
    d_cols=amap::Dist(t(hm_kapp_ndat[,cnames]), method="correlation")
    if (any(is.na(d))|| any(is.na(d_cols))) {
      #set distance of incomparables to more than maximum correlation distance
      warning("incomparables in data detected, setting correlation distance for these pairs to 2.1")
      d[is.na(d)]=2.1
      d_cols[is.na(d_cols)]=2.1
    }
    
    pdf(file.path(resdir, paste(fname,"_zs.pdf", sep="", collapse="")))
    p=pheatmap(hm_kapp_ndat[,cnames], clustering_distance_rows = d, clustering_method = "complete", clustering_distance_cols = d_cols,
               show_rownames = FALSE, scale = "row", treeheight_row = 0, labels_col = clabs)
    print(p)
    dev.off()
    outtab=t(scale(t(hm_kapp_ndat[,cnames])))
    write.table(outtab[p$tree_row$order], file =file.path(resdir, paste(fname,"_zs.tsv", sep="", collapse="")), row.names = FALSE, sep="\t" )
  }
  
  spars_hm =function(kapp_ndat, shrinkfac=8) {
    #function to plot a heatmap of kapp values for kapps with sparse entries
    #INPUT: 
    # - df kapp_ndat: a dataframe of kapp values as produced by NIDLE with unscaled raw kapp values
    # - num shrinkfac:  an integer giving the factor by which the dynamic area of the color range is shrinked 
    #                   in the heatmap with unscaled mean centere kapp values
    warning("ATTENTION for the heatmap of all reactions currently the minimum nonzero value times 10^-4 is added as pseudocount before log 10 transformation")
    kappmat=kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]
    #add pseodocount
    psco=min(kappmat[kappmat!=0 & !(is.na(kappmat))])*0.0001
    kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]+psco
    if (any(kapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]<=0, na.rm = TRUE)){stop("Values<=0 in kcat data detected. Aborting log transformation...")}
    #log transfrom
    logkapp_ndat=kapp_ndat
    logkapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=log10(logkapp_ndat[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")])
    
    
    abs_hm(logkapp_ndat, fname="hm_sp",shrinkfac)
    rel_hm(logkapp_ndat, fname="hm_sp")
    
    
  }
  
  spars_hm(kapp_ndat)
  
  nspars_hm = function(kapp_ndat, shrinkfac=8) {
    #function to plot a heatmap of kapp values for kapps with nonsparse entries (!=NA & !=0)
    #INPUT: 
    # - df kapp_ndat: a dataframe of kapp values as produced by NIDLE with unscaled raw kapp values
    # - num shrinkfac:  an integer giving the factor by which the dynamic area of the color range is shrinked 
    #                   in the heatmap with unscaled mean centere kapp values
    
    kapp_ndat_nons=kapp_ndat[apply(is.na(kapp_ndat[, colnames(kapp_ndat) %in% cnames]) | kapp_ndat[,colnames(kapp_ndat) %in% cnames]==0, 1, sum)==0,]
    if (any(kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]<=0, na.rm = TRUE)){stop("Values<=0 in kcat data detected. Aborting log transformation...")}
    #log transfrom
    kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")]=log10(kapp_ndat_nons[,c('viri_Kcats','euk_Kcats', 'prok_Kcats', cnames, "kappmax")])
    abs_hm(kapp_ndat_nons, fname="hm_nonsp", shrinkfac)
    rel_hm(kapp_ndat_nons, fname="hm_nonsp")
  } 
  
  nspars_hm(kapp_ndat)
  return(subsys_kapp_ndat)
}
subsys_kapp_ndat=plot_nidle_kcat(cnames, clabs, TRUE) 


plotexternalcomp = function(cnames, subsys_kapp_ndat) {
  #compare kcat values from GECKO and PFBA to NIDLE
  #INPUT:
  # - char cnames:  character vector with column names of conditions in NIDLE output
  # - df subsys_kapp_ndat:  Data frame of raw NIDLE values with kappmax column (maximum over conditions)
  #                         and cleaned up subsystem info from universal irreversible model
  comp_pfba = function (kapp_ndat, cnames) {
    #compare pFBA and NIDLE results
    #INPUT:
    #INPUT:
    # - df kapp_ndat: a dataframe of kapp values as produced by NIDLE with unscaled raw kapp values
    #- char cnames:  character vector with column names of conditions in NIDLE output
    require(VennDiagram)
    kapp_pdat=read.delim('Results/pFBA//kcat_p.tsv', na.strings = 'NaN')
    if (!all(cnames %in% colnames(kapp_pdat))) {
      stop("Not all cnames are contained as column names in kcat_p.tsv")
    }
    
    
    kapp_pdat$kappmax=apply(kapp_pdat[, colnames(kapp_pdat) %in% cnames],1,max, na.rm=TRUE)
    
    np_comptab=merge(kapp_ndat[,c("Reaction", "JGIgeneID", "UniprotgeneID", "ECNumber", "kappmax")], kapp_pdat[,c("Rxns", "kappmax")], by.x="Reaction", by.y="Rxns", all = TRUE)
    colnames(np_comptab)[colnames(np_comptab) %in% "kappmax.x"]="NIDLE"
    colnames(np_comptab)[colnames(np_comptab) %in% "kappmax.y"]="pFBA"
    pvnrho=cor(x=log10(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),'NIDLE']), y=log10(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),'pFBA']), method = "spearman")
    gvnrhotest=cor.test(x=log10(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),'NIDLE']), y=log10(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),'pFBA']), method = "spearman")
    print(paste("Pvalue for the spearman correlation of NIDLE and pFBA:", signif(gvnrhotest$p.value, 2)))
    pvn=ggplot(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),], aes(x=NIDLE, y=pFBA)) + geom_point(size=4, alpha=0.7) +
      scale_x_log10() + scale_y_log10() + ylab('max(kapp) from pFBA') + xlab('max(kapp) from NIDLE') +
      annotate(geom='text', x=min(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),"NIDLE"])+3, y=max(np_comptab[!(is.na(np_comptab$NIDLE))&!(is.na(np_comptab$pFBA)),"pFBA"])-1, size=8,  label=paste('\u03c1 = ', round(pvnrho, 2))) +
      geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text=element_text(size=25))
    ggsave("Results/NIDLE/kcat_pFBAvn.pdf",pvn, width=7, height=7, device=grDevices::cairo_pdf)
    #plot venn diagram
    venn.diagram(list(NIDLE=kapp_ndat$Reaction, pFBA=kapp_pdat$Reaction), filename="Results/NIDLE/kcat_pFBAvn_venn.png", imagetype = "png")
  }
  comp_pfba(subsys_kapp_ndat, cnames)
  
  comp_GKO = function (kapp_ndat, cnames) {
    #INPUT:
    # - df kapp_ndat: a dataframe of kapp values as produced by NIDLE with unscaled raw kapp values
    #- char cnames:  character vector with column names of conditions in NIDLE output
    #Import GECKO kcats - these is one value for each reaction in the irreversible model which is linked to a GPR rule
    GECKOkcat=read.csv('Program/deps/GECKOcre/models/geCre1355Autotrophic_Rep1/geCre1355Autotrophic_Rep1_kcats.tab')
    
    #match reaction specific kcats
    
    nidle_cols=c( cnames, 'kappmax', 'sparseSubSystem', 'Macrosystem')
    comp_mat=as.data.frame(matrix(nrow=nrow(GECKOkcat), ncol=length(nidle_cols)+2))
    colnames(comp_mat)=c('GECKO', 'GKO_Match' , nidle_cols)
    for (i in 1:nrow(comp_mat)) {
      #take maximum kcat for GECKO 
      if (nchar(GECKOkcat$kcats[i])>0) {#
        max_idx=which.max(as.numeric(unlist(strsplit(GECKOkcat$kcats[i], ';'))))
        if (length(unlist(strsplit(GECKOkcat$kcats[i], ';')))==length(unlist(strsplit(GECKOkcat$MScores[i], ';')))){
          comp_mat[i,1]=as.numeric(unlist(strsplit(GECKOkcat$kcats[i], ';')))[max_idx]
          comp_mat[i,2]=as.numeric(unlist(strsplit(GECKOkcat$MScores[i], ';')))[max_idx]
        }
      }
      rexn=GECKOkcat$rxns[i]
      #check if reversible reaction
      if (grepl('_REV', rexn)) {
        rexn=gsub('_REV', '', rexn)
        nkcat=kapp_ndat[kapp_ndat$Reaction %in% paste(rexn, '_b', sep=""), nidle_cols]
      } else {
        if (rexn %in% kapp_ndat$Reaction) {
          nkcat=kapp_ndat[kapp_ndat$Reaction %in% rexn, nidle_cols]
        } else {
          nkcat=kapp_ndat[kapp_ndat$Reaction %in% paste(rexn, '_f', sep="", collapse = ""), nidle_cols]
        }
      }
      if (nrow(nkcat)>1) {
        stop("one to many mapping in reactionkcat table detected")
      }
      if (nrow(nkcat)>0) {comp_mat[i,3:ncol(comp_mat)]=nkcat}
    }
    
    sc_plot = function(comp_mat, cnames, resdir) {
      #Plot scatter plots of GECKO assigned kcats and NIDLE calculate kapp
      #overall and per GECKO match criteria
      #INPUT:
      #- df comp_mat: comparison matrix between GECKO and NIDLE raw max cat values
      #- char cnames:  character vector with column names of conditions in NIDLE output
      #- char resdir: character giving the directory in which scatter plots are to be saved
      if (!(dir.exists(resdir))) {
        dir.create(resdir)
      }
      comp_mat=as.data.frame(comp_mat)
      gvnrho=cor(x=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax)),'GECKO']), y=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax)),'kappmax']), method = "spearman")
      gvnrhotest=cor.test(x=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax)),'GECKO']), y=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax)),'kappmax']), method = "spearman")
      print(paste("Pvalue for the spearman correlation of NIDLE and GECKO:", signif(gvnrhotest$p.value, 2)))
      gvn=ggplot(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax)),], aes(x=GECKO, y=kappmax, color=as.factor(GKO_Match))) + geom_point(size=4, alpha=0.7) +
        scale_color_manual(name="Match Quality", values=cbPalette[1:length(unique(comp_mat$GKO_Match))], labels=c('Organism + substrate', 'Substrate', 'Organism', 'Any', 'SA organism', 'SA any')) + 
        scale_x_log10() + scale_y_log10() + ylab('max(kapp) from NIDLE') + xlab('kcat from GECKO') +
        annotate(geom='text', x=min(comp_mat$GECKO[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax))])+3, y=max(comp_mat$kappmax[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax))])-1, size=8,  label=paste('\u03c1 = ', round(gvnrho, 2))) +
        geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text=element_text(size=25))
      ggsave(file.path(resdir, "kcat_gvn_all.pdf"),gvn, width=10, height=7, device=grDevices::cairo_pdf)
      orgkcatrho= cor.test(x=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax))&comp_mat$GKO_Match %in% c(1,3), 'GECKO']),
                           y=log10(comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax))&comp_mat$GKO_Match %in% c(1,3), 'kappmax']), method="spearman")
      print(paste('spearman rho for all organism specific kcats:', round(orgkcatrho$estimate,2)))
      print(paste("p-value:", signif(orgkcatrho$p.value, 2)))
      for (m in unique(comp_mat$GKO_Match[!(is.na(comp_mat$GKO_Match))])) {
        mat_labels=c('org+subs', 'subs', 'org', 'any', 'SA org', 'SA_any')
        plot_mat=comp_mat[!(is.na(comp_mat$GECKO))&!(is.na(comp_mat$kappmax))&comp_mat$GKO_Match %in% m,c('GECKO', 'GKO_Match','kappmax')]
        if (nrow(plot_mat>2)){
          rho=cor(x=log10(plot_mat$GECKO), y=log10(plot_mat$kappmax), method="spearman")
          p=ggplot(plot_mat, aes(x=GECKO, y=kappmax)) + geom_point(color=m) + labs(title = mat_labels[m]) +
            scale_x_log10() + scale_y_log10() + ylab('max(kapp) from NIDLE') + xlab('kcat from GECKO') +
            scale_color_discrete(name="kcat match", labels=c('org+subs', 'subs', 'org', 'any', 'SA org', 'SA_any')) +
            annotate(geom='text', x=min(plot_mat[,1:2])+(max(plot_mat[,1:2])-min(plot_mat[,1:2]))/3, y=max(plot_mat[,1:2])-1, size=8,  label=paste('\u03c1 = ', round(rho, 2))) +
            geom_abline(intercept = 0, slope = 1) + theme_bw() + theme(text=element_text(size=25))
          ggsave(file.path(resdir, paste('kcat_gvs_', mat_labels[m], '.pdf',collapse ='')), p, device=grDevices::cairo_pdf)
        }
      }
    }
    sc_plot(comp_mat, cnames, "Results/eccomp/gvn_scat")
    vio_plot = function(comp_mat, cnames, resdir){
      #Plot violin plots of the distribution of kapp/kcat for kapp max, and the conditions
      if (!(dir.exists(resdir))) {
        dir.create(resdir)
      }
      #calculate the ration of kappp divided by kcat for the sparse set
      ##TODO: set 0 values to NA to exclude them from plotting
      rat_comp_mat=comp_mat
      rat_comp_mat[!(is.na(rat_comp_mat)) & rat_comp_mat==0]=NA
      rat_comp_mat=rat_comp_mat[!is.na(comp_mat$GECKO) & apply(!is.na(comp_mat[,cnames]), 1, sum)!=0,]
      rat_comp_mat[,c(cnames, 'kappmax')]=rat_comp_mat[,c(cnames, 'kappmax')]/rat_comp_mat$GECKO
      #if zeros are included add pseudocount
      if (any(rat_comp_mat[,cnames]==0, na.rm=TRUE)) {
        error("0 values in comp_mat detected this is an  internal error :(")
      }
      
      mkplot_cond= function(rat_comp_mat, cnames, fp) {
        #Function to plot the violin plot or kapp/kcat with dots for observations subsdetted by condition
        #INPUT:
        #  - rat_comp_mat: the matrix of raw ratios between calculated kapp values and kcat values having cnames as column names
        #- char cnames:  character vector with column names of conditions in NIDLE output
        # - char fp:  charcter vector containing the file path to which a pdf of the violin plot is to be saved
        
        source("Program/deps/utilities/expandcolv3.r", local =TRUE)
        #transform to long format
        plotdat=expandcol(rat_comp_mat[,c(cnames, "sparseSubSystem",  "Macrosystem")], 1:length(cnames))
        #remove any NA values
        plotdat=plotdat[!is.na(plotdat$x1),]
        colnames(plotdat)[colnames(plotdat) %in% "x1"]="Ratio"
        colnames(plotdat)[colnames(plotdat) %in% "x2"]="Condition"
        plotdat$Ratio=log10(plotdat$Ratio)
        #function to produce the violin plot
        
        p= ggplot(data=plotdat, aes(x=Condition, y=Ratio)) + geom_violin(alpha=0.3, fill="grey30", color=NA) + geom_dotplot(binaxis = "y", stackdir ="center", dotsize = 0.1, stackratio = 0.8) + 
          theme_light() + ylab(expression(log[10]~kapp[NIDLE]/kcat[GECKO])) + scale_x_discrete(labels=clabs) +
          theme(text = element_text(size=25), axis.title.x = element_blank(), axis.text.x=element_text(angle=60, hjust=1)) 
        ggsave(fp, p, width = length(unique(plotdat$Condition))*1.5, height = 8)
        #save sourcedata file
        write.table(plotdat[, c("Ratio", "Condition")], file = gsub("pdf$", "tsv", fp), row.names=FALSE, sep="\t")
      }
      # plot for different sets
      mkplot(rat_comp_mat, cnames, "Results/eccomp/kappkcatratios_diffset.pdf")
      # plot for samesets
      mkplot(rat_comp_mat[!apply(is.na(rat_comp_mat[,cnames]), 1, any),], cnames, "Results/eccomp/kappkcatratios_sameset.pdf")
      
      
      mkplot_macsys= function(rat_comp_mat, cnames, fp) {
        #Function to plot the violin plot of max(kapp)/kcat with dots for observations subsdetted by macrosystem
        #INPUT:
        #  - rat_comp_mat: the matrix of raw ratios between calculated kapp values and kcat values having cnames as column names
        #- char cnames:  character vector with column names of conditions in NIDLE output
        # - char fp:  charcter vector containing the file path to which a pdf of the violin plot is to be saved
        
        #source("Program/deps/utilities/expandcolv3.r", local =TRUE)
        #filter for kappmax
        plotdat=rat_comp_mat[, c("sparseSubSystem", "Macrosystem", "kappmax")]
        #check for any NA values
        if (any(is.na(plotdat$kappmax))) {
          stop("Encounter NA values in kappmax column of comp_mat. Check the indexing of comp_mat")
        }
        plotdat$kappmax=log10(plotdat$kappmax)
        #function to produce the violin plot
        
        p= ggplot(data=plotdat, aes(x=Macrosystem, y=kappmax, fill=Macrosystem)) + geom_violin(alpha=0.3, color=NA) + geom_dotplot(binaxis = "y", stackdir ="center", dotsize = 0.1, stackratio = 0.8) + 
          theme_light() + ylab(expression(log[10]~kapp[NIDLE]/kcat[GECKO])) + scale_fill_manual(values=cbPalette) +
          theme(legend.position = "none", text = element_text(size=25), axis.title.x = element_blank(), axis.text.x=element_text(angle=60, hjust=1)) 
        ggsave(fp, p, width = length(unique(plotdat$Macrosystem))*1.5, height = 8)
        #save sourcedata file
        write.table(plotdat, file = gsub("pdf$", "tsv", fp), row.names=FALSE, sep="\t")
      }
      mkplot_macsys(rat_comp_mat, cnames, "Results/eccomp/kappcatratiospermacsys.pdf")
    }
  }
  comp_GKO(subsys_kapp_ndat, cnames)
}
plotexternalcomp(cnames, subsys_kapp_ndat) 

plotvdis = function(cnames, resdir) {
  #Function to plot the comulative distribution of flux values for the different conditions
  #- char cnames:  character vector with column names of conditions in NIDLE output
  # - char resdir: character vector giving the path to save figures to
  source('Program/deps/utilities/expandcolv3.r', local=TRUE)
  #Import fluxdata
  fluxdat=read.delim("Results/NIDLE/kcat_n_flux.tsv")

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
      scale_color_discrete(labels=clabs, name="Conditions") +theme(text=element_text(size=25), axis.text.x = element_blank()) + 
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
