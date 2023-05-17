plot_nidle_kapp= function(cnames, clabs, same_set) {
  #Code to plot the internal comparison of NIDLE kapp
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
    write.table(table(max_idx), file = file.path(resdir, 'maxpercond.tsv'), sep = "\t", row.names=FALSE)
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
    glob_kcats=merge(glob_kcats, kapp_ndat[, c("Rxns",cnames, "kappmax")],by.x="Reaction", by.y="Rxns", all.x=TRUE)
    kapp_ndat=merge(glob_kcats[,1:2], kapp_ndat, by.x="Reaction", by.y="Rxns", all.y = TRUE)
    #plot a stacked barplot with the kcat information available from viridiplantae
    subsyslevels=unique(unlist(strsplit(glob_kcats$SubSystem, split = "\t")))
    plotdat6=lapply(subsyslevels, function(x) {
      no_kcat=sum(!is.na(glob_kcats$viri_Kcats[grepl(x, glob_kcats$SubSystem)])) #Number of reactions with kcats from viridplantae in sabioRK and BRENDA
      no_kapp=sum(is.na(glob_kcats$viri_Kcats[grepl(x, glob_kcats$SubSystem)]) & #Number of reactions WITHOUT kcats available from sabio RK and BRENDA BUT with kappmax
                    !is.na(glob_kcats$kappmax[grepl(x, glob_kcats$SubSystem)])) 
      no_kappfull=sum(!apply(is.na(glob_kcats[grepl(x, glob_kcats$SubSystem), cnames]),1, any))#Total number of reactions that have kapp in all conditiosn
      if (any(c(no_kcat, no_kapp, no_kappfull)!=0)) {
        return(data.frame(subsys=rep(x,3), type=c("Viridiplantae kcat", "NIDLE kapp", "core NIDLE kapp"),
                          value=c(no_kcat, no_kapp, no_kappfull)))
      } else {
        return(NULL)
      }})
    plotdat6=do.call(rbind, plotdat6)
    plotdat7=plotdat6[plotdat6$type %in% c("Viridiplantae kcat", "NIDLE kapp"),]
    plotdat7$subsys=factor(plotdat7$subsys, levels=plotdat7$subsys[which(plotdat7$type %in% "NIDLE kapp")[order(plotdat7$value[plotdat7$type %in% "NIDLE kapp"])]])
    ggsave(file.path(resdir,'determinedkmax.pdf'), ggplot(plotdat7, aes(x=subsys, y=value, fill=type)) + geom_bar(position = "stack", stat="identity") +
             scale_fill_manual(values=cbPalette[c(2,1)], name="Kinetic data source") + ylab('No of Reactions') + xlab("Subsystem") + theme_classic() + 
             theme(axis.text.x = element_blank(), text=element_text(size=20)), width = 12, height = 4)
    ggsave(file.path(resdir,'determinedkmaxlabels.pdf'), ggplot(plotdat7, aes(x=subsys, y=value, fill=type)) + geom_bar(position = "stack", stat="identity") +
             scale_fill_manual(values=cbPalette[c(2,1)], name="Kinetic data source") + ylab('No of Reactions') + xlab("Subsystem") + theme_classic() + 
             theme(axis.text.x = element_text(angle=60, hjust=1), text=element_text(size=20)), width = 14, height = 10)
    write.table(plotdat7, file = file.path(resdir, 'determinedkmax.tsv'),sep = "\t", row.names=FALSE)
    #Plot a figure only with core values
    plotdat8=plotdat6[plotdat6$type %in% "core NIDLE kapp",]
    plotdat8=plotdat8[plotdat8$value>0,]
    plotdat8$subsys=factor(plotdat8$subsys, levels=as.character(plotdat8$subsys)[order(plotdat8$value)])
    ggsave(file.path(resdir,'determinedcorekmaxlabels.pdf'), ggplot(plotdat8, aes(x=subsys, y=value)) + geom_bar(position = "stack", stat="identity") +
             ylab('No of Reactions') + xlab("Subsystem") + theme_classic() + 
             theme(axis.text.x = element_text(angle=60, hjust=1), text=element_text(size=20)), width = 14, height = 10)
    write.table(plotdat8, file = file.path(resdir, 'determinedcorekmaxlabels.tsv'),sep = "\t", row.names=FALSE)
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
    # gnerate color scale
    colorval=scales::hue_pal()(length(cnames)+1)
    names(colorval)=c(cnames, 'kappmax')
    pca1 <- prcomp(t(kapp_ndat_nons[,colnames(kapp_ndat) %in% c(cnames, "kappmax")]))
    bp1 <- autoplot(pca1, data=data.frame(sample=factor(c(cnames, "kappmax"), levels =c(cnames, "kappmax"))),colour='sample', size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') +
      scale_color_manual(values=colorval,labels = c(clabs, 'maxKapp'))
    ggsave(file.path(resdir, 'nonsparce_pca.pdf'), bp1, width = 7, height=5, useDingbats=FALSE)
    write.table(pca1[["x"]][,c("PC1", "PC2")], file = file.path(resdir, 'nonsparce_pca.txt'), sep="\t")
    #plot PCA plot of log transformed fluxes for nonsparse nonzero entries
    nflux=read.delim(file.path(resdir,'kcat_n_flux.tsv'), na.strings = 'NaN')
    nflux_nons=nflux[apply(is.na(nflux[,colnames(nflux) %in% cnames])|nflux[,colnames(nflux) %in% cnames]==0, 1, sum)==0,]
    nflux_nons[,colnames(nflux) %in% cnames]=log10(nflux_nons[,colnames(nflux) %in% cnames])
    
    pca2= prcomp(t(nflux_nons[,colnames(nflux) %in% cnames]))
    bp2=  autoplot(pca2, data=data.frame(sample=factor(cnames, levels = cnames)),colour='sample', size=5) + 
      theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') +
      scale_color_manual(values=colorval[cnames],labels =clabs)
    ggsave(file.path(resdir, 'nonsparseflux_pca.pdf'), bp2, width=7, height=5, useDingbats=FALSE)
    write.table(pca1[["x"]][,c("PC1", "PC2")], file = file.path(resdir, 'nonsparseflux_pca.txt'), sep="\t")
    
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
    #calculate sd as annotation column
    sdanno=data.frame(sd=apply(hm_kapp_ndat, 1, sd, na.rm=TRUE))
    print(paste(sum(sdanno$sd<2)/nrow(sdanno), "reactions have a sd below 2. With the maximum being", max(sdanno$sd)))
    #generate row annotation color scheme
    #row_anno_color=list(sd=colorRampPalette(rev(brewer.pal(n = 7, name ="Greens")))(101))
    print(paste("The range of values for the absolute heatmap", fname, "is", min(hm_kapp_ndat[,cnames], na.rm=TRUE), 
                "to", max(hm_kapp_ndat[,cnames], na.rm=TRUE)))
    #calculate focussed range of breaks actual values for color code
    mybreaks=(1:101)*(max(hm_kapp_ndat[,cnames], na.rm=TRUE)-min(hm_kapp_ndat[,cnames], na.rm=TRUE))/(101*shrinkfac)+min(hm_kapp_ndat[,cnames], na.rm=TRUE)/shrinkfac
    pdf(file.path(resdir, paste(fname,"_abs.pdf", sep="", collapse="")))
    p=pheatmap(hm_kapp_ndat[,cnames], breaks = mybreaks, annotation_row = sdanno, 
               clustering_distance_rows = d, clustering_method = "single",
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