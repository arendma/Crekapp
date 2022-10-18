# Script for the summary statistics of Davids QConCat data 20220124
read_qcc=function() {
  #Reads in  data from QCC experiments and cleans of ambigous gene names and proteins only observed in a single biological replicate
  require(ggplot2)

  expandcol <- function(df, cols) {
    #takes a data frame and expands it by truning columns in key value pairs
    if(length(cols)<dim(df)[2]) {
      oldcols <- df[,-cols]
    }
    resdf <- data.frame()
    for (i in cols) {
      newdf=data.frame(oldcols,x1=df[,i], x2=rep(colnames(df)[i], dim(df)[1]))
      resdf=rbind(resdf,newdf)
    }
    return(resdf)
  }
  recur_union=function(ls) {
    #function that recursively builds the union of all elements in ls
    if(length(ls)==1) {
      return(ls[[1]])
    } else {
      ls[[2]]=union(ls[[1]], ls[[2]])
      ls[[1]]=NULL
      recur_union(ls)
    }
  }
SBP_QCC=read.delim('Data/AbolsuteQuantUVM4STOP1STOP2.txt')
UVM_QCC=SBP_QCC[, c('ProteinId', 'UVM4_1', 'UVM4_2', 'UVM4_3', "Synonyms")]
UVM_QCC=UVM_QCC[apply(!(is.na(UVM_QCC[, 2:4])), 1, any), ]
colnames(UVM_QCC)[colnames(UVM_QCC) %in% "Synonyms"]="Synonym"
stop1_QCC=SBP_QCC[, grep('Protein|Stop1\\_|Synonyms', colnames(SBP_QCC))]
stop1_QCC=stop1_QCC[apply(!(is.na(stop1_QCC[, 2:4])), 1, any), ]
colnames(stop1_QCC)[colnames(stop1_QCC) %in% "Synonyms"]="Synonym"
stop2_QCC=SBP_QCC[, grep('Protein|Stop12\\_|Synonyms', colnames(SBP_QCC))]
stop2_QCC=stop2_QCC[apply(!(is.na(stop2_QCC[, 2:4])), 1, any), ]
colnames(stop2_QCC)[colnames(stop2_QCC) %in% "Synonyms"]="Synonym"
colnames(stop2_QCC)[2:4]=c('Stop2_1', 'Stop2_2', 'Stop2_3')
cont_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantControlTypeAnnotated.txt')
colnames(cont_QCC)=gsub('Replicate', 'control', gsub('_.amol.cell.', '', colnames(cont_QCC), fixed = TRUE))
dark_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantDarkAnnotated.txt')
colnames(dark_QCC)=gsub('Replicate', 'dark', gsub('_.amol.cell.', '', colnames(dark_QCC), fixed = TRUE))
hcell_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantHighCellDensityAnnotated.txt')
colnames(hcell_QCC)=gsub('Replicate', 'highcell', gsub('_.amol.cell.', '', colnames(hcell_QCC), fixed = TRUE))
hsalt_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantHighSaltAnnotated.txt')
colnames(hsalt_QCC)=gsub('Replicate', 'highsalt', gsub('_.amol.cell.', '', colnames(hsalt_QCC), fixed = TRUE))
htemp_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantHighTemperatureAnnotated.txt')
colnames(htemp_QCC)=gsub('Replicate', 'hightemp', gsub('_.amol.cell.', '', colnames(htemp_QCC), fixed = TRUE))
noshak_QCC=read.delim('Data/QconCAT_David20220124/abs_abundance/AbsoluteQuantNoShakingAnnotated.txt')
colnames(noshak_QCC)=gsub('Replicate', 'noshaking', gsub('_.amol.cell.', '', colnames(noshak_QCC), fixed = TRUE))

#Cleaning the transcripts and keeping a statistic
QCClist=list(cont_QCC, dark_QCC, hcell_QCC, hsalt_QCC, htemp_QCC, noshak_QCC, UVM_QCC, stop1_QCC, stop2_QCC)
condNames=c('control', 'dark', 'highcell', 'highsalt', 'hightemp', 'noshaking', 'UVM', 'SBP_stop1', 'SBP_stop2')
rawQCClist=QCClist
#remove proteins that are NA in all entries
QCClist=lapply(QCClist, function(df) {return(df[apply(!is.na(df[2:4]), 1,any),])})
#all reported protein names
cleaning_stat=data.frame(condNames, raw=sapply(QCClist, nrow))

#removing entries without Cre protein Id
QCClist=lapply(QCClist, function(df) {return(df[grepl('Cre', df[, "ProteinId"]), ])})
cleaning_stat$only_Cre=sapply(QCClist, nrow)

#clean up entries

clean_ambig= function(df) {
    #STEP 1 find concatamer entries and only keep these entries for the respective proteins
    concat_entry=grep('CpQconCAT', df$ProteinId)
    rem_idx=numeric()
    for (i in concat_entry) {
      #if there is onlyone other entry and this contains cre
      entry=gsub(' ', '', unlist(strsplit(df$ProteinId[i], split=';')))
      if (length(entry)==2 && any(grepl('Cre', entry))){
        #remove all other entries of this protein from the data frame
        rem_idx=c(rem_idx, setdiff(grep(entry[grepl('Cre', entry)], df$ProteinId, fixed = TRUE), i))
        #only keep the Cre nam3e
        df$ProteinId[i]=entry[grepl('Cre', entry)]
      } else {
        #remove this entry
        rem_idx=c(rem_idx, i)
      }
    }
    if (length(rem_idx)!=0){
    df=df[-rem_idx,]}
    
    #STEP 2 remove measured abundance and protein names iteratively from ambigous entries
    unfinished=TRUE
    while (unfinished) {
      sav_nrow=nrow(df)
    protnames=gsub(' ', '', unlist(strsplit(df$ProteinId, split=';')))
    dups=protnames[duplicated(protnames)]
    for (dup in dups) {
      #function that checks if their is an unambigous entry for a protein and substract the 
      #unambigous amount form all ambigous entries and removes the protein id from those entries
      ##ATTENTION as the for loop progresses entries with only one duplicated Cre ID are created
      #These have to be resolved later
      
      #if unambigous entry exists and these do not contain NA values
      if (sum(df$ProteinId %in% dup)==1 && !(any(is.na(df[df$ProteinId %in% dup, 2:4])))) {

        protn=df[df$ProteinId %in% dup, 2:4]
        #remove protein name from ambigous entries
        ambig_idx=grepl(dup, df$ProteinId, fixed=TRUE) & !(df$ProteinId %in% dup)
        if (sum(ambig_idx)>1) {
          print(dup)
        }
        ##ATTENTION During the course of the for loop previously duplicated entries may become unique
        if (sum(ambig_idx)>0) {
        df$ProteinId[ambig_idx]=sapply(df$ProteinId[ambig_idx], function (x) {
          return(paste(unlist(strsplit(x, split = ';'))[!(unlist(strsplit(x, split = ';')) %in% dup)], collapse = ';'))
        })
        #remove protein amount from ambigous entries
        df[ambig_idx,2:4]= df[ambig_idx,2:4]-protn[rep(1, sum(ambig_idx)),]
        #remove entries which do not contain Chlamy ids anymore or contain negative molecular amount
        if (any(!(grepl('Cre', df$ProteinId[ambig_idx])) | apply(df[ambig_idx, 2:4]<=0,1, any, na.rm=TRUE))) {
          print('removing entries')
          print(df[(which(ambig_idx)[!(grepl('Cre', df$ProteinId[ambig_idx])) | apply(df[ambig_idx, 2:4]<=0,1, any, na.rm=TRUE)]),])
          df=df[-(which(ambig_idx)[!(grepl('Cre', df$ProteinId[ambig_idx])) | apply(df[ambig_idx, 2:4]<=0,1, any, na.rm=TRUE)]),]
        }}
      } else if (sum(df$ProteinId %in% dup)>1) {
        print(paste('Multiple unambigous entries for gene', dup, 'detected'))
      }
    }
    ##for newly introduced duplicates that are unambiguous exclude all except the maximum
    excl=unlist(sapply(df$ProteinId[duplicated(df$ProteinId) & !(grepl(';', df$ProteinId))], function(x) {
        #take mean over replicates
        mu_entry=apply(df[df$ProteinId %in% x, 2:4], 1, mean, na.rm=TRUE)
        #return data frame indes of all other entries
        return(which(df$ProteinId %in% x)[!(mu_entry %in% max(mu_entry))])
      }))
    if(length(excl)>0) {df=df[-excl,]}
    #repeat until no entries are removed anymore
    if(nrow(df)==sav_nrow) {unfinished=FALSE}
    }
    ##STEP3 assign the total amount of transcription isoforms to the first one
    for (i in grep(';', df$ProteinId)) {
      if (length(unique(gsub('\\.t.*$', '', unlist(strsplit(df$ProteinId[i], split=';')))))==1) {
        print('Changed')
        print(df$ProteinId[i])
        print(unlist(strsplit(df$ProteinId[i], split=';'))[1])
        df$ProteinId[i]=unlist(strsplit(df$ProteinId[i], split=';'))[1]
      }
    }
    ##for newly introduced duplicates that are unambiguous exclude all except the maximum
    excl=unlist(sapply(df$ProteinId[duplicated(df$ProteinId) & !(grepl(';', df$ProteinId))], function(x) {
      #take mean over replicates
      mu_entry=apply(df[df$ProteinId %in% x, 2:4], 1, mean, na.rm=TRUE)
      #return data frame indes of all other entries
      return(which(df$ProteinId %in% x)[!(mu_entry %in% max(mu_entry))])
    }))
    if(length(excl)>0) {df=df[-excl,]}
    ##remove all other ambigous entries  
    df=df[-(grep(';', df$ProteinId)),]
    return(df)
}
QCClist=lapply(QCClist, clean_ambig)
cleaning_stat$nodup=sapply(QCClist, nrow)

#remove any entries with more than 1 nan value

QCClist=lapply(QCClist, function(df) {return(df[rowSums(is.na(df[,2:4]))<2, ])})
cleaning_stat$rep=sapply(QCClist, nrow)

#plot the cleaning statistics
c_stat_plot=cleaning_stat
for (i in 2:(ncol(c_stat_plot)-1)) {
  c_stat_plot[,i]=c_stat_plot[,i]-c_stat_plot[,(i+1)]
}
colnames(c_stat_plot)[2:5]=c(colnames(c_stat_plot)[3:5], 'final')
c_stat_plot=expandcol(c_stat_plot, 2:5)
colnames(c_stat_plot)=c('Condition', 'value', 'step')
c_stat_plot$step=factor(c_stat_plot$step, levels=c('only_Cre', 'nodup', 'rep', 'final'))

filtplot=ggplot(c_stat_plot, aes(fill=step, y=value, x=Condition)) + geom_bar(position='stack', stat='identity') +
  theme_bw()+theme(text=element_text(size=20), axis.text.x=element_text(angle=45, hjust = 1))
dir.create('Results/QconCat20220124', recursive = TRUE, showWarnings = FALSE)
ggsave('Results/QconCat20220124/filter_plot.pdf', width=10, height = 6,useDingbats=FALSE)
return(list(QCClist, rawQCClist))
}
qccdat=read_qcc()[[1]]



QCC_sum=function(qccdat) {
  #Function to calculate and plot overview statistics for QCC data
  #Input:
  # - qccdat: data frame with mearused trancript/protein names as rownames
  library(ggplot2)
  library(ggfortify)
  #Import Cre1355 GPR matrix (the set of genes in Cre1355 is a subset of CreMora)
  Cre1355gxn=as.matrix(read.delim('Data/Cre1355/Cre1355_transcripts.txt', row.names = 1))
  if (!file.exists("Data/Cre1355/Organellar_geninf.txt")){
    org_id=colnames(Cre1355gxn)[grepl( "Chre[[:alpha:]]", colnames(Cre1355gxn))]
    query=paste(paste("cre:", org_id,sep=""), sep="", collapse="+")
    url='https://rest.kegg.jp/list/'
    print("Paste this into browser to obtain organellar genen info")
    print(paste(url, query, sep=""))
    stop(paste('Could not find info on Organellar gene names in Data/Cre1355/Organellar_geninf.txt'))
  } else {
    org_inf=read.delim("Data/Cre1355/Organellar_geninf.txt", header = FALSE)
  }
    #remove KEGG organism ID 
    org_inf[,1]=gsub('cre\\:', '', org_inf[,1])
    #only kee first entry of genen synonyms
    org_inf[,2]=sapply(strsplit(org_inf[,2], spli=";"), '[', 1)
    #remove ambigous entries
    org_inf=org_inf[!(org_inf[,2] %in% org_inf[duplicated(org_inf[,2]),2]),]
    adapt_CreCp = function(df){
    #function to convert the organellar genome proteins into the model IDs
      for (i in 1:nrow(org_inf)) {
        idx=grep(org_inf[i, 2], df$Synonym,ignore.case = TRUE)
        #check if grep hit is unique and has chloroplast gene_id
        if (length(idx)==1 && grepl("Cre-", df$ProteinId[idx])) {
         df$ProteinId[idx]=org_inf[i,1] 
        }
      }
      return(df)
    }
    qccdat=lapply(qccdat, adapt_CreCp)

recurmerge=function(QCClist) {
  #funcion to recursively merge the QCC data frames
  if (length(QCClist)==1) {
    return(QCClist[[1]])
  } else {
    QCClist[[2]]=merge(QCClist[[1]], QCClist[[2]], by=1, all=TRUE)
    QCClist[[1]]=NULL
    recurmerge(QCClist)
  }
}
#remove synonyme column
qccdat=lapply(qccdat, function(df) {return(df[, !(colnames(df)%in% "Synonym")])})
#merge together all cleaned data frames THBERE ARE NEW AMBIGOUTIES INTRODUCED HERE but for each sample no ambiguities are garantied (NA for ambigous value)
QCCagg=recurmerge(qccdat)
help_f2=function (x,df){
  dfslice=df[, grep(x, colnames(df), fixed=TRUE)]
  if (any(rowSums(is.na(dfslice))==2)) {
  stop(paste('Cases with only one maesurement in triplicates detected for condition', x))
  }
  return(apply(dfslice, 1, median, na.rm=TRUE))
}
med_QCCagg=sapply(unique(gsub('_[[:digit:]]', '', colnames(QCCagg)[2:ncol(QCCagg)])), help_f2, df=QCCagg)
med_QCCagg=data.frame(ProteinId=QCCagg$ProteinId, med_QCCagg)
#export the median abundance for use with NIDLE 
write.table(med_QCCagg, 'Data/QconCAT_David20220124/abs_abundance/med_abun_all.tsv', row.names=FALSE, sep='\t')
write.table(med_QCCagg[!(is.na(med_QCCagg$control)), c('ProteinId', 'control')],'Data/QconCAT_David20220124/abs_abundance/JGIproteinabundances.tsv',row.names = FALSE, col.names=c('prot_ID', 'Abund_amolpercell'))
sum_stat=data.frame(sample=gsub('_[[:digit:]]', '', colnames(QCCagg)[2:ncol(QCCagg)]) , tot_amol=colSums(QCCagg[,2:ncol(QCCagg)], na.rm = TRUE))

#get proteins included in the model (here a single match in gene names is ok)
cre_idx=rowSums(sapply(colnames(Cre1355gxn), grepl, x=QCCagg$ProteinId, fixed=TRUE))>0
sum_stat=cbind(sum_stat, cre_amol=colSums(QCCagg[cre_idx,2:ncol(QCCagg)], na.rm = TRUE))
cre_QCCagg=QCCagg[cre_idx,]
plot1=ggplot(data=sum_stat, aes(x=sample, y=tot_amol)) + stat_summary(fun.y=mean, geom='bar', alpha=0.4) + geom_point() + labs(y='Total Protein content [amol/cell]') +
  stat_summary(aes(y=cre_amol), fun.y=mean, geom='bar', fill='lightpink3') +  scale_x_discrete(labels = c('Control CC1690', 'Dark CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690', 'SDP OE1 UVM4', 'SDP OE2 UVM4', 'Control UVM4')) +
  geom_point(aes(y=cre_amol), color='red') +  theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=50, hjust = 1), text=element_text(size=20))  
ggsave('Results/QconCat20220124/tot_prot.pdf', plot1, useDingbats=FALSE)
#get number of proteins measured
sum_stat=cbind(sum_stat, tot_no=colSums(!(is.na(QCCagg[,2:ncol(QCCagg)]))))
sum_stat=cbind(sum_stat, cre_no=colSums(!(is.na(QCCagg[cre_idx,2:ncol(QCCagg)]))))
plot2=ggplot(data=sum_stat, aes(x=sample, y=tot_no)) + stat_summary(fun.y=mean, geom='bar', alpha=0.4) + geom_point() + labs(y='Proteins quantified') +
  stat_summary(aes(y=cre_no), fun.y=mean, geom='bar', fill='lightpink3') +  scale_x_discrete(labels = c('Control CC1690', 'Dark CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690', 'SDP OE1 UVM4', 'SDP OE2 UVM4', 'Control UVM4')) +
  geom_point(aes(y=cre_no), color='red') +  theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=50, hjust = 1), text=element_text(size=20))  
ggsave('Results/QconCat20220124/no_prot.pdf', plot2, useDingbats=FALSE)
#limit on fully measured proteins
ns_QCCagg=QCCagg[rowSums(is.na(QCCagg))==0,]
ns_cre_QCCagg=cre_QCCagg[rowSums(is.na(cre_QCCagg))==0, ]
pca1 <- prcomp(t(log2(ns_QCCagg[,2:ncol(QCCagg)])))
bp1 <- autoplot(pca1, data=sum_stat,colour='sample' ,size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') +
  scale_color_discrete(labels = c('Control CC1690', 'Dark CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690', 'SDP OE1 UVM4', 'SDP OE2 UVM4', 'Control UVM4'))
ggsave('Results/QconCat20220124/pca_full.pdf', bp1, useDingbats=FALSE)
pca2 <- prcomp(t(log2(ns_cre_QCCagg[,2:ncol(QCCagg)])))
bp2 <- autoplot(pca2, data=sum_stat,colour='sample' ,size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample') + 
  scale_color_discrete(labels = c('Control CC1690', 'Dark CC1690', 'High Cell CC1690', 'High Salt CC1690', 'High Temp CC1690', 'No Shaking CC1690', 'SDP OE1 UVM4', 'SDP OE2 UVM4', 'Control UVM4'))
ggsave('Results/QconCat20220124/pca_full_cre.pdf', bp2, useDingbats=FALSE)
}
QCC_sum(qccdat)
# #There aree aggrgated transcript names as labels in QCC data 
# QCC_tr=unlist(strsplit(qccdat[,1], ';'))
# 
# #Calculate how many transcripts from the model are quantified
# overlap=intersect(QCC_tr, colnames(Cre1355gxn))
# qcc_cre=qccdat[unlist(sapply(overlap, function(x) {return(grep(x, qccdat[,1], fixed=TRUE))})),2:(ncol(qccdat)-1)]
# sum_stat=cbind(sum_stat, cre_amol=colSums(qcc_cre))
# print(paste(length(overlap), 'of', ncol(Cre1355gxn), '( ', length(overlap)/ncol(Cre1355gxn), '% ) Proteins in the model iCre1355 have been quantified in QCC data', collapse=''))
# 
# #Calculate the number of reaction that can be bonded (not taking into account isoenzymes) by measurements)
# bonded=sum(rowSums(Cre1355gxn[, colnames(Cre1355gxn) %in% overlap])>0)
# print(paste('Corresponding to', bonded, 'of', nrow(Cre1355gxn), '(', bonded/nrow(Cre1355gxn), '% ) Reactions that can be bounded by this measurements if isoenzymes are not taken into account'))
# dircreater('Results/QconCat20201223/')
# plot1=ggplot(data=sum_stat, aes(x=sample, y=tot_amol)) + stat_summary(fun.y=mean, geom='bar', alpha=0.4) + geom_point() + labs(y='Total Protein content [amol/cell]') +
#   stat_summary(aes(y=cre_amol), fun.y=mean, geom='bar', fill='lightpink3') + geom_point(aes(y=cre_amol), color='red') +  theme_bw(base_size=14)
# ggsave('Results/QconCat20201223/tot_prot_full.pdf',plot1, height = 7, width = 6)
# 
# #Quick and dirty analysis usually you want to normalize your data before pca and annova
# pca1 <- prcomp(t(qccdat[,2:(ncol(qccdat)-1)]))
# bp1 <- autoplot(pca1, data=sum_stat,colour='sample' ,size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample')
# ggsave('Results/QconCat20201223/pca_full.pdf', bp1)
# pca2 <- prcomp(t(qcc_cre))
# bp2 <- autoplot(pca2, data=sum_stat,colour='sample' ,size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample')
# ggsave('Results/QconCat20220124/pca_full_cre.pdf', bp2)
# 
# anova_rep = function (normlogdat, sample) {
#   #reproduce the annova data in the fill
#   
#   mod=apply(normlogdat, 1, function(x) aov(x~sample))
#   pval=sapply(mod, function (x) {return(summary(x)[[1]][1,5])})
#   #values correspond to values in the excel sheet
#   padj=p.adjust(pval, "BH")
#   return(padj)
# }
# 
# padj=anova_rep(qccdat[,2:(ncol(qccdat)-1)], sum_stat$sample)
# DEnormdat=scale(t(qccdat[padj<0.1, 2:(ncol(qccdat)-1)]))
# pca3=prcomp(DEnormdat)
# bp3 <- autoplot(pca3, data=sum_stat,colour='sample' ,size=5) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Sample')
# ggsave('Results/QconCat20220124/pca_full_denorm.pdf', bp1)
# }
# QCC_sum(full_QCC)
# sparse_QCC=read.delim('Data/QconCAT_David20220124/Proteins_Sparse_AmolPerCell.txt')
# #only keep entries with at least one entr in each condition
# conditions=gsub('_[[:digit:]]', '', colnames(sparse_QCC)[2:(ncol(sparse_QCC)-1)])
# cond_includ= sapply(unique(conditions), function(x){
#   idx=which(conditions %in% x)+1
#   res=apply(sparse_QCC[,idx],1, function(x) {return(!(all(is.na(x))))})
#   return(res)
#   })
# filt_s_QCC=sparse_QCC[apply(cond_includ, 1, all),]
# 
# small_summy = function(qccdat) {
#   #tweak the upper function to be usable for inputs with NA values
# Cre1355gxn=as.matrix(read.delim('Data/Data_S2/Cre1355_transcripts.txt', row.names = 1))
# #There aree aggrgated transcript names as labels in QCC data 
# QCC_tr=unlist(strsplit(qccdat[,1], ';'))
# #Calculate how many transcripts from the model are quantified
# overlap=intersect(QCC_tr, colnames(Cre1355gxn))
# qcc_cre=qccdat[unlist(sapply(overlap, function(x) {return(grep(x, qccdat[,1], fixed=TRUE))})),2:(ncol(qccdat)-1)]
# print(paste(length(overlap), 'of', ncol(Cre1355gxn), '( ', length(overlap)/ncol(Cre1355gxn), '% ) Proteins in the model iCre1355 have been quantified in QCC data', collapse=''))
# #Calculate the number of reaction that can be bonded (not taking into account isoenzymes) by measurements)
# bonded=sum(rowSums(Cre1355gxn[, colnames(Cre1355gxn) %in% overlap])>0)
# print(paste('Corresponding to', bonded, 'of', nrow(Cre1355gxn), '(', bonded/nrow(Cre1355gxn), '% ) Reactions that can be bounded by this measurements if isoenzymes are not taken into account'))
# }
# small_summy(filt_s_QCC)
