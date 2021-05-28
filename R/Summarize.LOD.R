#' Summarize the LOD scores from parentage inference with FRANz
#'
#' This function summarizes the LOD scores from parentage inference of PseudoBabies simulations with the software FRANz.
#' @param AnalysisName Name of the Analysis
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param SavePlot Do you want to create plots? Defaults to FALSE
#' @param PlotName What do you want to name your plots? Defaults to LODsum
#' @param Plots What plots do you want to create? Options are 'LCA','PLCA','PLCA','PLICA'.
#' @param ErrorLab Label for Error in plots
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @keywords FRANz parentage LOD
#' @export
#' @examples
#' summarize.LOD()

Summarize.LOD<-function(AnalysisName,nSim,Miss_data=FALSE,MissingVec=c(0,2.5),Geno_Error=FALSE,ErrorVec=c(0,2.0),SavePlot=FALSE,PlotName='LODsum',Plots=c('LCA','PLCA','PLCA','PLICA'),ErrorLab,Markerfile,...){
  wd<-getwd()
  markers<-read.csv(Markerfile)
  Panel<-paste(markers[,2],' uSat & ',markers[,3],' SNPs',sep='')
  
  if(Miss_data==F){MissingVec<-0}
  if(Geno_Error==F){ErrorVec<-0
  }else{
    ErrorVec<-seq(1,ErrorVals,1)
  }  
  
  #LOD for each assignment type
  c<-1
  Res<-list()
  
  for (z in 1:nrow(markers)){
    mark<-markers[z,-1]
    for(e in (1:length(ErrorVec))){
      for (qi in (1:length(MissingVec))){
        LCA<-tryCatch(read.table(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/MissDat_',MissingVec[qi],'/','LODofCorrAss.txt',sep=""),sep=" "), error=function(e) NA)
        LICA<-tryCatch(read.table(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/MissDat_',MissingVec[qi],'/','LODofInCorrAss.txt',sep=""),sep=" "), error=function(e) NA)
        PLCA<-tryCatch(read.table(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/MissDat_',MissingVec[qi],'/','PairLODofCorrAss.txt',sep=""),sep=" "), error=function(e) NA)
        PLICA<-tryCatch(read.table(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/MissDat_',MissingVec[qi],'/','PairLODofInCorrAss.txt',sep=""),sep=" "), error=function(e) NA)
        
        LODinfo<-rbind(LCA,LICA,PLCA,PLICA)
        Res[[c]]<-cbind(LODinfo,
                        (c(rep('LCA',nrow(LCA)),rep('LICA',ifelse(any(is.na(LICA)),1,nrow(LICA))),
                           rep('PLCA',nrow(PLCA)),rep('PLICA',ifelse(any(is.na(PLICA)),1,nrow(PLICA))))),
                        rep(MissingVec[qi],length(LODinfo)),
                        rep(ErrorVec[e],length(LODinfo)),
                        rep(paste(markers[z,2],' uSats & ',markers[z,3],' SNPs',sep=""),length(LODinfo))) 
        
        
        c<-c+1
      }
    }
  }
  
  Res<-do.call(rbind,Res)
  colnames(Res)<-c('Val','Type','MissDat','Error','Panel')
  Res<-Res[-which(is.na(Res[,1])),]
  
  mapError<-matrix(data=c(unique(Res$Error),ErrorLab),ncol=2,nrow=length(unique(Res$Error)),
                   byrow=F)
  Res$ErrorMap<-mapvalues(Res$Error,mapError[,1],mapError[,2])
  
  addline_format <- function(x,...){
    gsub('\\s\\s','\n',x)
  }
  
  
  mf_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="Type") { 
      value[value=="LCA"] <- "LOD correct"
      value[value=="LICA"]   <- "LOD incorrect"
      value[value=="PLCA"] <- "Pair LOD  correct parent"
      value[value=="PLICA"]   <- "Pair LOD  incorrect parent"
    }
    return(value)
  }
  
if('LCA'%in% Plots){  
p <- ggplot(Res[Res$Type=='LCA',], aes(factor(MissDat), Val)) + 
    geom_boxplot(aes(fill = factor(Error))) +
    xlab('Percent Missing Data')+
    ylab('LOD') +
    scale_fill_discrete(name="Genotyping\nError",
                        breaks=factor(ErrorVec),
                        labels=ErrorLab)+
    facet_wrap( ~ Panel,ncol=4)+
    theme(legend.justification=c(1,0.15), legend.position=c(1,0))+
    labs(title="LOD of individuals assigned to correct parents")
  
  if(SavePlot==T){                            
    ggsave(plot=p,filename=paste(PlotName,'_LCA.pdf',sep=''),dpi=300,width=6,height=6,units='in')
  }
  p
}

if('LICA'%in% Plots){  
  q <- ggplot(Res[Res$Type=='LICA',], aes(factor(MissDat), Val)) + 
    geom_boxplot(aes(fill = factor(Error))) +
    xlab('Percent Missing Data')+
    ylab('LOD') +
    scale_fill_discrete(name="Genotyping\nError",
                        breaks=factor(ErrorVec),
                        labels=ErrorLab)+
    facet_wrap( ~ Panel,ncol=4)+
    theme(legend.justification=c(1,0.15), legend.position=c(1,0))+
    labs(title="LOD of individuals assigned to incorrect parents")
  
  if(SavePlot==T){                            
    ggsave(plot=q,filename=paste(PlotName,'_LICA.pdf',sep=''),dpi=300,width=8.5,height=11,units='in')
  }
  plot(p)
}
if('PLCA'%in% Plots){  
  r<- ggplot(Res[Res$Type=='PLCA',], aes(factor(MissDat), Val)) + 
    geom_boxplot(aes(fill = factor(Error))) +
    xlab('Percent Missing Data')+
    ylab('LOD') +
    scale_fill_discrete(name="Genotyping\nError",
                        breaks=factor(ErrorVec),
                        labels=ErrorLab)+
    facet_wrap( ~ Panel,ncol=4)+
    theme(legend.justification=c(1,0.15), legend.position=c(1,0))+
    labs(title="Pair LOD of individuals assigned to correct parents")
  
  if(SavePlot==T){                            
    ggsave(plot=r,filename=paste(PlotName,'_PLCA.pdf',sep=''),dpi=300,width=8.5,height=11,units='in')
  }
  #plot(r)
}

if('PLICA'%in% Plots){  
  s<- ggplot(Res[Res$Type=='PLICA',], aes(factor(MissDat), Val)) + 
    geom_boxplot(aes(fill = factor(Error))) +
    xlab('Percent Missing Data')+
    ylab('LOD') +
    scale_fill_discrete(name="Genotyping\nError",
                        breaks=factor(ErrorVec),
                        labels=ErrorLab)+
    facet_wrap( ~ Panel,ncol=4)+
    theme(legend.justification=c(1,0.15), legend.position=c(1,0))+
    labs(title="Pair LOD of individuals assigned to incorrect parents")
  
  if(SavePlot==T){                            
    ggsave(plot=s,filename=paste(PlotName,'_PLICA.pdf',sep=''),dpi=300,width=8.5,height=11,units='in')
  }
  #plot(s)
}
  LODsummary<-Res
  return(LODsummary)
}#summarize.LOD

