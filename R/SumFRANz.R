#' Summarize the results from FRANz analysis for parentage on muliple PseudoBabies simulations
#'
#' This function summarizes the results from multiple runs of Colony2 for parentage inference from PseudoBabies simulations.
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Cutoff cutoff value for probability value accepting parentage.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @keywords FRANz parentage
#' @export
#' @examples
#' SumFRANz()

#### Summarize FRANz ####
SumFRANz<-function(nSim,Miss_data=FALSE,MissingVec=c(0,2.5),Geno_Error=FALSE,ErrorVals,Cutoff,Markerfile){
  library(ggplot2)
  wd<-getwd()
  
  if(file.exists(Markerfile)){
    markers<-read.csv(Markerfile)
  }else{
    stop("Marker file not found.")
  }
  wd<-getwd()
  
  Panel<-paste(markers[,2],' uSats & ',markers[,3],' SNPs',sep='')
  
  if(Miss_data==F){MissingVec<-0}
  if(Geno_Error==F){ErrorVec<-0
  ErrorVals<-length(ErrorVec)
  }else{
    ErrorVec<-seq(1,ErrorVals,1)
  }
  
  c<-1
  
  Res<-matrix(data="NA",ncol=(nSim)*6+3,nrow=ErrorVals*length(MissingVec)*nrow(markers))
  Res[,1]<-rep(Panel,each=ErrorVals*length(MissingVec))
  Res[,2]<-rep(ErrorVec,each=length(MissingVec))
  Res[,3]<-if(length(MissingVec)>1) rep(MissingVec,nrow(markers)) else rep(MissingVec)
  colnames(Res)<-c('Panel','Error','MissDat',paste('Sim',seq(1,nSim,1),sep=''),
                   as.vector(sapply(1:2,function(x) paste(paste("Sim",seq(1,nSim,1),sep=""),paste("Type",c("I","II"),sep="")[x],sep=""))),
                   as.vector(sapply(1:3,function(x) paste(paste("Sim",seq(1,nSim,1),sep=""),c("LOD.Cor.Both","LOD.InCor.One","LOD.InCor.Both")[x],sep=""))))
  
 
  #Set up some lists for all the loops....
  MarkerLOD.CorrBoth<-rep(list(list()),nrow(markers))
  MarkerLOD.InCorrBoth<-rep(list(list()),nrow(markers))
  MarkerLOD.InCorrOne<-rep(list(list()),nrow(markers))
  ErrorLOD.CorrBoth<-rep(list(list()),length(ErrorVec))
  ErrorLOD.InCorrBoth<-rep(list(list()),length(ErrorVec))
  ErrorLOD.InCorrOne<-rep(list(list()),length(ErrorVec))
  MissDataLOD.CorrBoth<-rep(list(list()),length(MissingVec))
  MissDataLOD.InCorrBoth<-rep(list(list()),length(MissingVec))
  MissDataLOD.InCorrOne<-rep(list(list()),length(MissingVec))
  SimLOD.CorrBoth<-rep(list(list()),nSim)
  SimLOD.InCorrBoth<-rep(list(list()),nSim)
  SimLOD.InCorrOne<-rep(list(list()),nSim)
  
  for (z in 1:nrow(markers)){
    mark<-markers[z,-1]
    for(e in (1:length(ErrorVec))){
      for(MS in (1:length(MissingVec)) ){
        Pmd<-MissingVec[MS]
        for (s in 1:nSim){
          
          #system(paste("sed 's/,<//g;s/Posterior Parent 2,/Posterior Parent 2/' ",paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/SimIn_FRANz/SimNum",s,"/parentage.csv > ",sep=""),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/SimIn_FRANz/SimNum",s,"/parentage2.csv",sep=""),sep="")) #clean up parentage.csv so it imports nicely
          Franz.results<- readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/SimNum",s,"/parentage.csv",sep=""))
          Franz.results<- gsub(pattern=",<",replacement="",x=Franz.results)
          Franz.results<- gsub(pattern=",$",replacement="",x=Franz.results)
          Franz.results<-matrix(data=unlist(stringr::str_split(Franz.results,",")),nrow=length(Franz.results),ncol=length(unlist(stringr::str_split(Franz.results[1],","))),byrow=T)
          colnames(Franz.results)<-Franz.results[1,]
          Franz.results<-Franz.results[-1,]
          
          Answer<-read.csv( file = paste(wd,'/SimParents/','Truth_Sim',s,".csv",sep=""),stringsAsFactors=F)
          Answernames<-Answer[,1]
          Answer<-stringr::str_split(Answer[,2],'&')
          names(Answer)<-Answernames
          
          Num.Correct<-vector()
          LOD.Cor.Both<-vector()
          LOD.InCor.Both<-vector()
          LOD.InCor.One<-vector()
          PairLOD.Cor<-vector()
          PairLOD.InCor<-vector()
          TypeI<-vector()
          TypeII<-vector()
          #These are just found placement counting
          CountCorr.Both<-0
          CountInCorr.Both<-0
          CountInCorr.One<-0
          PairCountCor<-1
          PairCountInCor<-1
          Files2Write<-c('PropCorr','LOD.Cor.Both','LOD.InCor.Both','LOD.InCor.Both')#,'PairLOD.Cor','PairLOD.InCor')
          
          for (r in 1:nrow(Franz.results)){
            Offspring_F<-Franz.results[r,1]
            Parents_F<-Franz.results[r,c(3,5)]
            TRUTH.S<-eval(parse(text=paste("Answer$'",Offspring_F,"'",sep='')))
           
            if(Franz.results[r,8] >= Cutoff){ #what about ties to this one...
              Num.Correct[r]<-sum(Parents_F%in%TRUTH.S)
              TypeI[r]<-sum(!sum(Parents_F%in%TRUTH.S))
              TypeII[r]<-sum(Parents_F%in%c(""))
            }else{
              Num.Correct[r]<-0
              TypeI[r]<-sum(!sum(Parents_F%in%TRUTH.S))
              TypeII[r]<-sum(Parents_F%in%c(""))
            }
          
             
            #What are the LOD of the correct and incorrect assignments?
            if(Num.Correct[r]==2){
              LOD.Cor.Both[CountCorr.Both+1]<-as.numeric(Franz.results[r,7])
              CountCorr.Both<-CountCorr.Both+1
            }
            if(Num.Correct[r]==0){
              LOD.InCor.Both[(CountInCorr.Both+1)]<-as.numeric(Franz.results[r,7])
              CountInCorr.Both<-CountInCorr.Both+1
            }
            if(Num.Correct[r]==1){
              LOD.InCor.One[CountInCorr.One+1]<-as.numeric(Franz.results[r,7])
              CountInCorr.One<-CountInCorr.One+1
            }
         
             }#over r assignments
          
          SimLOD.CorrBoth[[s]]<-LOD.Cor.Both
          SimLOD.InCorrBoth[[s]]<-LOD.InCor.Both
          SimLOD.InCorrOne[[s]]<-LOD.InCor.One
          
          
          PropCorr<-(sum(Num.Correct)/((length(Num.Correct))*2))*100
          Res[c,s+3]<-PropCorr
          Res[c,s+3+nSim]<-if(sum(TypeI)>0) (sum(TypeI)/(length(TypeI)*2))*100 else 0
          Res[c,s+3+nSim*2]<-if(sum(TypeII)>0) (sum(TypeII)/(length(TypeII)*2))*100 else 0
          Res[c,s+3+nSim*3]<-paste(sprintf(mean(LOD.Cor.Both),fmt='%#.2f'),sprintf(sd(LOD.Cor.Both),fmt='%#.2f'),collapse="_")
          Res[c,s+3+nSim*4]<- paste(sprintf(mean(LOD.InCor.Both),fmt='%#.2f'),sprintf(sd(LOD.InCor.Both),fmt='%#.2f'),collapse="_")
          Res[c,s+3+nSim*5]<- paste(sprintf(mean(LOD.InCor.One),fmt='%#.2f'),sprintf(sd(LOD.InCor.One),fmt='%#.2f'),collapse="_")
          
        }#do it over s simulations
        c<-c+1
        MissDataLOD.CorrBoth[[MS]]<-SimLOD.CorrBoth
        MissDataLOD.InCorrBoth[[MS]]<-SimLOD.InCorrBoth
        MissDataLOD.InCorrOne[[MS]]<-SimLOD.InCorrOne
        
      }#over MS
      ErrorLOD.CorrBoth[[e]]<-MissDataLOD.CorrBoth
      ErrorLOD.InCorrBoth[[e]]<-MissDataLOD.InCorrBoth
      ErrorLOD.InCorrOne[[e]]<-MissDataLOD.InCorrOne
    }#over e error rates
    MarkerLOD.CorrBoth[[z]]<-ErrorLOD.CorrBoth
    MarkerLOD.InCorrBoth[[z]]<-ErrorLOD.InCorrBoth
    MarkerLOD.InCorrOne[[z]]<-ErrorLOD.InCorrOne
  }#over z panels
  write.csv(Res[,seq(1,(nSim*3)+3,1)],'FRANz_Summary.csv',row.names = F)
  
  #Plot the Accuracy
  
  Res<-as.data.frame(Res,stringsAsFactors=F)
  ResNum<-eval(parse(text=paste('cbind(Res[,(1:3)],',paste('as.numeric(Res[,',seq(4,nSim+3,1),'])',sep='',collapse=','),')',sep='',collapse=',')))
  colnames(ResNum)<-colnames(Res)[seq(1,nSim+3,1)]
  Res2<-reshape2::melt(ResNum, id.vars=c("Panel","Error", "MissDat")) 
  Res2[,1]<-factor(Res2[,1],levels(factor(Res2[,1]))[as.integer(factor(unique(Res2[,1])))])
  Res<-Res2
  
  
  return(Res)
}#SumFRANz
