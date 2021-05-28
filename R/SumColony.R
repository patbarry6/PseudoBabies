#' Summarize the results from Colony2 analysis for parentage on muliple PseudoBabies simulations
#'
#' This function summarizes the results from multiple runs of Colony2 for parentage inference from PseudoBabies simulations.
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVals How many error values were specified in Loci_error.csv? Required if Geno_Error=TRUE. Defaults to 1 percent.
#' @param Cutoff cutoff value for probability value accepting parentage.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @keywords Colony parentage
#' @export
#' @examples
#' SumColony()

####Summarize Colony####
SumColony<-function(nSim,Miss_data=T,MissingVec=c(0,2.5),Geno_Error=T,ErrorVals=1,Cutoff=.8,Markerfile='MarkerPanels.csv'){
  library(ggplot2)  
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
    
    Res<-matrix(data="NA",ncol=(nSim)*3+3,nrow=ErrorVals*length(MissingVec)*nrow(markers))
    Res[,1]<-rep(Panel,each=ErrorVals*length(MissingVec))
    Res[,2]<-rep(ErrorVec,each=length(MissingVec))
    Res[,3]<-if(length(MissingVec)>1) rep(MissingVec,nrow(markers)) else rep(MissingVec)
    colnames(Res)<-c('Panel','Error','MissDat',paste('Sim',seq(1,nSim,1),sep=''),as.vector(sapply(1:nSim,function(x) paste(paste("Sim",seq(1,nSim,1),sep="")[x],paste("Type",c("I","II"),sep=""),sep=""))))
    
    for (z in 1:nrow(markers)){
        mark<-markers[z,-1]
        for(e in (1:length(ErrorVec))){
            for(MS in (1:length(MissingVec)) ){
                Pmd<-MissingVec[MS]
                for (s in 1:nSim){
                    
                    Colony.results<- read.csv(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/SimNum",s,'/Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'_','Error',ErrorVec[e],"_MissDat",Pmd,"_SimNum",s,".ParentPair",sep=""),stringsAsFactors=F)
                    
                    Answer<-read.csv( file = paste(wd,'/SimParents/','Truth_Sim',s,".csv",sep=""),stringsAsFactors=F)
                    Answernames<-Answer[,1]
                    Answer<-stringr::str_split(Answer[,2],'&')
                    names(Answer)<-Answernames
                    
                    Num.Correct<-vector()
                    TypeI<-vector()
                    TypeII<-vector()
                    #These are just found placement counting
                    CountCorr<-0
                    CountInCorr<-0
                    Files2Write<-c('PropCorr')
                    #Colony gives multiple parent pairs if there is ambiguity in the
                    #assignment of offspring. We can either take the maximum likelihood
                    #assignement or look to see if it included the correct parents.
                    #here we take the maximum likelihood pair inferred
                    for (r in 1:length(Answer)){
                      Offspring_F<-names(Answer)[r]
                      W<-which(Colony.results[,1]==Offspring_F)
                      Parents_F<-Colony.results[W,c(2,3)][which(max(Colony.results[W,4])==Colony.results[W,4]),]
                      TRUTH.S<-eval(parse(text=paste("Answer$'",Offspring_F,"'",sep="")))
                      #if the Probability is >cutoff value then store number correct
                      if(Colony.results[W,4][which(max(Colony.results[W,4])==Colony.results[W,4])][1] >= Cutoff){ #what about ties to this one...
                        Num.Correct[r]<-max(sapply(1:nrow(Parents_F),function(x) sum(Parents_F[x,]%in%TRUTH.S)))
                        TypeI[r]<-sum(!(grep(pattern="[0-9]+",x=Parents_F,value=T)%in%TRUTH.S))
                        TypeII[r]<-sum(Parents_F%in%c("*","#"))
                      }else{
                        Num.Correct[r]<-0
                        TypeI[r]<-0
                        TypeII[r]<-0
                      }
                      #num of parents right?
                    }#over r rows
                    
                    PropCorr<-(sum(Num.Correct)/((length(Num.Correct))*2))*100
                    TypeIProp<-(sum(TypeI)/((length(TypeI))*2))*100
                    TypeIIProp<-(sum(TypeII)/((length(TypeII))*2))*100
                    Res[c,s+3]<-PropCorr
                    Res[c,s+3+nSim]<-TypeIProp
                    Res[c,s+3+nSim*2]<-TypeIIProp
                    
                    
                    
                }#do it over s simulations
                c<-c+1
                
            }#over MS
        }#over e error rates
    
      write.csv(x=Res,file='ColonySummary.csv',row.names = F)
      }#over z panels
    
    
    Res<-as.data.frame(Res,stringsAsFactors=F)
    ResNum<-eval(parse(text=paste('cbind(Res[,(1:3)],',paste('as.numeric(Res[,',seq(4,(3*nSim)+3,1),'])',sep='',collapse=','),')',sep='',collapse=',')))
    colnames(ResNum)<-colnames(Res)
    Res2<-reshape2::melt(ResNum[,1:(nSim+3)], id.vars=c("Panel","Error", "MissDat"))
    Res2[,1]<-factor(Res2[,1],levels(factor(Res2[,1]))[as.integer(factor(unique(Res2[,1])))])
    
    
    AssignmentSummaryColony<-Res2
    return(AssignmentSummaryColony)
}#SumColony

