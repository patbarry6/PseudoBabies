#' Run Colony2 through R for parentage on muliple PseudoBabies simulations
#'
#' This function facilitates the analysis of PseudoBabies simulations with the program Colony
#' @param ColonyDir Directory where the Colony program can be called from
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @keywords Colony parentage
#' @export
#' @examples
#' Rn.Colony()

Rn.Colony<-function(ColonyDir,nSim,Miss_data=FALSE,MissingVec=c(0,2.5),Geno_Error=FALSE,ErrorVals,Markerfile,ShowProgress=FALSE, ...){
                    wd<-getwd()
                    
                    if(file.exists(Markerfile)){
                    markers<-read.csv(Markerfile,stringsAsFactors = F) 
                    }else{stop(paste(Markerfile,"not found.",sep=" "))
                      }
                    
                    # some housekeeping for numbers of reps to expect
                    if(Miss_data==F) MissingVec<-0
                    if(Geno_Error==F){
                      ErrorVec<-0
                    }else{ErrorVec<-seq(from=1,to=ErrorVals,by=1)
                    }
                    
                    if(ShowProgress==T){
                      total<-nrow(markers)*length(ErrorVec)*length(MissingVec)*nSim
                      pb <- tcltk::tkProgressBar(title = "Colony progress", min = 0,
                                                 max = total, width = 300)
                      pg<-1
                    }
                    
                    for (z in 1:nrow(markers)){  
                      for (e in (1:length(ErrorVec))){
                        Type_Error<-read.csv("Loci_error.csv",header=F,stringsAsFactors = F)[,c(1,e+1)]
                        for(MS in (1:length(MissingVec)) ){
                          Pmd<-MissingVec[MS]
                          
                          for (s in 1:nSim){
                            
                            #run some commands
                            if(Sys.info()[['sysname']]=="Windows"){
                            cmd<-paste("C:\\WINDOWS\\system32\\cmd.exe /C ",ColonyDir,"/colony2s.exe IFN:",
                                             wd,"/","Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                                             "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                                             "/SimNum",s,"/ColonyInput.dat",sep="")
                            system(cmd)
                            }else{
                              cmd<-paste(ColonyDir,"/colony2s.out IFN:",
                                         wd,"/","Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                                         "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                                         "/SimNum",s,"/ColonyInput.dat",sep="")
                              system(cmd)
                              }
                            
                           
                            #move the output
                            files2copy<-list.files()[grep(pattern=paste("Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                                          "_Error",sep=""),x=list.files())]
                            file.copy(from=files2copy,to=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],
                                       "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                                       "/SimNum",s,'/',sep=""),recursive = T)
                            file.remove(files2copy)
                            
                            if(ShowProgress==T){
                              tcltk::setTkProgressBar(pb, pg, label=paste( round(pg/total*100, 0),
                                                                    "% done"))
                              pg<-pg+1
                            }
                            
                        }#over S
                        }#over MS
                      }#over error
                    }#over panels
                    if(ShowProgress==T) close(pb) #close the progress bar
                     }#Rn.Colony


#' Resume running Colony2 through R for parentage on muliple PseudoBabies simulations
#'
#' This function facilitates the analysis of PseudoBabies simulations with the program Colony
#' @param ColonyDir Directory where the Colony program can be called from
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @keywords Colony parentage
#' @export
#' @examples
#' Resume.Rn.Colony()
Resume.Rn.Colony<-function(ColonyDir,nSim,Miss_data=FALSE,MissingVec=c(0,2.5),Geno_Error=FALSE,ErrorVals,Markerfile,ShowProgress=FALSE){
  wd<-getwd()
  
  
  if(file.exists(Markerfile)){
    markers<-read.csv(Markerfile,stringsAsFactors = F) 
  }else{stop(paste(Markerfile,"not found.",sep=" "))
  }
  
  # some housekeeping for numbers of reps to expect
  if(Miss_data==F) MissingVec<-0
  if(Geno_Error==F){
    ErrorVec<-0
  }else{ErrorVec<-seq(from=1,to=ErrorVals,by=1)
  }
  
  
  if(ShowProgress==T){
    total<-nrow(markers)*length(ErrorVec)*length(MissingVec)*nSim
    pb <- tcltk::tkProgressBar(title = "Colony progress", min = 0,
                               max = total, width = 300)
    pg<-1
  }
  
  #now lets figure out where we need to pickup the analysis
  Files2Search<-matrix(data=NA,ncol=4,nrow=nrow(markers)*length(MissingVec)*length(ErrorVec)*nSim)
  Files2Search[,1]<-rep(sapply(1:nrow(markers),function(x) paste("Panel_uSats",markers[x,2],"_SNPs",markers[x,3],sep="")),each=length(MissingVec)*length(ErrorVec)*nSim)
  Files2Search[,2]<-rep(paste("Error_",rep(ErrorVec,each=length(MissingVec)*nSim),sep=""),nrow(markers))
  Files2Search[,3]<-rep(paste("MissDat_",rep(MissingVec,each=nSim),sep=""),nrow(markers)*length(ErrorVec))
  Files2Search[,4]<-paste("SimNum",rep(1:nSim,times=nrow(markers)*length(MissingVec)*length(ErrorVec)),sep="")
  #now just paste each row for the files to search!
  Files2SearchVec <- sapply(1:nrow(Files2Search),function(x) paste(Files2Search[x,],collapse="/"))
  IndexT<-min(which(unlist(sapply(1:length(Files2SearchVec),function(x) length(grep(pattern="ParentPair",x=list.files(path=file.path(paste("./",Files2SearchVec[x],sep="")))))==0))==T))
  StartVals<-stringr::str_split(Files2SearchVec[IndexT],"/") 
  
  MarkerNum.Res <- which(Files2Search[,1]%in%StartVals[[1]][1])
  ErrorVec.Res <- which(ErrorVec%in%unlist(stringr::str_split(StartVals[[1]][2],"_"))[2])
  MissingNum.Res <- which(MissingVec%in%unlist(stringr::str_split(StartVals[[1]][3],"_"))[2])
  Sim.Res <- which(1:nSim%in%gsub("SimNum","",StartVals[[1]][4]))
  
  for (z in MarkerNum.Res:nrow(markers)){  
    for (e in (ErrorVec.Res:length(ErrorVec))){
      Type_Error<-read.csv("Loci_error.csv",header=F,stringsAsFactors = F)[,c(1,e+1)]
      for(MS in (MissingNum.Res:length(MissingVec)) ){
        Pmd<-MissingVec[MS]
        
        for (s in Sim.Res:nSim){
          
          #run some commands
          if(Sys.info()[['sysname']]=="Windows"){
            cmd<-paste("C:\\WINDOWS\\system32\\cmd.exe /C ",ColonyDir,"/colony2s.exe IFN:",
                       wd,"/","Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                       "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                       "/SimNum",s,"/ColonyInput.dat",sep="")
            }else{
            cmd<-paste(ColonyDir,"/colony2s.out IFN:",
                       wd,"/","Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                       "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                       "/SimNum",s,"/ColonyInput.dat",sep="")
          }
          
          system(cmd)
          #move the output
          files2copy<-list.files()[grep(pattern=paste("Panel_uSats",markers[z,2],"_SNPs",markers[z,3],
                                                      "_Error",sep=""),x=list.files())]
          file.copy(from=files2copy,to=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],
                                             "/Error_",ErrorVec[e],"/MissDat_",MissingVec[MS],
                                             "/SimNum",s,'/',sep=""),recursive = T)
          file.remove(files2copy)
          
          if(ShowProgress==T){
            tcltk::setTkProgressBar(pb, pg, title = "Colony progress", label=paste( round(pg/total*100, 0),
                                                  "% done"))
            pg<-pg+1
          }
        }#over S
        Sim.Res<-1
      }#over MS
      MissingNum.Res<-1
    }#over error
    ErrorVec.Res<-1
  }#over panels
  MarkerNum.Res<-1
  if(ShowProgress==T) close(pb) #close the progress bar
}#Resume.Rn.Colony
