#' Run FRANz through R for parentage on muliple PseudoBabies simulations
#'
#' This function facilitates the analysis of PseudoBabies simulations with the program FRANz
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @param ShowProgress Would you like a progress bar to track the analysis. Defaults to FALSE
#' @keywords FRANz parentage
#' @export
#' @examples
#' Rn.FRANz()

##### RUN FRANz #####   
Rn.FRANzHatch<-function(nSim,Miss_data=T,MissingVec,Geno_Error=T,ErrorVals,Markerfile='MarkerPanels.csv',ShowProgress=FALSE,...){
  wd<-getwd()
  if(file.exists(Markerfile)){
    markers<-read.csv(Markerfile, stringsAsFactors = F) 
  }else{
    stop(paste(Markerfile,"not found.",sep=" "))
  }
  
  
  # some housekeeping for numbers of reps to expect
  if(Miss_data==F) MissingVec<-0
  if(Geno_Error==F){
    ErrorVec<-0
  }else{
    ErrorVec<-seq(1,ErrorVals,1)
  }
  
  
  if(ShowProgress==T){
  total<-nrow(markers)*length(ErrorVec)*length(MissingVec)*nSim
  pb <- tcltk::tkProgressBar(title = "FRANz progress", min = 0,
                      max = total, width = 300)
  pg<-1
  }
  
  for (z in 1:nrow(markers)){  
    for (e in (1:length(ErrorVec))){
      Type_Error<-read.csv('Loci_error.csv',header=F,stringsAsFactors = F)[,c(1,e+1)]
      for(MS in (1:length(MissingVec)) ){
        Pmd<-MissingVec[MS]
        
        for (s in 1:nSim){
            for (a in 1:4){
          # Running FRANz#
          ParSampN<-length(grep(pattern="2000 ?",x= readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/',"MissDat_",Pmd,"/SimNum", s,"/","FRANzSim",s,"_",a,".dat",sep=""))))
          
          FRANzCmds=paste('--femrepro 1:1 --malerepro 1:1 --typingerror " ',max(Type_Error[,2]),'" --mintyped 1 --N ',ParSampN,' --n ', ParSampN, " --geofile", paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/',paste("MissDat_",Pmd,
          "/SimNum", s,"/","Geo.dat",sep="")," --parentsmaxdist 1",sep=''))
          cmd<- paste("FRANz", paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/',paste("MissDat_",Pmd,
                                            "/SimNum", s,"/","FRANzSim",s,"_",a,".dat ",sep=""),FRANzCmds,sep=""))
          if(Sys.info()[['sysname']]=="Windows"){
            shell(cmd,wait = TRUE)
            if(ShowProgress==T){
              tcltk::setTkProgressBar(pb=pb, value=pg,title = "FRANz progress",label=paste( round(pg/total*100, 0),
                                                   "% done"))
              pg<-pg+1
            }
            }else{
              system(cmd,wait = TRUE)
              if(ShowProgress==T){
                tcltk::setTkProgressBar(pb=pb, value=pg,title = "FRANz progress",label=paste( round(pg/total*100, 0),
                                                                                              "% done"))
                pg<-pg+1
              }
              }
          
          #Move output from FRANz to simulation folders
          filesOfInterest<-paste("^",c("parentage.csv","locisummary.txt","allelefreqs.dat","mcmc.log","mismatches.txt","pedigree.dot","simulation.txt","summary.txt","pedigree.dat"),"$",sep="")
          files2copy<-list.files()[unlist(sapply(1:length(filesOfInterest),function(x) grep(pattern=filesOfInterest[x],x=list.files())))]
          sapply(1:length(files2copy),function(y) file.rename(from=file.path(files2copy)[y],to=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/","SimNum",s,"/A",a,"_",files2copy[y],sep="")))
          sapply(1:length(files2copy),function(y) file.remove(file.path(files2copy[y])))
            }#over a
            }#over s
      }#over MS
    }#over error
  }#over panels
  if(ShowProgress==T) close(pb) #close the progress bar
}#Rn.FRANz

#' Resume running FRANz through R for parentage on muliple PseudoBabies simulations
#'
#' This function facilitates the analysis of PseudoBabies simulations with the program FRANz
#' @param nSims Number of replicate simulations
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @param ShowProgress Would you like a progress bar to track the analysis. Defaults to FALSE
#' @keywords FRANz parentage resume
#' @export
#' @examples
#' Resume.Rn.FRANz()
Resume.Rn.FRANzHatch<-function(nSim,Miss_data=T,MissingVec,Geno_Error=T,ErrorVals,Markerfile='MarkerPanels.csv',ShowProgress=FALSE){
  wd<-getwd()
  if(file.exists(Markerfile)){
    markers<-read.csv(Markerfile, stringsAsFactors = F) 
  }else{stop(paste(Markerfile,"not found.",sep=" "))
  }
  
  
  # some housekeeping for numbers of reps to expect
  if(Miss_data==F) MissingVec<-0
  if(Geno_Error==F){
    ErrorVec<-0
  }else{
    ErrorVec<-seq(1,ErrorVals,1)
  }
  
  if(ShowProgress==T){
    total<-nrow(markers)*length(ErrorVec)*length(MissingVec)*nSim
    pb <- tcltk::tkProgressBar(title = "FRANz progress", min = 0,
                               max = total, width = 300)
    pg<-1
  }
  
  
  #now lets figure out where we need to pickup the analysis
  Files2Search<-matrix(data=NA,ncol=5,nrow=nrow(markers)*length(MissingVec)*length(ErrorVec)*nSim)
  Files2Search[,1]<-rep(sapply(1:nrow(markers),function(x) paste("Panel_uSats",markers[x,2],"_SNPs",markers[x,3],sep="")),each=length(MissingVec)*length(ErrorVec)*nSim)
  Files2Search[,2]<-rep(paste("Error_",rep(ErrorVec,each=length(MissingVec)*nSim),sep=""),nrow(markers))
  Files2Search[,3]<-rep(paste("MissDat_",rep(MissingVec,each=nSim),sep=""),nrow(markers)*length(ErrorVec))
  Files2Search[,4]<-paste("SimNum",rep(1:nSim,times=nrow(markers)*length(MissingVec)*length(ErrorVec)),sep="")
  Files2Search[,5]<-paste("A",1:4,"_parentage.csv",sep="")
  #now just paste each row for the files to search!
  Files2SearchVec <- sapply(1:nrow(Files2Search),function(x) paste(Files2Search[x,],collapse="/"))
  FilesPres<-grep(pattern="parentage.csv",x=list.files(path=file.path(wd),recursive=T),value=T)
  IndexT<-min(which(unlist(sapply(1:length(Files2SearchVec),function(x) Files2SearchVec[x]%in%FilesPres)==0)))
  StartVals<-stringr::str_split(Files2SearchVec[IndexT],"/")
  
  MarkerNum.Res <- which(paste("Panel_uSats",markers[,2],"_SNPs",markers[,3],sep="")%in%StartVals[[1]][1])
  ErrorVec.Res <- which(ErrorVec%in%unlist(stringr::str_split(StartVals[[1]][2],"_"))[2])
  MissingNum.Res <- which(MissingVec%in%unlist(stringr::str_split(StartVals[[1]][3],"_"))[2])
  Sim.Res <- which(1:nSim%in%gsub("SimNum","",StartVals[[1]][4]))
  a.Res<- which(1:4%in%unlist(stringr::str_split(StartVals[[1]][5],"A|_"))[2])
  
  for (z in MarkerNum.Res:nrow(markers)){  
    for (e in (ErrorVec.Res:length(ErrorVec))){
      Type_Error<-read.csv('Loci_error.csv',header=F,stringsAsFactors = F)[,c(1,e+1)]
      for(MS in (MissingNum.Res:length(MissingVec)) ){
        Pmd<-MissingVec[MS]
        
        for (s in Sim.Res:nSim){
            for (a in a.Res:4){
          # Running FRANz#
          ParSampN<-length(grep(pattern="2000 ?",x= readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/',"MissDat_",Pmd,"/SimNum", s,"/","FRANzSim",s,"_",a,".dat",sep=""))))
          
          FRANzCmds=paste('--femrepro 1:1 --malerepro 1:1 --typingerror " ',max(Type_Error[,2]),'" --mintyped 1 --N ',ParSampN,' --n ', ParSampN,' -q',sep='')
          cmd<- paste("FRANz", paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],'/',paste("MissDat_",Pmd,
                                                                                                                               "/SimNum", s,"/","FRANzSim",s,"_",a,".dat ",sep=""),FRANzCmds,sep=""))
          if(Sys.info()[['sysname']]=="Windows"){
            shell(cmd,wait = TRUE)
            if(ShowProgress==T){
              tcltk::setTkProgressBar(pb, pg, label=paste( round(pg/total*100, 0),
                                                    "% done"))
              pg<-pg+1
            }
          }else{
            system(cmd,wait = TRUE)
            if(ShowProgress==T){
              tcltk::setTkProgressBar(pb, pg, label=paste( round(pg/total*100, 0),
                                                    "% done"))
              pg<-pg+1
            }
          }
          
          #Move output from FRANz to simulation folders
          filesOfInterest<-paste("^",c("parentage.csv","locisummary.txt","allelefreqs.dat","mcmc.log","mismatches.txt","pedigree.dot","simulation.txt","summary.txt","pedigree.dat"),"$",sep="")
          files2copy<-list.files()[unlist(sapply(1:length(filesOfInterest),function(x) grep(pattern=filesOfInterest[x],x=list.files())))]
          sapply(1:length(files2copy),function(y) file.rename(from=file.path(files2copy)[y],to=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],'/','Error_',ErrorVec[e],"/MissDat_",Pmd,"/","SimNum",s,"/A",a,"_",files2copy[y],sep="")))
          sapply(1:length(files2copy),function(y) file.remove(file.path(files2copy[y])))
          
          
        }#over a
        a.Res<-1
        }#over S
        Sim.Res<-1
      }#over MS
      MissingNum.Res<-1
    }#over error
    ErrorVec.Res<-1
  }#over panels
  if(ShowProgress==T) close(pb) #close the progress bar
  MarkerNum.Res<-1
}#Resume.Rn.FRANzHatch
