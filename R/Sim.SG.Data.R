#' Simulate genotypes over a single year
#'
#' This function allows you to simulate multilocus genotypes under the assumption of completely random mating for a single year
#' @param Founders .csv file containing the genotypes to use for initiating the simulation. Use function CompleteGenotypes to remove missing data before performing simulations.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @param NumSNPs number of SNP markers.
#' @param NumuSats number of microsatellite markers.
#' @param NB number of crosses to perform. Defaults to 100.
#' @param OffspringDist What distribution should we use for the number of offspring produced for each cross. Options are Uniform (specify a min and a max), NegativeBinomial (specify a prob for the scale parameter see \code{\link[stats]{NegBinomial}} options), and Poisson (specify lambda). 
#' @param SexInfo Does your dataset have information on sex?
#' @param MatingStr Should the simulation allow for polygamy, polygyny, polandry, or monogamy?
#' @param StartYear What year should the simulation start? Defaults to current year.
#' @param nSim How many replicate simulations would you like to do? Defaults to 10.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVals How many different values of genotyping error for simulations. Required if Geno_Error=TRUE. Values specified in Loci_error.csv
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE.
#' @param Programs Which programs would you like to create input files for? Defaults to c('FRANz','Colony')
#' @param SexInf Use the sex of individuals in the parentage inference program
#' @param perl.dir Full path to .csv.pl perl script distributed with FRANz for formatting input files. Required when FRANz is specified as a Program argument
#' @param GPHomogeneity A logical to indicate if genepop should be run to compare the allele frequency distribution for all simulated datasets with the founders dataset. Defaults to FALSE
#' @importFrom magrittr "%>%"
#' @keywords simulation parentage genetics
#' @export
#' @examples
#' Sim_SG_Data()
##### Simulate Single Generation of Data ####
Sim.SG.Data<-function(Founders='FoundersCompleteGenotypes.csv', Markerfile='MarkerPanels.csv', NumSNPs=93,NumuSats=12,NB=650,StartYear=as.numeric(format(Sys.Date(), "%Y")),nSim=10, Geno_Error=F,ErrorVals,Miss_data=F ,MissingVec ,Programs=c('FRANz','Colony'),SexInf,
    perl.dir,GPHomogeneity=FALSE,OffspringDist="Poisson",SexInfo=F,lambda,Umin=0,Umax=1,nBprob,nBsize,
    ztnbprob,ztnbsize,MatingStr="Polygyny"){
  
  #check for the founders file
  if(!file.exists(Founders)){
      stop("Founders file not found.")
  }
  
  if(!file.exists(Markerfile)){
      stop("Marker file not found.")
     
  }
  
  if(Sys.info()[['sysname']]=="Windows"){
    EOLcharacter<-"\r\n"
  }else{EOLcharacter<-"\n"
  }
  
  wd<-getwd()
  dat<-read.csv(Founders,stringsAsFactors = F)
  if(SexInfo==F){dat$Sex<-sample(c("M","F"),size=nrow(dat),replace = T)}

  markers<-read.csv(Markerfile,stringsAsFactors = F)
  nloci<-NumSNPs+NumuSats
  loci<-colnames(dat)[seq(from=2,to=nloci*2+1,by=2)]
  Count<-paste('(',seq(2,nloci*2,2),':',seq(3,nloci*2+1,2),')',sep='')
  
  if(GPHomogeneity==T){
    GPin<-list()
    GPin[[1]]<-dat[,-c(nloci*2+2,nloci*2+3)]
  }
  
  if(Geno_Error==T){ErrorVec<-seq(1,ErrorVals,1)
  }else{ErrorVec<-0
  }
  if(Miss_data==F){MissingVec<-0}
  
  Nb<-NB
  ####Set up some folders to hold our results#### 
  dir.create(file.path("SimParents"))
  
  Panel<-paste("Panel_uSats",markers[,2],"_SNPs",markers[,3],sep="")
  sapply(1:length(Panel),function(x) dir.create(file.path(Panel[x])))
  path1<-paste("Error_",ErrorVec,sep="")
  sapply(1:length(path1),function(x) dir.create(file.path(path1[x])))
  path2<-paste("MissDat_",MissingVec,sep="")
  sapply(1:length(path2),function(x) dir.create(file.path(path2[x])))
  path3<-paste("SimNum",seq(1,nSim,1),sep="")
  sapply(1:length(path3),function(x) dir.create(file.path(path3[x])))

  
  #Now lets move the folders where they belong
  sapply(1:length(path2),function(x) sapply(1:length(path3),function(y) file.copy(from=path3[y],to=path2[x],recursive=T)))
  sapply(1:length(path1),function(x) sapply(1:length(path2),function(y) file.copy(from=path2[y],to=path1[x],recursive=T)))
  sapply(1:length(Panel),function(x) sapply(1:length(path1),function(y) file.copy(from=path1[y],to=Panel[x],recursive=T)))
  
  #clean up the working directory
  unlink(path1,recursive=T)
  unlink(path2,recursive=T)
  unlink(path3,recursive=T)
  
  #### Simulate Data ####
  for (s in 1:nSim){
    
    DataList.names <- seq(StartYear,StartYear+1)
    DataList <- vector("list", length(DataList.names))
    names(DataList) <- DataList.names
    
   
    DataList[[1]]<-dat  # Let's put the Founders in a list with year entry being the death year!
    
    parents<-as.data.frame(DataList[[1]],stringsAsFactors=F) #Pull out the year we are using for parents
    
    Nb_sire<-if(nrow(parents[which(parents$Sex=="M"),]) < Nb ) nrow(parents[which(parents$Sex=="M"),]) else Nb
    Nb_dam<-if(nrow(parents[which(parents$Sex=="F"),]) < Nb ) nrow(parents[which(parents$Sex=="F"),]) else Nb
    
    if(MatingStr=="Polygamy"){
      #How many offspring will each set of parents produce?
      if(OffspringDist == "Uniform"){
        OffspringNums<-round(runif(Nb,min=Umin,max=Umax))
      }else if(OffspringDist == "NegativeBinomial") {
        OffspringNums<-rnbinom(n=Nb, size=nBsize, prob=nBprob)
      }else if(OffspringDist == "ZeroTruncatedNegativeBinomial") {
          OffspringNums<-actuar::rztnbinom(n=Nb, size=ztnbsize, prob=ztnbprob)
      }else if(OffspringDist == "Poisson") {
        OffspringNums<-rpois(n=Nb,lambda=lambda)
      }
      
      Sires<-sample(parents[which(parents$Sex=="M"),1],Nb,replace=T)
      Dams<-sample(parents[which(parents$Sex=="F"),1],Nb,replace=T)
      PAR<-cbind(Sires,Dams)
      
    } else if(MatingStr=="Polygyny"){ 
      #How many offspring will each set of parents produce?
      if(OffspringDist == "Uniform"){
        OffspringNums<-round(runif(Nb_sire,min=Umin,max=Umax))
      }else if(OffspringDist == "NegativeBinomial") {
        OffspringNums<-rnbinom(n=Nb_sire, size=nBsize, prob=nBprob)
      }else if(OffspringDist == "ZeroTruncatedNegativeBinomial") {
          OffspringNums<-actuar::rztnbinom(n=Nb_sire, size=ztnbsize, prob=ztnbprob)
      }else if(OffspringDist == "Poisson") {
        OffspringNums<-rpois(n=Nb_sire,lambda=lambda)
      }
      
      Sires<-sample(parents[which(parents$Sex=="M"),1],Nb_sire,replace=F)
      Dams<-sample(parents[which(parents$Sex=="F"),1],Nb_sire,replace=T)
      PAR<-cbind(Sires,Dams)
      
    } else if(MatingStr=="Polyandry"){
      #How many offspring will each set of parents produce?
      if(OffspringDist == "Uniform"){
        OffspringNums<-round(runif(Nb_dam,min=Umin,max=Umax))
      }else if(OffspringDist == "NegativeBinomial") {
        OffspringNums<-rnbinom(n=Nb_dam, size=nBsize, prob=nBprob)
      }else if(OffspringDist == "ZeroTruncatedNegativeBinomial") {
          OffspringNums<-actuar::rztnbinom(n=Nb_dam, size=ztnbsize, prob=ztnbprob)
      } else if(OffspringDist == "Poisson") {
        OffspringNums<-rpois(n=Nb_dam,lambda=lambda)
      }
      
      Sires<-sample(parents[which(parents$Sex=="M"),1],Nb_dam,replace=T)
      Dams<-sample(parents[which(parents$Sex=="F"),1],Nb_dam,replace=F)
      PAR<-cbind(Sires,Dams)
      
    } else if(MatingStr=="Monogomy"){
      Nb<-if(Nb_sire < Nb | Nb_dam < Nb) min(Nb_sire,Nb_dam) else Nb
      #How many offspring will each set of parents produce?
      if(OffspringDist == "Uniform"){
        OffspringNums<-round(runif(Nb,min=Umin,max=Umax))
      }else if(OffspringDist == "NegativeBinomial") {
        OffspringNums<-rnbinom(n=Nb, size=nBsize, prob=nBprob)
      }else if(OffspringDist == "ZeroTruncatedNegativeBinomial") {
          OffspringNums<-actuar::rztnbinom(n=Nb, size=ztnbsize, prob=ztnbprob)
      }else if(OffspringDist == "Poisson") {
        OffspringNums<-rpois(n=Nb,lambda=lambda)
      }
      
      Sires<-sample(parents[which(parents$Sex=="M"),1],Nb,replace=F)
      Dams<-sample(parents[which(parents$Sex=="F"),1],Nb,replace=F)
      PAR<-cbind(Sires,Dams)
    } else cat("Houston we have a problem.")
  
    
    
    #who are the parents going to be?
    
    PAR<-PAR[which(OffspringNums!=0),]
    OffspringNums<-OffspringNums[which(OffspringNums!=0)]
    
    PRO<- matrix(data=NA, ncol=(nloci*2), nrow=sum(OffspringNums)) #Progeny
    CorrectPAR<-vector()
    CorrectPAR<-rep(paste(PAR[,1],'&',PAR[,2],sep=''),OffspringNums)
    
    # Fill in offspring genotypes from PAR
    Samp<-lapply(1:nrow(PAR), function(x) eval(parse(text=paste('cbind(',paste(paste('sample(parents[which(parents[,1]%in%PAR[x,]),',Count,'],OffspringNums[x],replace=T)',sep=''),collapse=','),')',sep=''))))
    
    for (i in 1:nrow(PAR)){
      Samp.trans<-t(Samp[[i]])
      PRO[((sum(OffspringNums[1:i])-OffspringNums[i]+1):(sum(OffspringNums[1:i]))),seq(1,nloci*2,2)]<-Samp.trans[,1]
      PRO[((sum(OffspringNums[1:i])-OffspringNums[i]+1):(sum(OffspringNums[1:i]))),seq(2,nloci*2,2)]<-Samp.trans[,2]
    }
    
    # When will each reproduce and die
    # Single generation doesn't pull from maturity file, just one year in the future
    AgeDeath<-DataList.names[1]+1
    
    #make some new snazzy names for each so we can keep track of when they were born
    Names<-paste(DataList.names[1],'_',seq(1,nrow(PRO),1),sep='') #This used to be paste('BY',DataList.names[1],'_',seq(1,nrow(PRO),1),sep='')
    NewDat<-as.matrix(cbind(Names,PRO,sample(x=c("M","F"),nrow(PRO),replace=T),rep(DataList.names[1],nrow(PRO)),AgeDeath,CorrectPAR))
    
    colnames(NewDat)<-c(colnames(dat),'Parents')
    
    # now I need to store each individual in their own year according to their death year
    DataList[[2]]<-NewDat
    
    # Now to get all the data we need for the analysis
    Input<-rbind(DataList[[1]],(DataList[[2]][,(-ncol(DataList[[2]]))]))
    #now to move last two columns to the front of the input file
    InputFull<-cbind(as.character(Input[,1]),c(rep(2000,nrow(DataList[[1]])),rep(2001,nrow(DataList[[2]]))),Input[,seq(2,nloci*2+1)],Input[,(ncol(Input)-2)])
    
    colnames(InputFull)<-c(colnames(dat)[1],'Year',colnames(dat)[c(-1,-(length(colnames(dat))-2):-(length(colnames(dat))))],"Sex")
    
    TRUTH<-DataList[[2]][,c(1,ncol(DataList[[2]]))]
    
    if(GPHomogeneity==T){GPin[[s+1]]<-NewDat[,-c((nloci*2+2):(nloci*2+4))]}
    
    write.table(TRUTH, file = paste(wd,"/SimParents/","Truth_","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
    
    if(Geno_Error==T){
      for (e in 1:length(ErrorVec)){
        Input2<-as.matrix(InputFull)  
        # Adding Genotyping Error #
        #Lets add some genotyping error to Input2
        #read in error file
        Type_Error<-read.csv('Loci_error.csv',header=F,stringsAsFactors = F)[,c(1,e+1)] #this gives me error by locus
        Type_Error$Num<-round(nrow(Input2)*Type_Error[,2],0) #how many individuals have bad genotypes?
        #sample genotypes by locus that get typing error
        #lets just apply a simple +/-4bp error
        
        for (LR in 1:nloci){
          if(Type_Error[LR,3]>0){
            temp<-Input2[,((LR*2+1):(LR*2+2))]
            Rplc<-sample(length(temp),Type_Error[LR,3],replace=F)
            #if the locus is fixed skip it
            if(length(unique(as.vector(temp)))==1){
              writeLines(paste('Locus',LR,'is fixed for allele',unique(as.vector(temp)),'\n',sep=' '))
            }else{  
              #if it is a SNP call it the opposite type
              if(length(unique(as.vector(temp)))<3){
                nal<-unique(as.vector(temp))
                repVec<-vector()
                temp[is.na(temp)]<-0
                temp[Rplc]<-unlist(lapply(1:length(Rplc), function (x) temp[sample(which(nal!=temp[Rplc][x]),1,F)]) )#
                Input2[,((LR*2+1):(LR*2+2))]<-temp 
              }#if the locus is a SNP
              
              
              if(length(unique(as.vector(temp)))>2){
                #if it is a uSAt add/subtract a repeat
                as<-sample(x=c(-4,4),size=length(Rplc),replace=T)
                temp[Rplc]<-as.numeric(temp[Rplc])+as #Super simple error - misscoring results in -4/+4 bases
                Input2[,((LR*2+1):(LR*2+2))]<-temp
              }
            }#if typing error
          }#else for fixed alleles
        }#over loci
        InputErr<-as.matrix(Input2)
        #if we have genotyping error we might have missing data we might not    
        if(Miss_data==T){
          for (M in 1:length(MissingVec)){
            Input4<-InputErr
            # if the values are missing it is most likely
            # a complete genotype and not just one allele of a genotype
            Pmd<-MissingVec[M]
            Nmd<-round(((nrow(Input4)*(ncol(Input4)-2))/2)*(0.01*Pmd),0) # of full genotypes that should get NAs
            temp<-Input4[,seq(from=3,to=nloci*2+2,by=2)] #subset our data
            Rplc<-sample(length(temp),Nmd,replace=F) #which individuals?
            temp[Rplc]<-888888
            Input4[,seq(from=3,to=nloci*2+2,by=2)]<-temp
            Input4[which(Input4==888888)+nrow(Input4)]<-888888
            Input4[Input4==888888]<-0  
            if(SexInf==F){
              for (z in 1:nrow(markers)){
                mark<-markers[z,-1]
                if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)])}  
                if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)])}  
                if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)])}  
                
                if('FRANz'%in% Programs){
                  write.table(DatExp, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                  }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                  }
                }# if FRANz
                
                if('Colony' %in% Programs){
                  Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                  Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                  Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  
                  #running colony from cmd line
                  Num_Offspring<-nrow(Col_Offspring)
                  Num_Canidates<-nrow(Col_Parents)
                  Num_Loci<-sum(markers[z,2],markers[z,3])
                  
                  Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                  Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                  Col_Markers[2,]<-0#codominant
                  Col_Markers[3,]<-0#dropout
                  Col_Markers[4,]<-0#error rate
                  Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                  
                  AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec[e],"_MissDat",MissingVec[M],"_SimNum",s,sep='')
                  
                  readLines(system.file("extdata/colony2_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                    gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                    gsub(pattern="NUM_CANDIDATES1",replacement=Num_Canidates,x=.)%>%
                    gsub(pattern="NUM_CANDIDATES2",replacement=Num_Canidates,x=.)%>%
                    gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                    gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                    gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                    gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                    gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                    gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                    gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                    writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput.dat",sep=""))
                }
              }#over z markers
            } else if (SexInf==TRUE){
              
              for (z in 1:nrow(markers)){
                mark<-markers[z,-1]
                if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)])}  
                if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)])}  
                if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input4)])}  
                colnames(DatExp)[ncol(DatExp)]<-"Sex"
                if('FRANz'%in% Programs){
                  DatExp2<-DatExp[,c(1,2,ncol(DatExp),seq(3,ncol(DatExp)-1,1))]
                  write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                  }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                  }
                }# if FRANz
                
                if('Colony' %in% Programs){
                  Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                  Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                  Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  Col_Sires<-Col_Parents%>%
                    filter(.,Sex=="M")%>%
                    select(.,-Sex)%>%
                    {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                  Col_Dams<-Col_Parents%>%
                    filter(.,Sex=="M")%>%
                    select(.,-Sex)%>%
                    {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                  #running colony from cmd line
                  Num_Offspring<-nrow(Col_Offspring)
                  Num_Sires<-nrow(Col_Parents[Col_Parents$Sex=="M",])
                  Num_Dams<-nrow(Col_Parents[Col_Parents$Sex=="F",])
                  Num_Loci<-sum(markers[z,2],markers[z,3])
                  
                  Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                  Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                  Col_Markers[2,]<-0#codominant
                  Col_Markers[3,]<-0#dropout
                  Col_Markers[4,]<-0#error rate
                  Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                  
                  AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec[e],"_MissDat",MissingVec[M],"_SimNum",s,sep='')
                  
                  readLines(system.file("extdata/colony2_sex_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE))%>%
                    gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                    gsub(pattern="NUM_CANDIDATES1",replacement=Num_Sires,x=.)%>%
                    gsub(pattern="NUM_CANDIDATES2",replacement=Num_Dams,x=.)%>%
                    gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                    gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                    gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Sires,x=.)%>%
                    gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Dams,x=.)%>%
                    gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                    gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                    gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                    writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput.dat",sep=""))
                }#if Colony is in program files
              }#over z markers  
            } #if SexInf==T
          }#over MS Miss_data  
          
        }else{ #if missing data = F
          Input4<-InputErr
          MissingVec<-0
          if(SexInf==F){
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)])}  
              
              if('FRANz'%in% Programs){
                write.table(DatExp, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                
                #running colony from cmd line
                Num_Offspring<-nrow(Col_Offspring)
                Num_Canidates<-nrow(Col_Parents)
                Num_Loci<-sum(markers[z,2],markers[z,3])
                
                Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                Col_Markers[2,]<-0#codominant
                Col_Markers[3,]<-0#dropout
                Col_Markers[4,]<-0#error rate
                Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                
                AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec[e],"_MissDat",MissingVec,"_SimNum",s,sep='')
                
                readLines(system.file("extdata/colony2_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                  gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES1",replacement=Num_Canidates,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES2",replacement=Num_Canidates,x=.)%>%
                  gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                  gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                  gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                  gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                  gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                  gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                  gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput.dat",sep=""))
              }
            }#over z markers
          }else if (SexInf==TRUE){
            
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input4)])}  
              colnames(DatExp)[ncol(DatExp)]<-"Sex"
              if('FRANz'%in% Programs){
                DatExp2<-DatExp[,c(1,2,ncol(DatExp),seq(3,ncol(DatExp)-1,1))]
                write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Sires<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                Col_Dams<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                #running colony from cmd line
                Num_Offspring<-nrow(Col_Offspring)
                Num_Sires<-nrow(Col_Parents[Col_Parents$Sex=="M",])
                Num_Dams<-nrow(Col_Parents[Col_Parents$Sex=="F",])
                Num_Loci<-sum(markers[z,2],markers[z,3])
                
                Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                Col_Markers[2,]<-0#codominant
                Col_Markers[3,]<-0#dropout
                Col_Markers[4,]<-0#error rate
                Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                
                AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec[e],"_MissDat",MissingVec,"_SimNum",s,sep='')
                
                readLines(system.file("extdata/colony2_sex_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                  gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES1",replacement=Num_Sires,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES2",replacement=Num_Dams,x=.)%>%
                  gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                  gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                  gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Sires,x=.)%>%
                  gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Dams,x=.)%>%
                  gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                  gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                  gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput.dat",sep=""))
              }
            }#over z markers
          } else {
            cat("Issue with SexInf")
          }
          
          
        } # else 
      }# over ErrorVec
    }else{ #if Geno_Error is false
      ErrorVec<-0
      
      #no genotyping error
      if(Miss_data==T){
        for (M in 1:length(MissingVec)){
          Input2<-as.matrix(InputFull)
          # if the values are missing it is most likely
          # a complete genotype and not just one allele of a genotype
          Pmd<-MissingVec[M]
          Nmd<-round(((nrow(Input2)*(ncol(Input2)-2))/2)*(0.01*Pmd),0) # of full genotypes that should get NAs
          temp<-as.matrix(Input2[,seq(from=3,to=nloci*2+2,by=2)]) #subset our data
          Rplc<-sample(length(temp),Nmd,replace=F) #which individuals?
          temp[Rplc]<-888888
          Input2[,seq(from=3,to=nloci*2+2,by=2)]<-temp#
          Input2[which(Input2==888888)+nrow(Input2)]<-888888
          Input2[Input2==888888]<-0  
          if(SexInf==F){
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)])}  
              
              if('FRANz'%in% Programs){
                write.table(DatExp, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                
                #running colony from cmd line
                Num_Offspring<-nrow(Col_Offspring)
                Num_Canidates<-nrow(Col_Parents)
                Num_Loci<-sum(markers[z,2],markers[z,3])
                
                Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                Col_Markers[2,]<-0#codominant
                Col_Markers[3,]<-0#dropout
                Col_Markers[4,]<-0#error rate
                Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                
                AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec,"_MissDat",MissingVec[M],"_SimNum",s,sep='')
                
                readLines(system.file("extdata/colony2_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                  gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES1",replacement=Num_Canidates,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES2",replacement=Num_Canidates,x=.)%>%
                  gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                  gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                  gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                  gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                  gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                  gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                  gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput.dat",sep=""))
              }
            }#over z markers
          }else if (SexInf==TRUE){
            
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)])}  
              colnames(DatExp)[ncol(DatExp)]<-"Sex"
              if('FRANz'%in% Programs){
                DatExp2<-DatExp[,c(1,2,ncol(DatExp),seq(3,ncol(DatExp)-1,1))]
                write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,".dat",sep="")))
                }
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Sires<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                Col_Dams<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                #running colony from cmd line
                Num_Offspring<-nrow(Col_Offspring)
                Num_Sires<-nrow(Col_Parents[Col_Parents$Sex=="M",])
                Num_Dams<-nrow(Col_Parents[Col_Parents$Sex=="F",])
                Num_Loci<-sum(markers[z,2],markers[z,3])
                
                Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                Col_Markers[2,]<-0#codominant
                Col_Markers[3,]<-0#dropout
                Col_Markers[4,]<-0#error rate
                Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                
                AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec,"_MissDat",MissingVec[M],"_SimNum",s,sep='')
                
                readLines(system.file("extdata/colony2_sex_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE))%>%
                  gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES1",replacement=Num_Sires,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES2",replacement=Num_Dams,x=.)%>%
                  gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                  gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                  gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Sires,x=.)%>%
                  gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Dams,x=.)%>%
                  gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                  gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                  gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput.dat",sep=""))
              }
            }#over z markers
          } else {
            cat("Issue with SexInf")
          }
          
          
        }# over MS miss_data   
      }else{ #if missing data = F then do all the stuff below
        Input2<-InputFull
        MissingVec<-0
        ErrorVec<-0
        
        if(SexInf==F){
          for (z in 1:nrow(markers)){
            mark<-markers[z,-1]
            if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)])}  
            if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)])}  
            if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)])}  
            
            if('FRANz'%in% Programs){
              write.table(DatExp, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
              # Convert to FRANz file using Perl Script
              if(Sys.info()[['sysname']]=="Windows"){  
                shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep=""))) 
              }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep=""))) 
              }
            }# if FRANz
            
            if('Colony' %in% Programs){
              Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
              Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
              Col_Parents<-DatExp[DatExp[,2]==2000,-2]
              Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
              
              #running colony from cmd line
              Num_Offspring<-nrow(Col_Offspring)
              Num_Canidates<-nrow(Col_Parents)
              Num_Loci<-sum(markers[z,2],markers[z,3])
              
              Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
              Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
              Col_Markers[2,]<-0#codominant
              Col_Markers[3,]<-0#dropout
              Col_Markers[4,]<-0#error rate
              Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
              
              AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec,"_MissDat",MissingVec,"_SimNum",s,sep='')
              
              readLines(system.file("extdata/colony2_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                gsub(pattern="NUM_CANDIDATES1",replacement=Num_Canidates,x=.)%>%
                gsub(pattern="NUM_CANDIDATES2",replacement=Num_Canidates,x=.)%>%
                gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Parents2,x=.)%>%
                gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput.dat",sep=""))
            }
          }#over z markers
          }else if (SexInf==TRUE){
            
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)])}  
              colnames(DatExp)[ncol(DatExp)]<-"Sex"
              if('FRANz'%in% Programs){
                DatExp2<-DatExp[,c(1,2,ncol(DatExp),seq(3,ncol(DatExp)-1,1))]
                write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInput","Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep=""))) 
                }else{system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputSim",s,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,".dat",sep=""))) 
                }
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-2]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-2]
                Col_Parents2<-paste(sapply(1:nrow(Col_Parents),function(x) paste(as.character(Col_Parents[x,1]),paste(Col_Parents[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Sires<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                Col_Dams<-Col_Parents%>%
                  filter(.,Sex=="M")%>%
                  select(.,-Sex)%>%
                  {paste(sapply(1:nrow(.),function(x) paste(as.character(.[x,1]),paste(.[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)}
                #running colony from cmd line
                Num_Offspring<-nrow(Col_Offspring)
                Num_Sires<-nrow(Col_Parents[Col_Parents$Sex=="M",])
                Num_Dams<-nrow(Col_Parents[Col_Parents$Sex=="F",])
                Num_Loci<-sum(markers[z,2],markers[z,3])
                
                Col_Markers<-matrix(data=NA,nrow=4,ncol=Num_Loci)
                Col_Markers[1,]<-colnames(DatExp[,seq(3,1+Num_Loci*2,2)])#marker names
                Col_Markers[2,]<-0#codominant
                Col_Markers[3,]<-0#dropout
                Col_Markers[4,]<-0#error rate
                Col_Markers2<-paste(sapply(1:nrow(Col_Markers),function(x) paste(Col_Markers[x,],collapse=", ")),collapse=EOLcharacter)
                
                AnalysisName<-paste('Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"_Error",ErrorVec,"_MissDat",MissingVec,"_SimNum",s,sep='')
                
                readLines(system.file("extdata/colony2_sex_skeleton.dat",package="PseudoBabies",lib.loc=NULL,mustWork=TRUE)) %>%
                  gsub(pattern="NUM_OFFSPRING",replacement=Num_Offspring,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES1",replacement=Num_Sires,x=.)%>%
                  gsub(pattern="NUM_CANDIDATES2",replacement=Num_Dams,x=.)%>%
                  gsub(pattern="NUM_LOCI",replacement=Num_Loci,x=.)%>%
                  gsub(pattern="OFFSPRING_GENOTYPES",replacement=Col_Offspring2,x=.)%>%
                  gsub(pattern="PATERNAL_GENOTYPES",replacement=Col_Sires,x=.)%>%
                  gsub(pattern="MATERNAL_GENOTYPES",replacement=Col_Dams,x=.)%>%
                  gsub(pattern="MARKER_ERROR",replacement=Col_Markers2,x=.)%>%
                  gsub(pattern="DATASETNAME",replacement=AnalysisName,x=.)%>%
                  gsub(pattern="OUTPUTNAME",replacement=AnalysisName,x=.)%>%
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput.dat",sep=""))
              }
            }#over z markers
            } else {
              cat("Issue with SexInf")
            }
            
         
    } # else if there is no missing data  
  
} # else if genotyping error = F
    
}  #over s simulations
  
  if(GPHomogeneity==T){
  #now run homogeneity test on all simulations to see if they are ~ drawn from the same distribution
  GPin2<-lapply(1:length(GPin),function(y) paste(formatC(GPin[[y]][,1],width=max(nchar(GPin[[y]][,1])),flag="-")," , ", unlist(lapply(1:nrow(GPin[[y]]),function(x) paste(paste(formatC(as.numeric(GPin[[y]][x,seq(2,nloci*2+1,2)]),width=3,flag="0"),formatC(as.numeric(GPin[[y]][x,seq(3,nloci*2+1,2)]),width=3,flag="0"),sep=""),collapse=" "))),sep=""))
  
  #lets write out my great looking Genepop file
  OutputFile<-"GP_Simulations.txt"
  write.table("Founders & Simulations as Populations",OutputFile,quote=F,row.names=F,col.names=F)
  write.table(loci,OutputFile,quote=F,row.names=F,col.names=F,append=T)
  for (s in 1:(nSim+1)){
  write.table("POP",OutputFile,quote=F,row.names=F,col.names=F,append=T)
    if(GPHomogeneity==T){write.table(GPin2[[s]],OutputFile,quote=F,row.names=F,col.names=F,append=T)}
  }
  
  genepop::test_diff(OutputFile, genic = TRUE, pairs = FALSE, outputFile = "SimGenicDiff.txt",
            settingsFile = "", dememorization = 10000, batches = 100,
            iterations = 5000, verbose = F)
  }#if genepop
}#function SimulateSingleGenData
