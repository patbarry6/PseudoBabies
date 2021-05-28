#' Simulate genotypes over many generations
#'
#' This function allows you to simulate multilocus genotypes through time
#' @param AnalysisName name of the analysis
#' @param Founders .csv file containing the genotypes to use for initiating the simulation. Use function CompleteGenotypes to remove missing data before performing simulations.
#' @param AOMfile .csv of Age of maturity estimates.
#' @param Markerfile .csv for marker panels that will be evaluated.
#' @param NumSNPs number of SNP markers.
#' @param NumuSats number of microsatellite markers.
#' @param NB number of crosses to perform each year. Defaults to 100.
#' @param OffspringDist What distribution should we use for the number of offspring produced for each cross. Options are Uniform (specify a min and a max), NegativeBinomial (specify a prob for the scale parameter see \code{\link[stats]{NegBinomial}} options), and Poisson (specify lambda). 
#' @param SexInfo Does your dataset have information on sex?
#' @param MatingStr Should the simulation allow for polygamy, polygyny, polandry, or monogomy?
#' @param NumYears How many years should the simulation last? Defaults to 50.
#' @param StartYear What year should the simulation start? Defaults to current year.
#' @param nSim=10 How many replicate simulations would you like to do? Defaults to 10.
#' @param Geno_Error Would you like to add genotyping error to the simualated data? Defaults to FALSE
#' @param ErrorVec A vector of error values for genotyping error. Required if Geno_Error=TRUE. Defaults to 0 and 2 percent.
#' @param Miss_data Would you like to add missing data to the simulations? Defaults to FALSE
#' @param MissingVec A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.
#' @param Programs Which programs would you like to create input files for? Defaults to c('FRANz','Colony')
#' @param SexInf Use the sex of individuals in the parentage inference program
#' @param perl.dir File path for the perl directory for FRANz utility for formating. 
#' @param HatchYr In what year do we want to implement a hatchery experiment
#' @param SireHatch How many individuals to collect for the sires
#' @param DamHatch How many individuals to collect for the dams
#' @param HatchDist What distrbution do we want to use for the hatchery offspring?
#' @param BreedDesign What is the breeding design for the experiment? Defaults to 2:1.

#' @importFrom magrittr "%>%"
#' @keywords simulation parentage genetics
#' @export
#' @examples
#' Sim.HatchExp()

Sim.HatchExp<-function(AnalysisName,Founders, AOMfile,Markerfile='MarkerPanels.csv',NumuSats,NumSNPs,NB=600,NumYears=19,StartYear=format(Sys.Date(), "%Y"),nSim=10,
    Geno_Error=FALSE,Miss_data=FALSE,MissingVec=c(0,2.5),ErrorVals=1,Programs=c('FRANz','Colony'), SexInf=FALSE,
    perl.dir, GPHomogeneity=FALSE,OffspringDist="Poisson",SexInfo=FALSE,MatingStr="Polygamy",lambda=2,Umin,Umax,nBprob,
    HatchYear,SireHatch,DamHatch,Hatch_Dist,lambda_h=2,Umin_h,Umax_h,nBprob_h,BreedDesign="1:2"){
  
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
  FoundersFile<-Founders
  dat<-read.csv(FoundersFile,stringsAsFactors = F)
  aom<-read.csv(AOMfile,stringsAsFactors = F)
  markers<-read.csv(Markerfile,header=TRUE, stringsAsFactors = F)
  nloci<-NumSNPs+NumuSats
  loci<-colnames(dat)[seq(from=2,to=nloci*2+1,by=2)]
  if(SexInfo==F){dat$Sex<-sample(c("M","F"),size=nrow(dat),replace = T)}
  
  if(Geno_Error==T){ErrorVec<-seq(1,ErrorVals,1)
  }else{ErrorVec<-0
  }
  if(Miss_data==F){MissingVec<-0}
 
  if(GPHomogeneity==T){
    GPin<-list()
    GPin[[1]]<-dat[,-c(nloci*2+2,nloci*2+3)]
  }
  
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
 
  fixedloci<-rep(list(rep(NA,nloci)),nSim)
  
  PopSize<-list()
  
  
  #### Lets simulate some data!
  for (s in 1:nSim){
    
    #set up where I am going to save all my data
    DataList.names <- seq(StartYear,StartYear+NumYears+max(aom[,1]),1)
    DataList <- vector("list", length(DataList.names))
    names(DataList) <- DataList.names
    
    # Let's put the Founders in a list with year entry being the death year!
    DataList[[1]]<-dat
    
    Count<-paste('(',seq(2,nloci*2,2),':',seq(3,nloci*2+1,2),')',sep='')
    
    #LeslieLike<-matrix(data=NA,ncol=length(which(aom[,2]!=0))+3,nrow=NumYears)
    #colnames(LeslieLike)<-c("PopSize","Crosses","Offspring",aom[which(aom[,2]!=0),1]) 
    
    #Lets create a bunch of genotypes that are unrelated, but pulled from the empirical allele freq distributions
    SimInds<-floor((NB*2)*(1-(cumsum(aom[,2])/100)))
    SimGenos<-lapply(1:length(Count),function(x) as.vector(sample(x=unlist(DataList[[1]][,eval(parse(text=Count[x]))]),size=sum(SimInds)*2,replace=T)))
    SimGenosMat<-as.data.frame(matrix(data=c(paste("Sim",seq(1,sum(SimInds),1),sep="_"),
                               as.numeric(as.character(unlist(SimGenos))),
                               sample(c("M","F"),size=sum(SimInds),replace = T), #information on sex
                               rep(StartYear,sum(SimInds)),#birth year 
                               StartYear+rep(aom[,1],times=SimInds),#death year
                               rep("RandSim",sum(SimInds))),#parents
                               nrow=sum(SimInds),ncol=ncol(DataList[[1]])+1),stringsAsFactors = F)
    colnames(SimGenosMat)<-c(colnames(dat),'Parents')
    
      for (R in 1:nrow(aom)){
      DataList[[aom[R,1]+1]]<-as.data.frame(SimGenosMat[which((as.numeric(SimGenosMat[,ncol(SimGenosMat)-1])-StartYear)==aom[R,1]),],stringsAsfactors=F)
    } 
    
    
    
    for(Y in 1:NumYears){
      Nb<-NB
      #Pull out the year we are using for parents
      parents<-as.data.frame(DataList[[Y]],stringsAsFactors=F)
      
      #pull out broodstock if we are doing the hatchery experiment, or at least make sure we have adequate numbers of each sex
      if(Y==HatchYear){
        HatchSire<-sample(parents[,1],SireHatch,replace = F)
        HatchDam<-sample(parents[,1],DamHatch,replace = F)
        parents_h<-parents[which(parents[,1]%in%c(HatchSire,HatchDam)),]
        parents<-parents[-which(parents[,1]%in%c(HatchSire,HatchDam)),]
        }
      #LeslieLike[Y,1]<-nrow(parents)#save total number of parents in year Y
      
            Nb_sire<-if(nrow(parents[which(parents$Sex=="M"),]) < Nb ) nrow(parents[which(parents$Sex=="M"),]) else Nb
            Nb_dam<-if(nrow(parents[which(parents$Sex=="F"),]) < Nb ) nrow(parents[which(parents$Sex=="F"),]) else Nb
           
                  
          if(MatingStr=="Polygamy"){
            #How many offspring will each set of parents produce?
            if(OffspringDist == "Uniform"){
              OffspringNums<-round(runif(Nb,min=Umin,max=Umax))
            }else if(OffspringDist == "NegativeBinomial") {
              OffspringNums<-rnbinom(n=Nb, size=Nb, prob=nBprob)
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
              OffspringNums<-rnbinom(n=Nb_sire, size=Nb_sire, prob=nBprob)
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
              OffspringNums<-rnbinom(n=Nb_dam, size=Nb_dam, prob=nBprob)
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
              OffspringNums<-rnbinom(n=Nb, size=Nb, prob=nBprob)
            }else if(OffspringDist == "Poisson") {
              OffspringNums<-rpois(n=Nb,lambda=lambda)
            }
            
            Sires<-sample(parents[which(parents$Sex=="M"),1],Nb,replace=F)
            Dams<-sample(parents[which(parents$Sex=="F"),1],Nb,replace=F)
            PAR<-cbind(Sires,Dams)
            }
          
          # so now that we have determined how many progeny each cross will make, lets get rid of all the zeros
          PAR<-PAR[which(OffspringNums!=0),]
          OffspringNums<-OffspringNums[which(OffspringNums!=0)] 
          
          
          CorrectPAR<-vector()
          CorrectPAR<-rep(paste(PAR[,1],'&',PAR[,2],sep=''),OffspringNums)
            
          #if Y = HatchYr then do some hatchery crosses 
          if(Y==HatchYear){
          BD<-unlist(str_split(BreedDesign,":"))
          Nbh<-max(DamHatch,SireHatch)
          if(Hatch_Dist == "Uniform"){
            OffspringNumsH<-round(runif(Nbh,min=Umin_h,max=Umax_h))
          }else if(Hatch_Dist == "NegativeBinomial") {
            OffspringNumsH<-rnbinom(n=Nbh, size=Nbh, prob=nBprob_h)
          }else if(Hatch_Dist == "Poisson") {
            OffspringNumsH<-rpois(n=Nbh,lambda=lambda_h)
          }
          PARh<-cbind(HatchSire,HatchDam)
          PARh<-PARh[which(OffspringNumsH!=0),]
          OffspringNumsH<-OffspringNumsH[which(OffspringNumsH!=0)]
          
          PAR<-rbind(PAR,PARh)
          OffspringNums<-c(OffspringNums,OffspringNumsH)
          
          HatchPAR<-matrix(paste(PARh,"H",sep=""),nrow=nrow(PARh),ncol=2)%>%
            {rep(paste(.[,1],'&',.[,2],sep=''),OffspringNumsH)}
          
          CorrectPAR<-c(CorrectPAR,HatchPAR)
          parents<-rbind(parents,parents_h)
          }
          # end hatchery experiment
          
          PRO<- matrix(data=NA, ncol=(nloci*2), nrow=sum(OffspringNums)) #Progeny
          
        # Fill in offspring genotypes from PAR
        Samp<-lapply(1:nrow(PAR), function(x) eval(parse(text=paste('cbind(',paste(paste('sample(parents[which(parents[,1]%in%PAR[x,]),',Count,'],OffspringNums[x],replace=T)',sep=''),collapse=','),')',sep=''))))
        
        for (i in 1:nrow(PAR)){
          Samp.trans<-t(Samp[[i]])
          PRO[((sum(OffspringNums[1:i])-OffspringNums[i]+1):(sum(OffspringNums[1:i]))),seq(1,nloci*2,2)]<-as.numeric(as.character(Samp.trans[,1]))
          PRO[((sum(OffspringNums[1:i])-OffspringNums[i]+1):(sum(OffspringNums[1:i]))),seq(2,nloci*2,2)]<-as.numeric(as.character(Samp.trans[,2]))
        }
        
        
        
        # When will each reproduce and die
        # lets sample from the age of maturity scheduale
        Ages<-sample(x=aom[,1], size=nrow(PRO), replace=TRUE, prob=aom[,2]/sum(aom[,2]))
        AgeDeath<-DataList.names[Y]+Ages
        
        
        # LeslieLike[Y,2]<-nrow(PAR)#the number of crosses occuring
        # LeslieLike[Y,3]<-sum(OffspringNums) #the number of offspring produced
        # LeslieLike[Y,which(colnames(LeslieLike)%in%names(table(Ages)))]<-t(table(Ages)) #the number of offspring produced
        
        #make some new snazzy names for each so we can keep track of when they were born
        Names<-paste(DataList.names[Y],'_',seq(1,nrow(PRO),1),sep='')#This used to be paste('BY',DataList.names[Y],'_',seq(1,nrow(PRO),1),sep='')
        
        NewDat<-as.data.frame(cbind(Names,PRO,sample(x=c("M","F"),nrow(PRO),replace=T),rep(DataList.names[Y],nrow(PRO)),AgeDeath,CorrectPAR),stringsAsFactors = F)
        colnames(NewDat)<-c(colnames(dat),'Parents')
        
        # now I need to store each individual in their own year according to their death year
        for (R in 1:nrow(aom)){
          DataList[[aom[R,1]+Y]]<-rbind(DataList[[aom[R,1]+Y]],NewDat[NewDat[,(ncol(NewDat)-1)]==DataList.names[aom[R,1]+Y],])
        } 
      
    cat(Y)
      }#over Y years 
    
    PopSize[[s]]<-sapply(1:NumYears,function(x) nrow(DataList[[x]]))
    
    # Now to get all the data we need for the analysis
    #what year are we going to assign to parents?
    Ass.year<-sapply(1:4,function(x) StartYear+NumYears-x) #use the last year of the simulated genotypes
    #what are the years that could be potential parents for Ass.year
    Ages<-aom[aom[,2]>0,1]
    
    for (a in 1:4){
    
    Input<-eval(parse(text=paste('rbind(DataList[[(Ass.year[',a,']-StartYear+1)]],',paste(paste('DataList[[',(Ass.year[a]-StartYear+1-Ages),
                                                                                                        ']]',sep=''),collapse=','),')',sep='')))
    
    
    #now to move last two columns to the front of the input file
    InputFull<-cbind(Input[,1],c(rep(2001,nrow(DataList[[(Ass.year[a]-StartYear+1)]])),
                                 rep(2000,(sum(sapply(Ass.year[a]-StartYear+1-Ages,function(x) nrow(DataList[[x]])))))),
                     Input[,seq(2,nloci*2+4)])
    
    
    if(GPHomogeneity==T) GPin[[s+1]]<-InputFull[-which(InputFull[,2]==2000),-2]
    
    colnames(InputFull)<-c(colnames(dat)[1],'Year',colnames(dat)[c(-1,-(length(colnames(dat))-2):-(length(colnames(dat))))])
    
    
    TRUTH<-DataList[[Ass.year[a]-StartYear+1]][,c(1,ncol(DataList[[Ass.year[a]-StartYear+1]]))]
  
    # Write out each analysis  
    
    write.table(TRUTH, file = paste(wd,"/SimParents/","Truth",a,"Sim", s, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
    
    
#### Geno Error + Missing Data ####    
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
                if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)])}  
                if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)])}  
                if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input4)])}  
                colnames(DatExp)[ncol(DatExp)]<-"Pop"
                
                if('FRANz'%in% Programs){
                  for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                    DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,seq(3,ncol(DatExp)-1,1))]
                    write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                   # Convert to FRANz file using Perl Script
                    if(Sys.info()[['sysname']]=="Windows"){  
                      shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                    }else{
                      system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                    }
                  }# for p in pops
                
                DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                
                
                GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                diag(GeoMat)<-"0.000"
                GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                
                writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""))
                write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""),append=T)
                
                files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,sep="")),value=T)
                file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/",files2delete,sep=""))
              
                }#if FRANz
                
                if('Colony' %in% Programs){
                  Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                  Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                    writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
                }
              }#over z markers
            } else if (SexInf==TRUE){
              
              for (z in 1:nrow(markers)){
                mark<-markers[z,-1]
                if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)-2],Input4[,ncol(Input4)])}  
                if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)-2],Input4[,ncol(Input4)])}  
                if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input4)-2],Input4[,ncol(Input4)])}  
                colnames(DatExp)[ncol(DatExp)-1]<-"Sex"
                colnames(DatExp)[ncol(DatExp)]<-"Pop"
                
                if('FRANz'%in% Programs){
                  for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                    DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,ncol(DatExp)-1,seq(3,ncol(DatExp)-2,1))]
                    
                    write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                    # Convert to FRANz file using Perl Script
                    if(Sys.info()[['sysname']]=="Windows"){  
                      shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                    }else{
                      system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                    }
                  }#for p in pops
                  DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                  writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                  write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                  
                  
                  GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                  diag(GeoMat)<-"0.000"
                  GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                  GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                  
                  writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""))
                  write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""),append=T)
                  
                  files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,sep="")),value=T)
                  file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/",files2delete,sep=""))
                
                }# if FRANz
                
                if('Colony' %in% Programs){
                  Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                  Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                  Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                    writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
                }#if Colony is in program files
              }#over z markers  
            } #if SexInf==T
          }#over MS Miss_data  
#### GenoErrro + No Missing Data ####          
        }else{ #if missing data = F
          Input4<-InputErr
          MissingVec<-0
          if(SexInf==F){
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input4)])}  
              colnames(DatExp)[ncol(DatExp)]<-"Pop"
              
              if('FRANz'%in% Programs){
                for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                  DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,seq(3,ncol(DatExp)-1,1))]
                  write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }else{
                    system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }
                }# for p in pops
                
                DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                
                
                GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                diag(GeoMat)<-"0.000"
                GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                
                writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""))
                write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""),append=T)
                files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,sep="")),value=T)
                file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/",files2delete,sep=""))
              }#if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
              }
            }#over z markers
          }else if (SexInf==TRUE){
            
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input4[,ncol(Input4)-2],Input4[,ncol(Input4)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input4[,1:2],Input4[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input4[,ncol(Input4)-2],Input4[,ncol(Input4)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input4[,1:2],Input4[,3:as.numeric(mark[1]*2+2)],Input4[,ncol(Input2)-2],Input4[,ncol(Input4)])}  
              colnames(DatExp)[ncol(DatExp)-1]<-"Sex"
              colnames(DatExp)[ncol(DatExp)]<-"Pop"
              
              if('FRANz'%in% Programs){
                for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                  DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,ncol(DatExp)-1,seq(3,ncol(DatExp)-2,1))]
                  
                  write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }else{
                    system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }
                }#for p in pops
                DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                
                
                GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                diag(GeoMat)<-"0.000"
                GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                
                writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""))
                write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""),append=T)
                
                files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,sep="")),value=T)
                file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/",files2delete,sep=""))
              }# if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec[e],"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
              }
            }#over z markers
          } else {
            cat("Issue with SexInf")
          }
          
          
        } # else 
      }# over ErrorVec
    }else{ #if Geno_Error is false Below this point we should not index on e
      ErrorVec<-0
#### No Errro + Missing Data ####      
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
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)])}  
              colnames(DatExp)[ncol(DatExp)]<-"Pop"
              
              if('FRANz'%in% Programs){
                for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                  DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,seq(3,ncol(DatExp)-1,1))]
                  write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }else{
                    system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }
                }# for p in pops
                
                DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                
                
                GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                diag(GeoMat)<-"0.000"
                GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                
                writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""))
                write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""),append=T)
                files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,sep="")),value=T)
                file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/",files2delete,sep=""))
              }#if FRANz
              
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
              }
            }#over z markers
          }else if (SexInf==TRUE){
            
            for (z in 1:nrow(markers)){
              mark<-markers[z,-1]
              if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
              if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
              if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
              colnames(DatExp)[ncol(DatExp)-1]<-"Sex"
              colnames(DatExp)[ncol(DatExp)]<-"Pop"
              
              if('FRANz'%in% Programs){
                for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                  DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,ncol(DatExp)-1,seq(3,ncol(DatExp)-2,1))]
                  
                  write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                  # Convert to FRANz file using Perl Script
                  if(Sys.info()[['sysname']]=="Windows"){  
                    shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }else{
                    system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                  }
                }#for p in pops
                DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
                writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
                write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
                
                
                GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
                diag(GeoMat)<-"0.000"
                GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
                GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
                
                writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""))
                write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/Geo.dat",sep=""),append=T)
                
                files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,sep="")),value=T)
                file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/",files2delete,sep=""))
              }# if FRANz
                
              if('Colony' %in% Programs){
                Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
                Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
                Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                  writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec[M],"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
              }
            }#over z markers
          }else{
            cat("Issue with SexInf")
          }
          
          
        }# over MS miss_data   
      }else{ #if missing data = F then do all the stuff below
#### No Errro + No Missing Data ####        
        Input2<-InputFull
        MissingVec<-0
        ErrorVec<-0
        
        if(SexInf==F){
          for (z in 1:nrow(markers)){
            mark<-markers[z,-1]
            if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)])}  
            if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)])}  
            if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)])}  
            colnames(DatExp)[ncol(DatExp)]<-"Pop"
            
            if('FRANz'%in% Programs){
              for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
                DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,seq(3,ncol(DatExp)-1,1))]
                write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
                # Convert to FRANz file using Perl Script
                if(Sys.info()[['sysname']]=="Windows"){  
                  shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                }else{
                  system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 2 --birth_col 1 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
                }
              }# for p in pops
              
              DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
              writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
              write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
              
              
              GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
              diag(GeoMat)<-"0.000"
              GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
              GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
              
              writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""))
              write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""),append=T)
              files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,sep="")),value=T)
              file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/",files2delete,sep=""))
            }#if FRANz
          
            
            if('Colony' %in% Programs){
              Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
              Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
              Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
            }
          }#over z markers
        }else if (SexInf==TRUE){
          
          for (z in 1:nrow(markers)){
            mark<-markers[z,-1]
            if(mark[1]>0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,(2+NumuSats*2+1):(as.numeric(mark[2])*2+(NumuSats*2)+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
            if(mark[1]==0&mark[2]>0){DatExp<-cbind(Input2[,1:2],Input2[,(2+NumuSats*2+1):as.numeric(mark[2]*2+NumuSats*2+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
            if(mark[1]>0&mark[2]==0){DatExp<-cbind(Input2[,1:2],Input2[,3:as.numeric(mark[1]*2+2)],Input2[,ncol(Input2)-2],Input2[,ncol(Input2)])}  
            colnames(DatExp)[ncol(DatExp)-1]<-"Sex"
            colnames(DatExp)[ncol(DatExp)]<-"Pop"
            
            if('FRANz'%in% Programs){
              for (p in 1:length(unique(DatExp[,ncol(DatExp)]))){
              DatExp2<-DatExp[which(DatExp[,ncol(DatExp)]==unique(DatExp[,ncol(DatExp)])[p]),c(1,2,ncol(DatExp)-1,seq(3,ncol(DatExp)-2,1))]
              
              write.table(DatExp2, file = paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim", s,"_",p, ".csv",sep=""),sep = ",",row.names = FALSE, quote = FALSE)
              # Convert to FRANz file using Perl Script
              if(Sys.info()[['sysname']]=="Windows"){  
                shell(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
              }else{
                system(paste("perl ",perl.dir, " --in", paste (wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/","FRANzInputTemp",a,"Sim",s,"_",p,".csv",sep=""), "--data_col 3 --birth_col 1 --sex_col 2 --has_header --missing_allele 0", paste(">",wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",p,".dat",sep=""))) 
              }
              }#for p in pops
              DatLines<-unlist(lapply(1:length(unique(DatExp[,ncol(DatExp)])),function(x) readLines(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSimTemp",s,"_",a,"_",x,".dat",sep=""))[-1]))
              writeLines(paste(length(unique(DatExp[,ncol(DatExp)])),sum(mark),".","csv.pl",sep=" "),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""))
              write_lines(x=DatLines,path=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/FRANzSim",s,"_",a,".dat",sep=""),append=T)
              
              
              GeoMat<-matrix(data="2.000",nrow=length(Ages)+1,ncol=length(Ages)+1)
              diag(GeoMat)<-"0.000"
              GeoMat<-cbind(formatC(c(0,Ages),width=10),GeoMat)
              GeoFile<-sapply(1:(length(Ages)+1),function(x) paste(GeoMat[x,],collapse=" "))
              
              writeLines(as.character(length(Ages)+1),paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""))
              write_lines(GeoFile,paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/Geo.dat",sep=""),append=T)
              
              files2delete<-grep("Temp",list.files(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,sep="")),value=T)
              file.remove(paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/",files2delete,sep=""))
              }# if FRANz
            
            if('Colony' %in% Programs){
              Col_Offspring<-DatExp[DatExp[,2]==2001,-c(2,ncol(DatExp))]
              Col_Offspring2<-paste(sapply(1:nrow(Col_Offspring),function(x) paste(as.character(Col_Offspring[x,1]),paste(Col_Offspring[x,-1],collapse=" "),sep=" ")),collapse=EOLcharacter)
              Col_Parents<-DatExp[DatExp[,2]==2000,-c(2,ncol(DatExp))]
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
                writeLines(text=.,con=paste(wd,"/",'Panel_uSats',markers[z,2],'_SNPs',markers[z,3],"/Error_",ErrorVec,"/MissDat_",MissingVec,"/SimNum",s,"/ColonyInput",a,".dat",sep=""))
            }
          }#over z markers
        } else {
          cat("Issue with SexInf")
        }
        
        
      } # else if there is no missing data  
      
    } # else if genotyping error = F
    }#over a analyses
  }  #over s simulations

  if(GPHomogeneity==T){
  #now run homogeneity test on all simulations to see if they are ~ drawn from the same distribution
  GPin2<-lapply(1:length(GPin),function(y) paste(formatC(GPin[[y]][,1],width=max(nchar(as.character(GPin[[y]][,1]))),flag="-")," , ", unlist(lapply(1:nrow(GPin[[y]]),function(x) paste(paste(formatC(as.numeric(GPin[[y]][x,seq(2,nloci*2+1,2)]),width=3,flag="0"),formatC(as.numeric(GPin[[y]][x,seq(3,nloci*2+1,2)]),width=3,flag="0"),sep=""),collapse=" "))),sep=""))
  
  #lets write out my great looking Genepop file
  OutputFile<-"GP_Simulations.txt"
  write.table("Founders & Simulations as Populations",OutputFile,quote=F,row.names=F,col.names=F)
  write.table(loci,OutputFile,quote=F,row.names=F,col.names=F,append=T)
  for (s in 1:(nSim+1)){
    write.table("POP",OutputFile,quote=F,row.names=F,col.names=F,append=T)
    write.table(GPin2[[s]],OutputFile,quote=F,row.names=F,col.names=F,append=T)
  }
  
  genepop::test_diff(OutputFile, genic = TRUE, pairs = TRUE, outputFile = "SimGenicDiff.txt",
                     settingsFile = "", dememorization = 10000, batches = 100,
                     iterations = 5000, verbose = F)
  }#if running genepop
  
 PopSizeSims<-matrix(data=NA,nrow=nSim*NumYears,ncol=3)
 PopSizeSims[,1]<-rep(paste("Sim",1:nSim,sep=""),each=NumYears)
 PopSizeSims[,2]<-as.numeric(as.character(rep(1:NumYears,times=nSim)))
 PopSizeSims[,3]<-as.numeric(as.character(unlist(PopSize)))
 colnames(PopSizeSims)<-c("Simulation","Year","PopSize")
 PopSizeSims<-as.data.frame(PopSizeSims,stringsAsFactors=F)
 
 png("PopulationSize_Simulations.jpg", width = 600, height = 400)
 p<-ggplot2::ggplot(PopSizeSims, ggplot2::aes(x=as.numeric(Year),y=as.numeric(as.character(PopSize)),color=Simulation)) + 
 ggplot2::geom_line()+
   ggplot2::xlab("Year")+
   ggplot2::ylab("Population Size")
 print(p)
 dev.off()
  
  
   }#function SimulateData




