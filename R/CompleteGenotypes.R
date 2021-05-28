#' A function to remove missing data
#'
#' This function allows you to remove missing data from your file, by sampling from the allele frequency distributions at each locus
#' @param file .csv Input file name with genotypic information
#' @param StartYear Year that the simulation starts
#' @param nloci Number of loci in the data file
#' @keywords missing data
#' @export
#' @examples
#' CompleteGenotypes()

CompleteGenotypes<-function(file,nloci,StartYear=2008,...){
    #check that the file exists
    if(!file.exists(file)){
        warning("Input file not found.")
    }
  
  
  # read in data - may contain missing genotypes
  Filename1<-file
  dat<-read.csv(Filename1,stringsAsFactors = F)
  
  # Let's first fill in the missing data. If we want missing data we will
  # Basic approach would be to sample from the alleles present to fill in missing genotypes
  
  for (i in (2:(nloci*2+1))){
    if(any(dat[,-1]=='?')){
      stop('Please Code Missing values as 0s not ?s')
    }
    temp<-dat[,i]
    temp[which(temp==0)]<-sample(x=temp[-which(temp==0)],size=length(which(temp==0)),replace=F)
    if(0%in%temp){ #are there any other zeros?
      temp[which(temp==0)]<-sample(x=temp,size=length(which(temp==0)),replace=F)
    }
    dat[,i]<-temp
  }
  
  if(!any(colnames(dat)=="Sex")){dat$Sex<-NA}
  NamesNew<-c(colnames(dat),'BirthYear','DeathYear')
  dat<-cbind(dat,rep(NA,nrow(dat)),rep(2008,nrow(dat)))
  colnames(dat)<-NamesNew
  
  write.table(dat,'FoundersCompleteGenotypes.csv',row.names=F,col.names=T,sep=',')
  cat(paste('No more missing genotypes!!!!', '\n', 
            'Data are located in file FoundersCompleteGenotypes.csv',sep=''))
}
