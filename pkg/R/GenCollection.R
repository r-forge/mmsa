
## Collection of GenSequences in a classification tree
## Creator function creates an empty collection

#load functions
source("Classification.R")
source("GenSequences.R")
library("seqinr")
#end load functions


GenCollection <- function(classification=GenClass16S(), type=c("NSV", "sequence"), 	annotation=NA) 
{
    type <- match.arg(type)
    #x is the collection object, data contains sequences, 
    x<- list(GenSequences)
    #x <- list(data=list(),classification=list(classification), type= type, annotation=annotation)      
    class(x) <- "GenCollection"
    return(x)
    
}


readfile.GenCollection <- function(object,dir)
{
  #delete all .txt files in  the directory dir
  unlink(paste(dir, "*.txt", sep="/")) 
  for(f in dir(dir, full.names=T))
   {
     if(!is(try(sequences <- read.fasta(f)), "try-error")){
          for(i in 1:length(sequences))
          {
            object[[i]]<-readSequence.GenCollection(object,sequences[i])        
            
          }
       
       }
      
   }
  return(object)
}

#helper function
readSequence.GenCollection <- function(object,currSequence)
{
  seq <- currSequence[[1]]
  #start processing
  annot <- getAnnot(currSequence[[1]])
  # remove leading ">"
  annot <- sub(">", "", annot)
  #split based on ;
  x<-gregexpr(";",annot)
  #start1 for splitting the headers
  classification=list()
  classification[[1]]=c("kingdom","phylum","class","order","family","genus","species","otu")
  classification[[2]]=c("k__","p__","c__","o__","f__","g__","s__","otu_")
  classification[[3]]=""
  #for(j in 1:length(classification[[1]]))
  for(j in 1:8)
  {
    location <-regexpr(text=annot,pattern=classification[[2]][j])
    if(location[1]!=-1)
      {
          y<-which(x[[1]]>location[1])
          if (length(y) == 0)
            y<-nchar(annot)
      		end<- x[[1]][y[1]]
          if (is.na(end))
            end <- nchar(annot)
          if (j<8)
            classification[[3]][j] <- substr(annot,location[1]+3,end-1)
          else if (j==8)
            classification[[3]][j] <- substr(annot,location[1],end-1)
          #remove _
          sub("_","",classification[[3]][j])
          if (is.na(classification[[3]][j]))
          {
            classification[[3]][j]<- 'UNKNOWN'            
          }            
        }
        else
            classification[[3]][j]='UNKNOWN'
                
      } #for (j in 1:8)
        
        org_name<-strsplit(annot,split=" k__")[[1]][1]
        org_name<-sub(";","",org_name)
        org_name <- sub(" ","",org_name)
        kingdom = classification[[3]][1]
        phylum = classification[[3]][2]
        class = classification[[3]][3]
        order = classification[[3]][4]
        family = classification[[3]][5]
        genus = classification[[3]][6]
        species = classification[[3]][7]
        otu = classification[[3]][8]
        #end1
        gen16class<- GenClass16S(kingdom,phylum,class,order,family,genus,species,otu)
        gs <- GenSequences(getSequence(currSequence[[1]],as.string=T),gen16class,"sequence") 
        return(gs) 
        
        
  #end processing
  
}

listLevels.GenCollection <- function(object,level)
{
    
    levels <- c("kingdom","phylum","class","order","family","genus","species")
    m<- match(level,levels)
    if (is.na(m))
    #if (!(any(levels==level)))
      stop("error in level")
    l<-array()    
    for(i in 1:length(object))
    {
      l[i]<-object[[i]]$classification[m]
      
    }
        
    t<-table(l)
    #print(t)
    #l<- unique(l)
    return(t)
}


sequences.GenCollection <- function(object)
{
  for(s in object$data)
    cat(s,"\n")
}

print.GenCollection <- function(object) {
    
    cat("Object of class GenCollection\n")
    
    ## report some basic information
}
    
## adds a GenSequences object to the collection tree at the correct
## position


#length.GenCollection <- function(x) {}

## make count a generic
count.GenCollection <- function(x, level) {
    s
    stop("Not implemented!")
}

#subset.GenCollection <- function() {}

select.GenCollection <- function(x, level, value) {
    stop("Not implemented!")
}


## gets a GenCollection of type sequence and returns the same structure 
## with NSV (counts)
toNSV.GenCollection <- function(x) {
    stop("Not implemented")
}
