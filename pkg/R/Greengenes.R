## FIXME: make it a tree


## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA) {

    c(kingdom=kingdom, phylum=phylum, class=class,
	    order=order, family=family, genus=genus, species=species,
	    otu=otu)
}



read_Greengenes <- function(object, dir)
{
    for(f in dir(dir, full.names=T))
    {
	if(!is(try(sequences <- read.fasta(f)), "try-error")){
	    for(i in 1:length(sequences))
	    {
		object$data[[i]]<-readSequence_Greengenes(object,sequences[i])        
	    }

	}

    }
    return(object)
}

#helper function
readSequence_Greengenes <- function(object,currSequence)
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

