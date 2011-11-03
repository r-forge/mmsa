
## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA,org_name=NA) {

    c(kingdom=kingdom, phylum=phylum, class=class,
	    order=order, family=family, genus=genus, species=species,
	    otu=otu,org_name=org_name)
}



read_Greengenes <- function(object, dir, window=100, overlap=0, last_window=FALSE, word=3)
{
    #green_sequences<- list()
    #num_objects=length(object$data)
    #green_sequences<- list()
    for(f in dir(dir, full.names=T))
    {
	    if(!is(try(sequences <- read.fasta(f)), "try-error")){
	        for(i in 1:length(sequences))
	        {
		        tempObject <-readSequence_Greengenes(object,sequences[i],green_sequences)
            #this is to make sure the data is appended and not overwritten
            
            desc <- tempObject$classification
            object$data[[desc["kingdom"]]][[desc["phylum"]]][[desc["class"]]][[desc["order"]]][[desc["family"]]][[desc["genus"]]][[desc["species"]]][[desc["otu"]]][[desc["org_name"]]] <-tempObject

            
     #       object$data[[num_objects+i]]<-tempObject[[1]]
    #        object$classification[[num_objects+i]]<-tempObject[[2]]
            
            #make tree start
    #        desc <- as.vector(c(object$classification[[i]]))
    #        names(desc) <- as.vector(c(names(object$classification[[i]])))            
            
	        } #for(i in 1:length(sequences))
	    }

    }
    #object$data <- green_sequences
    return(object)
}

#helper function
readSequence_Greengenes <- function(object,currSequence,green_sequences)
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
  for(j in 1:length(classification[[1]]))
  #for(j in 1:8)
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
        cat(paste(classification[[3]][j],"\t"),file="summary.txt",append=T)
        
                
      } #for (j in 1:8)
        cat("\n",file="summary.txt",append=T)
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
        #make a GenClass16S_Greengenes class object
        gen16class<- GenClass16S_Greengenes(kingdom,phylum,class,order,family,genus,species,otu,org_name)
        sequence<- getSequence(currSequence[[1]],as.string=T)
        return(GenSequences(sequence,gen16class, type="sequence"))
              
        
  #end processing
  
}

toNSV <- function(object, window=100, overlap=0, last_window=FALSE, word=3)
{
    num_objects=length(object$data)
    cat("num objects = ",num_objects,"\n")
    green_sequences<- list()
    for(i in 1:length(object$data))
    {
      
      #call the function count_sequences in the file counter.R
      sequence <- object$data[[i]]
      desc <- as.vector(c(object$classification[[i]]))
      names(desc) <- as.vector(c(names(object$classification[[i]])))
      cnt <- count_sequences(sequence,
            window=window, overlap=overlap, word=word, 
            last_window=last_window)
      stream <- make_stream(cnt)
            
      ### save stream at appropriate node in tree                                 
      green_sequences[[desc["kingdom"]]][[desc["phylum"]]][[desc["class"]]][[desc["order"]]][[desc["family"]]][[desc["genus"]]][[desc["species"]]][[desc["otu"]]][[desc["org_name"]]] <- cnt
      #make tree end      
    }
    object$data <- green_sequences
    object$type<- "NSV"
    return(object)

}
