library("seqinr")
library("rEMM")

dir <- "greengenes"
source("GenCollection.R")
#call the constructor for GenCollection object which is an empty collection of sequences
gc <- GenCollection()
#read file and put into collection of sequences      
gc<-readfile.GenCollection(gc,"greengenes/")      
cat("Size of gc object is ",length(gc),"\n")
listLevels.GenCollection(gc,"kingdom")
listLevels.GenCollection(gc,"phylum")
listLevels.GenCollection(gc,"class")
listLevels.GenCollection(gc,"order")
listLevels.GenCollection(gc,"family")
#should produce error
listLevels.GenCollection(gc,"abcd")