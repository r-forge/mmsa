#create an object of class gencollection
gc<-GenCollection()
#read all fasta files from a directory
gc<- read_Greengenes(gc,"greengenes/")
#count at a particular level
count.GenCollection(gc,"kingdom")
#list at a particular level
list.GenCollection(gc,"class")
