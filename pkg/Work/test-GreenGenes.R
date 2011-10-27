library(MMSA)

#create an object of class gencollection
gc<-GenCollection()
#read all fasta files from a directory
gc<- read_Greengenes(gc,"greengenes/")

## gc has sequences

gc_nsv <- toNSV(gc)

gc_model <- GenModel(gc_nsv)


plot(gc_model)

##subset

gc_nsv_Thermo <- select(gc_nsv, level="class", name="Thermoplasmata")

#count at a particular level
count.GenCollection(gc,"kingdom")
#list at a particular level
list.GenCollection(gc,"class")
