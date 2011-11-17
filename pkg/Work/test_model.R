library(MMSA)

drv<-dbDriver("SQLite");

## create/open DB
db <- openGenDB("16S.db")
## create/open collection
sequ <- createGenCollection(db, "sequences")
#sequ <- openGenCollection(db, "sequences") 

addSequencesGreengenes_large(sequ, "greengenes/current_GREENGENES_gg16S_unaligned_300.fasta")

nSequences(sequ)


nsv <- createGenCollection(db, "nsv")

toNSV(sequ, nsv)

emm<-genModel(nsv,"Gen","Bac")

plot.GenCollection(emm, method="graph")


closeGenDB(db)
dbUnloadDriver(drv)

