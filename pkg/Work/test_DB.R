library(MMSA)

drv<-dbDriver("SQLite");


## create/open DB
db <- openGenDB("16S.db")

## create/open collection
sequ <- createGenCollection(db, "sequences")

#sequ <- openGenCollection(db, "sequences") 

addSequencesGreengenes(sequ, "greengenes/phylums/Firmicutes100.fasta")
addSequencesGreengenes(sequ, "greengenes/phylums/Bacteroidetes100.fasta")

nSequences(sequ)

getRank(sequ)
getRank(sequ, "gen")
getRank(sequ, "gen", whereRank="phy", whereName="Bacter")

nSequences(sequ, "Gen", "Bac")
d<-getSequences(sequ, "Gen", "Bac")
dim(d)

head(d$sequence)

#nsv <- createGenCollection(db, "nsv")

toNSV(sequ, "NSV")
d<-getSequences(sequ, "Gen", "Bac","NSV")
decodeSequence(d$NSV[1])
length(d$NSV)

emm<-genModel(sequ,"Gen","Bac","NSV")

plot.GenModel(emm)




closeGenDB(db)
dbUnloadDriver(drv)

