library(MMSA)

## create/open GenDB
db <- createGenDB("16S.db")
db <- openGenDB("16S.db")

db

getClassification(db)
nSequences(db)


addSequencesGreengenes(db, "greengenes/phylums/Firmicutes100.fasta")
addSequencesGreengenes(db, "greengenes/phylums/Bacteroidetes100.fasta")


getRank(db)
getRank(db, "gen")
getRank(db, "gen", whereRank="phy", whereName="Bacter")

nSequences(db, "Gen", "Bac")
d<-getSequences(db, "Gen", "Bac")
dim(d)

head(d$sequence)

createNSVTable(db, "NSV")
createNSVTable(db, "NSV2", window=50)
createNSVTable_large(db,"NSV2","phylum","Cyanobacteria") #convert only those sequences to NSV which have phylum=Cyanobacteria

d<-getSequences(db, "Gen", "Bac","NSV")

## FIXME: getSequence should decode automatically
decodeSequence(d$NSV[1])

emm<-genModel(db,"Gen","Bac","NSV")
emm

plot.GenModel(emm)




closeGenDB(db)
dbUnloadDriver(drv)

