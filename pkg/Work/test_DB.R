library(MMSA)

## create/open GenDB
db <- createGenDB("16S.sqlite")
db <- openGenDB("16S.sqlite")

db

getClassification(db)
nSequences(db)


addSequencesGreengenes(db, system.file("data/Firmicutes100.fasta",package="MMSA"))
#new scheme
addSequencesGreengenes(db, system.file("examples/Firmicutes100.fasta",package="MMSA"))

getRank(db)
getRank(db, "gen")
getRank(db, "gen", whereRank="phy", whereName="Firm")


nSequences(db, "Gen", "Des")
d<-getSequences(db, "Gen", "Bac")
dim(d)

head(d$sequence)

createNSVTable(db, "NSV")
createNSVTable(db, "NSV2", window=50)
createNSVTable(db,"NSV2","phylum","Cyanobacteria") #convert only those sequences to NSV which have phylum=Cyanobacteria

d<-getSequences(db, "Gen", "Des","NSV")



emm<-genModel(db,"Gen","Des","NSV")
emm

plot.GenModel(emm)




closeGenDB(db)
dbUnloadDriver(drv)

