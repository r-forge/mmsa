library(MMSA)

## create/open GenDB
db <- createGenDB("16S.sqlite")
db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))

db

getClassification(db)
nSequences(db)


#new scheme
addSequencesGreengenes(db, system.file("examples/Firmicutes100.fasta",package="MMSA"))

addSequencesGreengenes(db, system.file("examples/Proteobacteria100.fasta",package="MMSA"))

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

