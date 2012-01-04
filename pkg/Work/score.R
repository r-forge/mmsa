
db <- createGenDB("phylums.sqlite")
processSequencesGreengenes("scoreData/phylums", db)
createModels("scoreData/models/phylum", rank="phylum", db)
createModels("scoreData/models/class", rank="class", db)




proteoBacteria<-readRDS("scoreData/models/Proteobacteria.rds")
db<-openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
createNSVTable(db,"NSV")
d<-getSequences(db,table="NSV",limit=2)
sequence<-d[[1]]


score(proteoBacteria$model,sequence+1, plus_one=TRUE)
