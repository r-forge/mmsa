
db <- createGenDB("phylums.sqlite")
#adds all sequences from a directory to a db
processSequencesGreengenes("pkg/Work/scoreData/phylums", db)
#creates models at a specified level
createModels("pkg/Work/scoreData/models/phylum", rank="phylum", db)
createModels("pkg/Work/scoreData/models/class", rank="class", db)




proteoBacteria<-readRDS("pkg/Work/scoreData/models/phylum/Proteobacteria.rds")
db<-openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
createNSVTable(db,"NSV")
d<-getSequences(db,table="NSV",limit=2)
sequence<-d[[1]]


score(proteoBacteria$model,sequence+1, plus_one=TRUE)
classify("pkg/Work/scoreData/models/phylum/","classify.fasta")

