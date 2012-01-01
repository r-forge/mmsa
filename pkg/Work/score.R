proteoBacteria<-.readRDS(system.file("models/Proteobacteria.rds",package="MMSA"))
db<-openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
createNSVTable(db,"NSV")
d<-getSequences(db,table="NSV",limit=2)
sequence<-as.data.frame(d[1])
score(proteoBacteria,sequence) #always returns 0
