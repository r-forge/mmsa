library(MMSA)

gc<-GenCollection()
gc

gc<-read_Greengenes(gc,"../greengenes/genus/Clostridium100.fasta")
gc

showRank(gc, "Phy")
nNodesRank(gc, "Phy")

## how many sequences
nSequences(gc)

loc <- findLocation(gc, "Phylum", "Firmicutes" )
loc

## get only the leaves below loc
length(getSequences(gc, loc))

gc.NSV <- toNSV.GenCollection(gc)
gc.NSV

## this is how far I did it for now...

th <- 6*10
m <- genModel.GenCollection(gc.NSV, measure="Manhattan", 
	threshold=th, plus_one=FALSE)

m <- genModel.GenCollection(gc.NSV, loc, measure="Manhattan", 
	threshold=th, plus_one=FALSE)
plot.GenModel(m, method="graph")

m_p <- prune(m, 3, transitions=TRUE, clusters=FALSE)
plot.GenModel(m_p, method="graph")



save(m, file="model.Rda")

