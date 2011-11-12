library(MMSA)

gc<-GenCollection()
gc

## find rank as number
findRank(gc, "Phy")

gc<-read_Greengenes(gc,"./greengenes/")
gc

showRank(gc, "Phy")
nNodesRank(gc, "Phy")

## how many sequences
nSequences(gc)

loc <- findLocation(gc, "class", "Deinoco")
loc

## get all squence objects
length(getSequences(gc))
## get only the leaves below loc
length(getSequences(gc, loc))

gc.NSV <- toNSV.GenCollection(gc)
gc.NSV

## this is how far I did it for now...

m <- genModel.GenCollection(gc.NSV, loc, measure="Manhattan", 
	threshold=30, plus_one=FALSE)

plot.GenModel(m, method="graph")

save(m, file="model.Rda")

