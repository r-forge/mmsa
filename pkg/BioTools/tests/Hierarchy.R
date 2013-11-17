library(BioTools)

db <- openGenDB(system.file("examples/16S.sqlite",package="BioTools"))

getTaxonomyNames(db)

getRank(db, rank="genus")
getRank(db, rank="genus", count=TRUE)
getRank(db, rank="species", whereRank="genus", whereName="Desulfotomaculum")
getRank(db, rank="species", whereRank="genus", whereName="D%")
getRank(db, rank="species", whereRank="genus", 
	whereName=c("Desulfotomaculum", "Syntrophomonas"))
getRank(db, rank="S", whereRank="G", 
	whereName=c("De%", "Sy%"))
getRank(db, rank="species", whereRank="genus", whereName="Desulfotomaculum",
	all=TRUE)
getRank(db, rank="species", whereRank="genus", whereName="Doesnotexist")

getIDs(db)
ids <- getIDs(db, whereRank="genus", whereName="Desulfotomaculum")
ids


getHierarchy(db, rank="genus", name="Syntrophomonas")
getHierarchy(db, rank="genus", name="S%")
getHierarchy(db, rank="genus", name=c("Se%", "Sy%"))
getHierarchy(db, rank="genus", name="doesnotexist!")
getHierarchy(db, rank="id", name="13655")
getHierarchy(db, rank="id", name=ids)

#getHierarchy(db)

closeGenDB(db)

