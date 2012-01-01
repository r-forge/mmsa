
## classification hierarchy for 16S
GenClass16S <- function(domain=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, strain=NA) {

    c(domain=domain, phylum=phylum, class=class,
	    order=order, family=family, genus=genus, species=species,
	    strain=strain)
}


.createModels <- function(dir=system.file("phylums",package="MMSA"))
{
	db<-createGenDB("phylums.sqlite")
		
	for(f in dir(dir, full.names=T))
	{
		tmpName<-as.character(f)
		tmpName<-unlist(strsplit(tmpName,"/"))
		phylumName<-tmpName[length(tmpName)]
		phylumName<-sub(".fasta","",phylumName)	
		phylumName<-sub("100","",phylumName)	
		print(phylumName)
		cat("Creating model for ",phylumName,"\n")	
		addSequencesGreengenes(db, f)
		print("starting createNSV")
		createNSVTable(db, "NSVPhylum","phylum",phylumName)
		emm<-genModel(db,table="NSVPhylum")
		outfileName<-paste(phylumName,".rds",sep="")
		.saveRDS(emm,file=outfileName)
		dropNSVTable(db,"NSVPhylum")			
	}
	rm(db)
	unlink("phylums.sqlite")
	
}

classify<-function()
{

}
