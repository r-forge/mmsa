#	dirpath is the path of the directory that contains split up phylum files
#
createExampleSeq<-function(dirpath, outpath="./inst/examples")
{
	proteobacteria <- file.path(dirpath,"Proteobacteria.fasta")
	planctomycetes <- file.path(dirpath,"Planctomycetes.fasta")
	firmicutes  <- file.path(dirpath,"Firmicutes.fasta")

	sample1<-sample(read.DNAStringSet(proteobacteria),100)
	sample2<-sample(read.DNAStringSet(planctomycetes),100)
	sample3<-sample(read.DNAStringSet(firmicutes),100)

	write.XStringSet(sample1, filepath=file.path(outpath,"phylums/Proteobacteria.fasta"))
	write.XStringSet(sample2, filepath=file.path(outpath,"phylums/Planctomycetes.fasta"))
	write.XStringSet(sample3, filepath=file.path(outpath,"phylums/Firmicutes.fasta"))
}
 
createExampleDb <- function(outpath="./inst/examples"){
	unlink(file.path(outpath, "16S.sqlite"))
  db<-createGenDB(file.path(outpath, "16S.sqlite"))
	processSequences(dir=file.path(outpath, "phylums"), db = db)
}

