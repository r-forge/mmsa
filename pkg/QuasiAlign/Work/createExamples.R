#	dirpath is the path of the directory that contains split up phylum files
#
createExamples<-function(dirpath,outpath)
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
	
	db<-createGenDB("16S.sqlite")
	processSequences(dir=outpath,db)
	createModels(modelDir=file.path(outpath,"models"),rank="Phylum",db)

}
