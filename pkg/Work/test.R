# library("mmsa")
library("seqinr")
library("rEMM")



path <- "../R"
for(f in dir(path, full.names=TRUE)) source(f)

## custom classification hierarchy for NCBI
GenClass16S_NCBI <- function(phylum=NA, class=NA, order=NA,
	family=NA, genus=NA, species=NA, strain=NA) {

    c(phylum=phylum, class=class, species=species,
	    strain=strain)
}


sequences <- read.fasta("Zymomonas_mobilis_ZM416s.wri")
annot <- getAnnot(sequences)
annot <- sub(">", "", annot)
annot <- strsplit(annot, split=": ")
annot

## create a GenSequence
s <- lapply(getSequence(sequences), c2s)
gc <- GenClass16S_NCBI(
	strain="Mobilis ZM4",
	species="Zymomonas mobilis", 
	class="Alphaproteobacteria", 
	phylum="Proteobacteria" 
	)

gs <- GenSequences(s, class=gc, type="seq")

## tool chain
NSVs <- toNSV(gs)

model <- GenModel(NSVs)


## plot emm
plot(model)

## show sequence names (data.frame with the classification)
model$sequences

## plot 1st sequence red and 2nd green
plot(model, mark=c(1, 2), col=c("red", "green")) 


## implement tool chain for GenCollection
## selection for GenCollection
