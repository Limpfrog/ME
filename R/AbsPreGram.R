#' Plot the AbsPreGram
#' 
#' Plotting the matrix of the absense and the presence of the genes in the phylogeny 
#' @param tree the species tree
#' @param dir The directory where all the specieswise files of the genes are available. Note that all the files need to be in the fasta format with the '.fa' extention
#' @author Ali Amiryousefi
#' @export 
plot.apg<- function(tree, dir=getwd()){
fastaFiles<- unlist(list.files(pattern="*\\.fa"))
#genename<-sub(".fa", "", sub("result_", "", fastaFiles))
genename<-sub(".fa", "",fastaFiles)
genename[13]<-"infa" #this is something to remember that the 13 is a bad lock number (for some reasons sub(".f") is not replacing the infA gene!)
species <- sub("'","",sub("'", "", tree$tip.label))
NG<- unlist(list.files(path=paste(getwd(), "/nucleus_genomes", sep=""), pattern="*\\.fasta"))
ng<-sub(".fasta", "",NG)
n<-length(species)
m<-length(fastaFiles)
nrep<- length(NG) #Length of the representative list of nucleus genomes 
nm<-matrix(0, n, m)

for (i in 1:m){
  fasta<-read.fasta(fastaFiles[i])
  for (j in 1:length(names(fasta))){
    index<-which(species==names(fasta)[j])
    nm[index, i]<-1
  }
} #form a nm abspre matrix 1:presence, 0: absence

dist<-cophenetic(tree)
max<-max(dist)


nucl_nearest<-function(spec){
  r<-as.vector(rank(dist[spec,], ties.method="random"))
  q<-1
  while(sum(ng==species[which(r==q)])==0){
    q<-q+1
  }
  species[which(r==q)]
} #find the nearest species that has the nucleus genome


nmup<-nm #name updated for the found blast genes
nmspectrum<-nm #for the alignment similarity
nmd<-nm# For the Phylogenetic distance

system("cp nucleus_genomes/*.fasta ../AbsPre_nucleus")
for (i in 1:nrep){
  sp<-NG[i]
  system(paste(paste("makeblastdb -in", sp, "-dbtype nucl -out ", sep=" "), sp, "_db", sep=""))
} #making the blats databases

for (i in 1:n){
  sp<- species[i]
  sp_spaced<-sub("_", " ", sp)
  for (j in 1:m){
    if (nm[i,j]==0) {
      nearest<-find.nearest(i, j)
      gene<-read.fasta(fastaFiles[j])
      genespecial<-genename[j]
      write.fasta(gene[which(names(gene)==nearest)], file.out=paste(sp, "_", genespecial,"_", nearest, ".cut", sep=""), name=paste(genespecial,"_", nearest, sep=""))
    }
  }
}#cut the nearest gene

for (i in 1:n){
  for (j in 1:m){
    if (nm[i,j]==0) {
      nsp<-nucl_nearest(i)
      sp<- species[i]
      nearest<-find.nearest(i, j)
      genespecial<-genename[j]
      trans1<- paste("blastn -db", paste(nsp, ".fasta_db", sep=""))
      trans2<-paste(trans1, "-query", paste(sp, "_", genespecial, "_", nearest, ".cut", sep=""), "-out", paste(sp, "_", genespecial, "_", nearest, "_", nsp, ".out", sep=""), "-outfmt 6 -max_target_seqs 10 -evalue 0.5")
      system(trans2)
    }
  }
}#Doing the blast


for (i in 1:n){
  sp<- species[i]
  sp_spaced<-sub("_", " ", sp)
  for (j in 1:m){
    if (nm[i,j]==0) {
      nsp<-nucl_nearest(i)
      nearest<-find.nearest(i, j)
      gene<-read.fasta(fastaFiles[j])
      genespecial<-genename[j]
      if (file.info(paste(sp, "_", genespecial, "_", nearest, "_", nsp, ".out", sep=""))$size > 0){#this loop will fill in the nm and nmspectrum
        r<- read.table(paste(sp, "_", genespecial, "_", nearest, "_", nsp, ".out", sep=""))
        nmup[i,j]<-2 #2 for finding the blast hit
        nmspectrum[i,j]<-(r$V3[1]/101) # 101 to avoid absolute 1
        nmd[i,j]<-(1-dist[i, which(species==nsp)]/max) #dividing to scale the date to [0-1]
      }
    }
  }
}#loop for updating the nm, nmspectrum, and nmd based on the blast result

#save all the 4 matrices!

save(nm, file="nm.RData")
save(nmup , file="nmup.RData")
save(nmd, file="nmd.RData")
save(nmspectrum, file="nmspectrum.RData")

#or alternatively save all the material
#save.image("file.RData")

x<-character(n)
for (i in 1:n){x[i]<-nucl_nearest(i)} #a vector of the nearest nucleus for each taxa


NM<-nm


#Plotting the nm for the Alignment found
nm<- nmup
image(1:nrow(nm), 1:ncol(nm),nm, axes=F, col=terrain.colors(6)[c(6,1,2)], main="Abspregram/Hits", sub="", ylab="", xlab="")
axis(1, at=seq(1, n, by=1), labels=sub("_", " ", species), cex.axis=0.36, las=2)
axis(2, at=seq(1, m, by=1), labels=genename, cex.axis=0.36, las=2)
grid(NULL, NULL, lwd = 1)
abline(v=(seq(0,n,5)), col="lightgray", lty="dotted", lwd=0.57)
abline(v=(seq(0,n,1)), col="lightgray", lty="dotted", lwd=0.14)
abline(h=(seq(0,m,5)), col="lightgray", lty="dotted", lwd=0.57)

X<-as.data.frame(nm, row.names=species)
colnames(X)<-gene
cs<-round(colSums(X), 1)
rs<-round(rowSums(X), 1)
axis(4, at=seq(1, m, by=1), labels=cs, cex.axis=0.36, las=2)
axis(3, at=seq(1, n, by=1), labels=rs, cex.axis=0.36, las=2)
}


#' Find the nearest species
#' @param spec The number of the species in the set of the phylogenetic species
#' @param gen The number of the gene in the set of the genes read into the console
#' @return The scientific name of the phylogenetically nearest species that has that gene available
#' @author Ali Amiryousefi
#' @export
find.nearest <- function(spec, gen){#get the species name (number) and gene name(number) and return the closest species that has the gene available 
  #cordinate<-c(which(species==spec), which(genename==gen))
  #a<-cordinate[1]
  #b<-cordinate[2]
  #r<-as.vector(rank(dist[cordinate[1],], ties.method="random"))   
  r<-as.vector(rank(dist[spec,], ties.method="random"))
  q<-1
  while(nm[spec,gen]==0){
    spec<- which(r==q)
    q<-q+1
  }
  species[spec]
}#find the nearest species having the gene, output e.g. "Eucaliptus_epifagus"
