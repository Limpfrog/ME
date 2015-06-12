#' Reformat the file from genwise to specieswise
#' 
#' Function to convert the genewise fasta file of different species to specieswise file of different genes
#' @author Ali Amiryousefii
#' @param dir The path where the fasta files are deposited
#' @param duplicate The method to handle the duplicate genes
#' @export
Gw2Sw <- function(dir=getwd(),duplicate="random"){
  setwd(dir)
  fastaFiles<- unlist(list.files(pattern="*\\.fasta"))#vector of the names of species in .fasta format
  N<- length(fastaFiles)
  for (i in 1:N){
    a<-paste( "awk '{print $1, $2;}'", fastaFiles[i] , sep=" ")
    b<-paste("| sed 's/\\./\ /' | awk '{print $3, $1}'| sed 's/\\[gene=/</' | sed 's/\\]//' | sed 's/>nc_[0-9][0-9][0-9][0-9][0-9][0-9]/>/' | sed 's/\ //' | sed 's/\ //'| sed 's/\ //' | sed 's/>/\ /' | awk '{print $1;}' | sed 's/</>/' >")
    c<- paste(fastaFiles[i], ".genes", sep="")
    d<- paste(b, c, sep="")
    sys<- paste(a, d, sep="")
    system(sys)
    species<- read.fasta(paste(fastaFiles[i], ".genes", sep=""))#fasta file of the ith species
    genes<- names(species)#names of the genes of the ith species
    for (j in 1:length(genes)){#loop for all the genes in the ith species
      gene.add<- species[j]#a gene to be added to the file of genes if already exits or to form a new file if it doesn't
      if (file.exists(paste(genes[j], ".txt", sep=""))){
        gene<- read.fasta(paste(genes[j], ".txt", sep=""))#read the previously saved file of genes
        write.fasta(c(gene, gene.add), names=c(names(gene), sub(".fasta", "",fastaFiles[i] )), file.out=paste(genes[j], ".txt", sep=""))
      }else {
        write.fasta(gene.add, names=sub(".fasta", "",fastaFiles[i]), file.out=paste(genes[j], ".txt", sep=""))	
      }
    }
    
  }
  system("rm *.genes") 
  ##duplicate is the arguments that asignig the handling of the duplicated/triplicated genes 'random' for random choosing
  #for the cases that genes are conserved, 'keep' for the cases that user wants to have the both versions, and 'longer' for discarding the shorter genes 
  geneFiles<- unlist(list.files(pattern="*\\.txt"))#vector of the names of genes in .txt format
  if(duplicate=="random"){
    for (i in 1:length(geneFiles)){
      write.fasta(unique(read.fasta(geneFiles[i])), names=unique(names(read.fasta(geneFiles[i]))), file.out=sub(".txt", ".fasta", geneFiles[i]))
    }
  }
  else if (duplicate=="keep"){
    for (i in 1:length(geneFiles)){
      write.fasta((read.fasta(geneFiles[i])), names=make.unique(names(read.fasta(geneFiles[i]))),file.out=sub(".txt", ".fasta", geneFiles[i]))
    }
  }
  else {
    #This feature to be completed later	
  }
  system("rm *.txt") 
  
} 

#One addendum of the function would be the accomodiation for capital/small letter of the genes