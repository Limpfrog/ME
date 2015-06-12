#' Find the nearest species
#' 
#' It will find and do the job
#' @param spec Is the species number
#' @param gen Is the gene of the speceis
#' @author Ali Amiryousefi
#' @export
nucl_nearest<-function(spec){
  r<-as.vector(rank(dist[spec,], ties.method="random"))
  q<-1
  while(sum(ng==species[which(r==q)])==0){
    q<-q+1
  }
  species[which(r==q)]
} #find the nearest species that has the nucleus genome

