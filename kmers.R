library(stringi)
library(Biostrings)

DNAlst<-list("CAAACTGATTTT","GATGAAAGTAAAATACCG","ATTATGC","TGGA","CGCGCATCAA")
dna <- stri_paste(rep(c("A","C","G","T"),each=4),c("A","C","G","T"))
result <- t(sapply(DNAlst, stri_count_fixed,pattern=dna,overlap=TRUE))
colnames(result) <- dna
result

t(sapply(DNAlst, function(x){x1 <-  DNAString(x)
oligonucleotideFrequency(x1,2)}))


sequence_kmers <- function(sequence, k){
  k_mers <- lapply(sequence,function(x){
    seq_loop_size <- length(DNAString(x))-k+1
    
    kmers <- sapply(1:seq_loop_size, function(z){
      y <- z + k -1
      kmer <- substr(x=x, start=z, stop=y)
      return(kmer)
    })
    return(kmers)
  })
  
  uniq <- unique(unlist(k_mers))
  ind <- t(sapply(k_mers, function(x){
    tabulate(match(x, uniq), length(uniq))
  }))
  colnames(ind) <- uniq
  
  return(ind)
}


s <- unlist(DNAlst[1])
sequence_kmers(DNAlst[1],2)
