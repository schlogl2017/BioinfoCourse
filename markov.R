# sample() is employed to generate an element to be placed at position i + 1 given a particular letter at i.
# Note particularly how the probability distributions pi and P[i,], the rows of the transition matrix, are employed for each use of sample(). Each application of markov1 will produce a different string (check this), but the overall properties of each string should be similar.

markov1 <- function(x,pi,P,n){
  # x = vector [1 2 3 4] representing A, C, G, T, respectively
  # pi = the probability distribution for X0: (1x4 row vector)
  # P = transition matrix (4x4)
  # n = length of simulated sequence
  # Initialize vector to contain simulated sequence
  mg <- rep(0,n)
  # Produce initial element
  mg[1] <- sample(x,1,replace=TRUE,pi)
  for(k in 1:(n-1)){
    mg[k+1]<-sample(x,1,replace=T,P[mg[k],])
  }
  return(mg)
}

# To input the parameters in the simulation, we use:
x <- c(1:4)
pi <- c(.342,.158,.158,.342)
P <- matrix(scan(),ncol=4, nrow=4,byrow=T)
tmp<-markov1(x,pi,P,50000)

# We can check the base composition (remembering that C is represented by 2 and G is represented by 3) of the generated sequence:
length(tmp[tmp[]==1])
length(tmp[tmp[]==2])
length(tmp[tmp[]==4])
(8000+7905)/(16697+8000+7905+17398) # compute fr(G+C)

# This compares favorably with the value 31.6% G+C given in the transition matrix. Now we check whether tmp contains an appropriate fraction of CG dinucleotides:
count=0
# 49999 because i+1 undefined for 50000
for(i in 1:49999){ 
  if(tmp[i]==2 && tmp[i+1]==3)
    count<- count+1
  }
count
count/49999
# From (2.27), the relative abundance of the CG dinucleotide in M. genitalium is 0.010, whereas the string produced by the Markov model contains CG at a
# relative abundance 0.0096. This matches the observed data well.




