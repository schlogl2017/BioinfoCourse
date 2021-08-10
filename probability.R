## Initializations
p.letter <- vector() ## Initialize a vector to store residue probabilities
p.word <- vector() ## Initialize a vector to store word probabilities
W <- 'GATAAG' ## The word
k <- nchar(W) ## Word length
## Specify residue frequencies
p.letter["A"] <- 0.32
p.letter["C"] <- 0.17
p.letter["G"] <- 0.18
p.letter["T"] <- 0.33
## Check the result
print(p.letter)
sum(p.letter)
## calculate word probability
p.word[W] <- 1 # initialize
for (i in 1:nchar(W)) {
  letter <- substr(W,i,i) ## Select the letter at position i of word W
  p.word[W] <- p.word[W]*p.letter[letter] ## update word probability
}

## Calculate number of positions for a word of length k in a sequence of length L
L <- 1000
pos <- L - k + 1
print(pos)

E.occ <- pos*p.word[W]

# Probabilty to observe exactly 3 occurrences
dbinom(3,size=pos,prob=p.word[W])

## Probabilty to observe at least 3 occurrences
## Thus, in order to calculate the probability to observe at least x, 
## we need to ask the probability to observe 
## more than x-1 successes. In our case, x=3 and x-1=2.
pbinom(2,size=pos,prob=p.word[W],lower.tail=F)

x <- 0:pos ## x is a vector taking all possible positions
y <- dbinom(x,size=pos,prob=p.word[W])
plot(x,y,type="l",col="#0000BB",lwd=2)

# For this, we will select the 11 first values of the occurrence values (x[1:11]), and the corresponding 
# # probability values, which are found in the 11 first entries of the y 
# vector (y[1:11]).
## Plot frequency polygon
plot(x[1:11],y[1:11],type="l",col="#0000BB",lwd=2)
## Plot with histogram type and add a few legends
plot(x[1:10],y[1:10],
     type="h",
     col="#0000BB",
     lwd=2,
     main=paste("Probability of occurrences for word", W),
     xlab="occurrences",
     ylab="probability",
     panel.first=grid(col="#00BB00"))
## Do the same plot, with a logarithmic scale on Y axis.
plot(x[1:10],y[1:10],
     type="h",
     col="#0000BB",
     lwd=2,
     main=paste("Probability of occurrences for word", W),
     xlab="occurrences",
     ylab="probability",
     panel.first=grid(col="#00BB00"),
     log="y")


# Poisson approximation
# Evaluate the effect of using the Poisson distribution as an approximation of the binomial in the tests above.
# For the Poisson distribution, we need only one parameter (the expected mean) instead of two (the number of trials and the probability of success at each trial). Actually, we already calculated the expected mean above : it is the expected number of occurrences for the considered word
x <- 1:10
W <- 'GATAAG'
k <- 6
pos <- L-k+1
## Calculate the P-value of occurrences with a Poisson distribution
E.occ <- pos*p.word[W]
Pval.Poisson <- ppois(x-1, E.occ,lower=F)
Pval.binomial <- pbinom(x-1, pos, p.word[W], lower=F)
## Plot the histogram of the binomial P-value
plot(x[1:10],Pval.binomial,
     type="l",
     col="#0000BB",
     lwd=2,
     main=paste("Probability of occurrences for word", W),
     xlab="occurrences",
     ylab="probability",
     panel.first=grid(col="#00BB00"))
## Superimpose the Poisson approximation
lines(x[1:10],Pval.Poisson[1:10],col="red",lwd=1)
## Plot the histogram of the binomial P-value
plot(x[1:10],Pval.binomial,
     type="l",
     col="#0000BB",
     lwd=2,
     main=paste("Probability of occurrences for word", W),
     xlab="occurrences",
     ylab="probability",
     panel.first=grid(col="#00BB00"),
     log="y")
## Superimpose the Poisson approximation
lines(x[1:10],Pval.Poisson[1:10],col="red",lwd=1)


# In the tutorial on word probabilities, we assumed that the number of oligonucleotide occurrences observed would folow a binomial distribution.
# 
# We will now test this assumption by analyzing the distribution of word occurrences observed in biological sequences. We will fit a binomial distribution on the observed count distribution, and test the goodness of the fitting.
# 
# Tutorial: occurrence probabilities
# The table below shows the distribution of occurrences of the word GATAA in a set of 1000 sequences of 800 base pairs each.

distrib <- data.frame(occ=0:8,
                      obs=c(223, 337, 247, 119, 54, 13, 6, 0, 1)
)
print(distrib)

# Computing relative frequencies
# The relative frequencies are the fraction of sequences containing a given number of occurrences. These are obtained by dividing the number of sequences (obs) by the total number of sequences (N).
# 
# We will first calculate the number of sequences. Actually this number was given above (1000), but we will check that our data contains the right number of ovservations.

N <- sum(distrib$obs) ## Number of sequences
print(N) # sum(c(223,337, 247, 119,  54,  13,   6,   0,   1))

# We can now compute the relative frequencies.
distrib$f <- distrib$obs/N # 223 337 247 119  54  13   6   0   1
print(distrib)

# Fitting a binomial
# We need 2 parameters to characterize a binomial
# the number of trials
# the probability of success at each trial
# Number of trials
# In our case, a trial is a position in a sequence. Each sequence contains L=800bp and the considered word (GATAA) contains k=4 letters, there are n = L - k + 1 possible positions for a GATAA in each sequence (the three last position cannot contain a tetranucleotide).
L <- 800 ## sequence length
k <- 4 ## word length
pos.per.seq <- L - k + 1 ## Number of possible positions for the word in one sequence
print(pos.per.seq)
n <- pos.per.seq

# Probability of succcess (occurrence) at each trial (position)
# How to estimate the probability to find a GATAA at a given position of a sequence ? One possibility would be to estimate that all the tetranucleotides have the same probability, but we have no reason a priori to make such a strong assumption.
# Another possibility would be to estimate the word probability on the basis of the observed distribution itself: we can calculate the average number of occurrences per sequence, and divide this number by the number of positions per sequences to obtain the average number of occurrences per position, i.e. an estimator of the word probability. The total number of occurrences can be obtained by multiplying each value of occurrence (from 0 to 8) by he number of sequences where these occurrences have been observed. This can be done very easily with R, by multiplying the two vectors.
## A vector giving the product of each occurrence 
## value by the corresponding number of sequences
distrib$occ * distrib$obs

## The sum of this vector is the total number of occurrences
occ.total <- sum(distrib$occ*distrib$obs)
print(occ.total)

## Average occurrences per sequence
## multiplicou as
occ.per.seq <- occ.total/N
print(occ.per.seq)
## Average occurrences per position
occ.per.pos <- occ.per.seq/pos.per.seq
print(occ.per.pos) # [1] 0.001898369
# Notice that this number is clearly different from what would have been obtained under a model of equiprobable tetranucleotide (1/4^4=0.0039).

# Computing the theoretical distribution
# We have now the two parameters required for fitting a binomial distribution on our observed distribution.
## Calculate the binomial density for each value between 0 and 8
distrib$dbinom <- dbinom(0:8, size=pos.per.seq, prob=occ.per.pos)

## Calculate expected occurrences for 1000 sequences
distrib$exp <- N*dbinom(0:8, size=pos.per.seq, prob=occ.per.pos)

## Print the distribution and compare the observed and expected distributions
print(distrib)

# Drawing the observed and theoretical distribution
# We will plot the observed distribution as vertical bars, and superimpose the theoretical distribution as a frequency polygon. You can draw a very simple plot with the following command:
plot(distrib$occ,distrib$obs)

# But this is not very aesthetic. The plot can be improved in the following way (the comments can be copy-pasted to R, they will be ignored).
X11(width=8,height=6)

## Plot the observed distribution
plot(distrib$occ,   # X axis values
     distrib$obs,   # Y axis values
     type='h',      # Plot type: histogram-like
     lwd=5,         # Line width
     col="#0000BB", # Line color: blue
     main='GATAA occurrences', # graph title
     xlab='Occurrences', # Label for the X axis
     ylab='Number of sequences', # Label for the Y axis
     ylim=c(0,400),  # Limits of the Y axis
     panel.first=grid(col='#000000') # Grid
)

# Let us now draw the theoretical distribution on top of the observed distribution.
## Draw the binomial distribution over the observed distribution
lines(distrib$occ,
      distrib$exp,
      type='l',
      lwd=2,
      col='#00BB00'
)

# Testing the goodness of the fit
# We will compute the chi2 statistics manually, in order to get familiar with the details of the computation.
# 
# We will first apply the raw formula of the chi2: chi2obs = SUM [(obs -exp) 2 / exp]

distrib$diff <- distrib$obs - distrib$exp
distrib$diff2 <- distrib$diff^2
distrib$chi2 <- distrib$diff^2 / distrib$exp

print(distrib)

chi2.obs <- sum(distrib$chi2)
print(chi2.obs)

# BEWARE: this result is not what we are supposed to obtain. Why ?
#   
#   We applied the raw formula, but we forgot to check two essential conditions for being able to apply Pearson's chi-squared test.
# 
# The weighted sums of observed and expected frequencies should be equal.
# The chi-squared test relies on the assumption that all the expected values are "sufficiently" large. The classical threshold for this condition of applicability is to require exp(i) >= 5 for each class (in our case, for each possible number of occurrences of the word GATAA).
# We will see how to treat these two problems.
# Sums of obervation and expectation must be equal
# The domain of definition of the binomial distribution encompasses alll possible values from 0 to n (the number of trials). In the scripts above, we only computed the binomial density (dbinom) and the corresponding expected values (E(x) = N*dbinom())from 0 to 8. All the values from 9 to n where thus ignored so far.
# 
# A bad way to solve this problem would be to include in our computation of the chi-square all the possible values from 0 to n. This would fix the problem of the sum of expected values, but raise another problem, since at the right tail of the distribution, the values are decreasing very rapidly, so that expected values will be too small to meet the condition of applicability that all expected values are >= 5.
# 
# The correct way to solve this is to merge all the classes on the right tail (from x=9 to x=n) together with the last class where an observation was found (x=8). This ast class will thus from now on represent all the possible values from x=8 to x=n.
# 
# All expected values should be >= 5











