# We will use R to answer a series of questions on word counts. These analyses will allow us to put in practice several topics treated in the theoretical course. We will follow a problem-driven approach starting from a question, we will see how to solve it with a R script.
# We will start with a tutorial, where the solutions will be shown below each question. A few additional exercises are also proposed.
# 
# In the current tutorial, we will analyze the probabilities associated to word occurrences. In the next tutorial, we will fit various theoretical distributions on an observe distribution of word counts, and compare the significance calculated with these differenc distributions.
# 
# These tutorials assume that you already learned the following chapters of the course.
# 
# Probabilities
# Theoretical distributions
# Fitting
# Goodness of fit
# Test of significance
# Tutorial: occurrence probabilities
# Assuming that a 1000 bp DNA sequence has been generated with a Bernouilli schema, and the following residue probabilities:
#   P(A) = 0.32
# P(T) = 0.33
# P(C) = 0.17
# P(G) = 0.18
# Questions
# What is the probability to observe an occurrence of the word GATAAG at a given position of the sequence ?
#   What is the expected number of occurrences for the word GATAAG in a sequence of this length ?
#   What is the probability to observe exactly 3 occurrences of the word GATAAG ?
#   What is the probability to observe at least 3 occurrences of this word ?
#   Plot the theoretical distribution of probability of occurrences of this word.
# Evaluate the effect of using the Poisson distribution as an approximation of the binomial in the tests above.
# Solutions
# The problem can be formulated as a Bernouilli schema, where each position of the 1000 sequence corresponds to a trial. At each position, the word can either be found (success) or not (failure). We assume that the probability of occurrence (success) is constant along the sequence, and we can thus use the binomial distribution to calculate the probability to observe a certain number of occurrences (successes).
# We need to calculate
# 
# The probability of success at each trial, i.e. the probability to observe the word at a given position of the sequence.
# The number of trials, i.e. the number of possible positions for this word in a 1000 bp sequences.
# Probability of occurrence at a given position
# The probability of the word can be estimated as the product of probability of its letters.
## Initializations
p.letter <- vector() ## Initialize a vector to store residue  probabilities
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

# Number of possible positions
## Calculate number of positions for a word of length k in a sequence of length L
L <- 1000
pos <- L - k + 1
print(pos)

# Expected number of occurrences
# The expected number of occurrences is the product of the number of possible positions by the probability of occurrence at a given position.
E.occ <- pos*p.word[W]

# Occurrence probabilities
# The probability of occurrences can be obtained with the binomial function. Before using it, you should read the on-line help and learn how to use the R implementation of this distribution.

## Probabilty to observe exactly 3 occurrences
dbinom(3,size=pos,prob=p.word[W])

# In order to calculate the probability of observing ##at least 3 matches##, we could sum the results obtained with the same density function, for each valu between 3 and 6. This would however be inefficient in terms if calculation.
sum(dbinom(3:6,size=pos,prob=p.word[W]))

# The distribution function can return the same result in one operation, but we need to be ccareful about the parameters: by *"default, this function returns the lower tail"* of the distribution, i.e. the probability to observe at most x successes.
# We can use the option lower.tail=FALSE to specify that we want the "right tail" of the distribution. ** However, this will return the probability to observe more than x successes**.
# Thus, in order to calculate the probability to observe $$ at least x $$, we need to ask the probability to observe ## more than x-1 successes ##. In our case, x=3 and x-1=2.
# One possibility for calculating the * cumultative sum * would be to calculate each individual probability, with the density function dbinom, and sum them. However, this would be very inefficient, and R provides a * cumulative density function pbinom * to calculate directly the left (default) or the right tail of a distribution.

## Probabilty to observe ** at least 3 occurrences **
pbinom(2,size=pos,prob=p.word[W],lower.tail=FALSE)

# We will now plot the probability to observe exactly x occurrences, as a function of x. for this, we will come back to the density function, since we want the individual probability for each particular occurrence value.

## x is a vector taking all possible positions (N-k+1)
x <- 0:pos
y <- dbinom(x,size=pos,prob=p.word[W])
plot(x,y,type="l",col="#0000BB",lwd=2)

# This graphic is not very informative: we used the whole range of possible values (0:995) for the X axis, but, ** given the small value of the word probability **, the majority of the distribution is restricted to the smaller occurrence values (between 0 and 10). We wil thus plot the distribution over a shorter range.
# For this, we will select the 11 first values of the occurrence values (x[1:11]), and the corresponding probability values, which are found in the 11 first entries of the y vector (y[1:11]).

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
# For the Poisson distribution, we need only one parameter (the expected mean) instead of two (the number of trials and the probability of success at each trial). Actually, we already calculated the expected mean above : it is the expected number of occurrences for the considered word.


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

# Exercise : probability of occurrence with substitutions
# Assuming that a DNA sequence contains the same proportion of A, C, G, T, what is the probability to observe, at a given position of this sequence, the word CAGTGAT
# without substitution ?
#   with exactly one substitution ?
#   with at most one substitution ?
#   with at most 3 substitutions ?
#   with at most 6 substitutions ?
#   Plot the distribution of probabilities, as a function of the number of substitutions.

# Exercices
# A sequence of length 10,000 has the following residue frequencies
# q F(A) = F(T) = 0.325
# q F(C) = F(G) = 0.175
# n What is the probability to observe the word GATAAG at any position of a sequence
# (assuming a Bernouilli model)
# n What would be the probability to observe, in the whole sequence
# q 0 occurrences
# q at least one occurrence
# q exactly one occurrence
# q exactly 15 occurrences
# q at least 15 occurrences
# q less than 15 occurrences
kk = length('GATAAG')
positions = 10000 - kk + 1
p = (0.175)^2*(0.325)^3*0.325
ocb  = positions*p

zero = dbinom(0, positions, p)
print(zero)

uma = dbinom(1, positions, p)
print(uma)

quinze = dbinom(15, positions, p)
print(quinze)
