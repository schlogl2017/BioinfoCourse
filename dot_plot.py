import matplotlib.pyplot as plt
import numpy


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))
    
def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in xrange(m-k+1)] for i in xrange(n-k+1)]  
    
def plotMatrix(M,t, seq1, seq2, nonblank = unichr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)    
    
      
    
#wrapper
def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice
    
    
    
seqx = "ACCTGAGCTCACCTGAGTTA"
seqy = "ACCTGAGCTCACCTGAGTTA"
dotplot(seqx,seqy)










dotplot=plt.imshow(numpy.array(makeMatrix(seqx,seqy,1)))
xt=plt.xticks(numpy.arange(len(list(seqx))),list(seqx))
yt=plt.yticks(numpy.arange(len(list(seqx))),list(seqx))
plt.show()
