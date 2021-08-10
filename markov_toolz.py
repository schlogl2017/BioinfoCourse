#!usr/bin/env python
import toolz as tz
from toolz import curried as cur
from glob import glob
import itertools as it
import more_itertools as mit
from functools import partial



def is_sequence(line):
    return not line.startswith('>')

def is_nucleotide(letter):
    # doesn't care abou Ns
    return letter in LDICT1


def get_sequence(path_to_files):
    """Stream a genome, letter by letter, from a list of FASTA filenames."""
    return tz.pipe(path_to_files,
                  cur.map(fasta_reader),
                  tz.concat,
                  cur.filter(is_sequence),
                  # concatenate characters from all lines
                   tz.concat,
                   # discard newlines and 'N'
                   cur.filter(is_nucleotide))


def  genome_gz(file_pattern): 
    """Stream a genome, letter by letter, from a list of FASTA filenames.""" 
    return  tz.pipe(file_pattern,  
                    glob,  
                    sorted,   
                    # Filenames 
                    cur.map(gzopen(mode='rt')),   # lines 
                    # concatenate lines from all files: 
                    tz.concat, # drop header from each sequence 
                    cur.filter(is_sequence), 
                    # concatenate characters from all lines 
                    tz.concat, 
                    # discard newlines and 'N' 
                    cur.filter(is_nucleotide))


LDICT = dict(zip('ACGT', range(8)))
LDICT


# dicitonary of of pairs of letters
PDICT = {(a, b): (LDICT[a], LDICT[b]) for a, b in it.product(LDICT, LDICT)}
PDICT


@tz.curry
def increment_model(model, index):
    """  whenever  we  get  a  new  nucleotide  pair,  say,  (A, T), 
    we want to increment our Markov  model  (our  NumPy  matrix)  at  
    the  corresponding  position. """
    model[index] += 1


def markov(seq):
    """Get a 1st-order Markov model from a sequence of nucleotides."""
    model = np.zeros((4, 4))
    tz.last(tz.pipe(seq,
                    cur.sliding_window(2),
                    # each successive tuple
                    cur.map(PDICT.__getitem__),
                    # location in matrix of tuple
                    cur.map(increment_model(model)))) # increment matrix
    # convert counts to transition probability matrix
    model /= np.sum(model, axis=1)[:, np.newaxis]
    return model


def  plot_model(model,  labels,  figure=None): 
    fig  =  figure  or  plt.figure() 
    ax  =  fig.add_axes([0.1,  0.1,  0.8,  0.8]) 
    im  =  ax.imshow(model,  cmap='magma'); 
    axcolor  =  fig.add_axes([0.91,  0.1,  0.02,  0.8]) 
    plt.colorbar(im,  cax=axcolor) 
    for  axis  in  [ax.xaxis,  ax.yaxis]: 
        axis.set_ticks(range(4)) 
        axis.set_ticks_position('none') 
        axis.set_ticklabels(labels) 
    return  ax


gzopen  =  tz.curry(gzip.open)
model = tz.pipe('Data/Escherichia_coli_k_12.GCA_000800765.1.29.dna.genome.fa.gz', genome_gz, markov)
plot_model(model,  labels='ACGT');








