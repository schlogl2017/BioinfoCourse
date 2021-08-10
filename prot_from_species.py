import sys
import argparse
from Bio import Entrez

parser = argparse.ArgumentParser(description='Searches for protein sequences in the Title Word field ([TITL]) based on any provided key terms.\nSee here for further details: http://cbsu.tc.cornell.edu/resources/seq_comp/pb607_introductory/entrez/ncbi_entrez.html')
parser.add_argument('-e', action='store', dest='EmailAddress', required=True, help='Entrez requires your email address.')
parser.add_argument('-t', action='store', dest='SearchTerm', required=True, help='Requires a search term (wrap in double quotes).')

arguments = parser.parse_args()

Entrez.email = arguments.EmailAddress

SearchTerm = arguments.SearchTerm

#LookupCommand = "refseq[FILTER] AND txid9606[Organism] AND " + SearchTerm + "[TITL]"
LookupCommand = "refseq[FILTER] AND " + SearchTerm + "[TITL]"

handle = Entrez.esearch(db='protein', term=LookupCommand)

results = Entrez.read(handle)

handle.close()

#Lookup the FASTA sequence for each protein by its GeneInfo Identifier (GI) number
for gi in results['IdList']:
    handle = Entrez.efetch(db='protein', id=gi, rettype='fasta')

    print handle.read()

    handle.close()
