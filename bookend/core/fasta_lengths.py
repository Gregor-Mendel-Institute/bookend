import sys
from cython_utils._fasta_utils import import_genome
'''
Outputs a two-column tab-delimited file from a FASTA file:
    Feature name    Length of the feature sequence

Commandline arguments: [1] genome_FASTA
genome_FASTA - Full path of the FASTA file to search
'''

genome_FASTA = sys.argv[1]
split_on=' '

chrom='none'
genome=import_genome(genome_FASTA)
chromosomes = {}
for k,v in genome.items():
    chromosomes[k]=len(v)

for k,v in sorted(list(chromosomes.items())):
    print('\t'.join([k,str(v)]))

