#!/usr/bin/env python

import argparse, sys, collections
import numpy as np

###############################################################################
def parse_args(args):
    '''
    parses command line arguments
    '''
    parser = argparse.ArgumentParser(description ='Filters MAZTER-seq output suppl table')

    parser.add_argument('--bed', '-b', type = str, help = 'input sk1 bed12 file')

    parser.add_argument('--annotation', '-a', type = str, help = 'Schwartz transcriptome annotation table file (Supp Table 5)')

    parser.add_argument('--m6a', '-m', type = str, help = 'Schwartz m6a sites table file (Supp Table 1)')

    parser.add_argument('--mazter', '-i', type = str, help = 'Schwartz mazter m6a sites table file (Supp Table 4)')

    parser.add_argument('--outfile', '-o', type = str, help = 'output fasta file')

    return parser.parse_args()
###############################################################################

###############################################################################
class Annotation:
    def __init__(self, name, chrom, start, stop, strand):
        self.name = name
        self.start = start
        self.stop = stop
        self.strand = strand
        self.chrom = chrom
###############################################################################

###############################################################################
def parse_schwartz_annotation(infile):
    '''
    '''
    ref_coords = collections.defaultdict()
    for gene_symbol, chrom, start, stop, strand in read_schwartz_annotation_tsv(infile):
        ref_coords[gene_symbol] = Annotation(gene_symbol, chrom, start, stop, strand)

    return ref_coords
###############################################################################

###############################################################################
def read_schwartz_annotation_tsv(infile):
    '''
    Generator script used to parse the tsv version of Supp Table 5 from Schwartz et al Cell 2013
    input : tsv file
    yeilds the gene_symbol, chr name, CDS start coordinate, CDS end coordinate, and the strand
    for each annotated gene in the table
    '''
    with open(infile, 'r') as tsv:
        #Reads the header line
        line = tsv.readline()
        #Begins reading the data lines in the tsv
        for line in tsv:
            line = line.strip().split('\t')
            gene_symbol = line[0].strip()
            chrom = line[1].strip()
            strand = line[6].strip()
            try:
                cds_start = int(line[4].strip())
            except:
                cds_start = float('NaN')

            try:
                cds_end = int(line[5].strip())
            except:
                cds_end = float('NaN')

            yield gene_symbol, chrom, cds_start, cds_end, strand
###############################################################################

##############################################################################
class Sites:
    def __init__(self, name, chrom, site):
        self.name = name
        self.chrom = chrom
        self.site = site
###############################################################################

###############################################################################
def parse_m6a_tsv(infile, parse_type=''):
    '''
    A switch function that determines which of the two m6A site files to read (m6A-seq or MAZTER-seq)
    and uses the correct reading function. Then converts the site information from either file type into
    a single unified site object so that the data can be used in the same way for both file types
    '''
    m6a_coords = collections.defaultdict(list)

    if parse_type == 'sites':
        for gene_symbol, chrom, peak in read_m6a_tsv(infile):
            m6a_coords[gene_symbol].append(Sites(gene_symbol, chrom, peak))

    elif parse_type == 'mazter':
        for gene_symbol, chrom, group, site in read_mazter_tsv(infile):
            #Selects MAZTER-seq sites that have a confidence group greater than 1
            #This was done as in Leger et al Nature Com 2021
            if group > 1:
                m6a_coords[gene_symbol].append(Sites(gene_symbol, chrom, site))
    else:
        sys.stderr.write("No parse type was specified")

    return m6a_coords
###############################################################################

###############################################################################
def read_m6a_tsv(infile):
    '''
    Generator script used to parse the tsv version of Supp Table 1 from Schwartz et al Cell 2013
    input : tsv file
    yeilds the gene_symbol, chr name, and the m6A peak site for each site in the table
    '''
    with open(infile, 'r') as tsv:
        line = tsv.readline()
        for line in tsv:
            line = line.strip().split('\t')
            gene_symbol = line[0].strip()
            chrom = line[2].strip()

            try:
                peak = int(line[3].strip())
            except:
                peak = float('NaN')

            yield gene_symbol, chrom, peak
###############################################################################

###############################################################################
def read_mazter_tsv(infile):
    '''
    Generator script used to parse the tsv version of Supp Table 4 from Garcia-Campos et al Cell 2019
    input : tsv file
    yeilds the gene_symbol, chr name, confidence group, and the m6A peak site for each site in the table
    '''
    with open(infile, 'r') as tsv:
        line = tsv.readline()
        for line in tsv:
            line = line.strip().split('\t')
            chrom = line[1].strip()
            gene_symbol = line[23].strip()
            try:
                confidence_group = int(line[22].strip())
            except:
                confidence_group = float('NaN')

            try:
                peak = int(line[2].strip())
            except:
                peak = float('NaN')

            yield gene_symbol, chrom, confidence_group, peak
###############################################################################

###############################################################################
def calc_indexes(site_coords, annotations):
    '''
    Converts the m6A site from the reference coordinates used by Schwartz et al Cell 2013 
    to the sk1 MV0 reference genome coordinates based on the distance between the start codons
    '''
    indexes = collections.defaultdict(list)

    site_genes = set(site_coords.keys())
    annotations_genes = set(annotations.keys())
    genes = site_genes.intersection(annotations_genes)

    for gene in genes:
        site_data = site_coords[gene]
        ref_data = annotations[gene]

        for site in site_data:
            index = site.site - ref_data.start
            indexes[gene].append(index)

    return indexes
###############################################################################

###############################################################################
def convert_m6a_coords(m6a_indexes, mazter_indexes, inbed, outfile):
    '''
    '''

    with open(outfile, 'w') as out:
        with open(inbed) as bed:
            for line in bed:
                line = line.strip().split('\t')

                chrom = line[0].strip()
                start = int(line[1])
                gene = line[3].strip().rsplit('_', 1)[0]
                score = line[4].strip()
                strand = line[5].strip()

                if gene in m6a_indexes:
                    for site in m6a_indexes[gene]:
                        peak = start + site
                        m6a_start = str(peak-1)
                        m6a_stop = str(peak)
                        outstring = f'{chrom}\t{m6a_start}\t{m6a_stop}\t{gene}_m6a-seq\t{score}\t{strand}'
                        out.write(f'{outstring}\n')

                if gene in mazter_indexes:
                    for site in mazter_indexes[gene]:
                        #Remove the comment and tab the remaining code if you don't want to keep MAZTER-seq dublicates
                        #if site not in m6a_indexes[gene]:
                        peak = start + site
                        m6a_start = str(peak-1)
                        m6a_stop = str(peak)
                        outstring = f'{chrom}\t{m6a_start}\t{m6a_stop}\t{gene}_mazter\t{score}\t{strand}'
                        out.write(f'{outstring}\n')
###############################################################################

###############################################################################
def main(args):
    #Parse the inputs args/options
    options = parse_args(args)

    annotations = parse_schwartz_annotation(options.annotation)
    m6a_sites = parse_m6a_tsv(options.m6a, "sites")
    mazter_sites = parse_m6a_tsv(options.mazter, 'mazter')

    m6a_indexes = calc_indexes(m6a_sites, annotations)
    mazter_indexes = calc_indexes(mazter_sites, annotations)

    convert_m6a_coords(m6a_indexes, mazter_indexes, options.bed, options.outfile)
###############################################################################

if (__name__ == "__main__"):
    main(sys.argv)
    raise SystemExit

