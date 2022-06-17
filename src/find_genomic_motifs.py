import gzip
import os
import re
import numpy as np
import pandas as pd

from pybedtools import BedTool

from Bio import SeqIO, Seq

from settings import DATA_DIR, MOTIF_BED, SEQ_FILE, REFERENCE_GENOME
import argparse

parser = argparse.ArgumentParser(description="Extract motif locations and motif flanking sequences from the reference genome for a transcription factor of interest")
parser.add_argument('--factor', type=str, default="AhR", help="Name of the factor (e.g. AhR, ER, GR)")
parser.add_argument('--forward-motif', type=str, default='GCGTG', help="Forward motif sequence (as python regex)")
parser.add_argument('--reverse-complement-motif', type=str, default='CACGT', help="Reverse complement motif sequence (as python regex)")
parser.add_argument('--flank-width', type=int, default=7, help='Number of flanking nucleotides to extract')

script_args = parser.parse_args()

factor = script_args.factor
forward_motif = script_args.forward_motif
reverse_complement_motif = script_args.reverse_complement_motif
flank_width = script_args.flank_width

factor = 'AhR'
forward_motif = 'GCGTG'
reverse_complement_motif = 'CACGC'

genome_path = os.path.join(
    DATA_DIR,
    f'genome_reference/{REFERENCE_GENOME}.fa.gz'
)
genome_readme = os.path.join(
    DATA_DIR,
    'genome_reference/README.txt'
)

motif_bed = os.path.join(
    DATA_DIR,
    f'motifs/{factor}-{MOTIF_BED}',
)

blacklisted_region_files = [
    os.path.join(DATA_DIR, 'blacklist/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz'),
    os.path.join(DATA_DIR, 'blacklist/ENCFF001TDO.bed.gz'),
]

if not(os.path.exists(genome_path)):
    raise ValueError(
        f'Reference genome sequence file "{genome_path}" does not exist.\n' +
        f'Please refer to "{genome_readme}" to download the appropriate reference genome file first.'
    )

def check_chromosome_placement(chrom_len, motif_pos, min_distance_from_chromosome_end):
    # ignore motifs too close to the chromosome boundaries
    return motif_pos[(motif_pos>min_distance_from_chromosome_end) & (motif_pos<chrom_len-min_distance_from_chromosome_end)]

def all_genomic_motifs(min_distance_from_chromosome_end=1500, flank_width=7, forward_motif=forward_motif, reverse_complement_motif=reverse_complement_motif):
    all_motifs = []

    motif_len = len(forward_motif)
    # ignore unnecessary 
    chroms = [x for x in genome_dict.genome_dict.keys() if '_' not in x and x != 'chrM']
    
    for chrom in chroms:
        print(f'Processing chromosome: {chrom}')

        # sequence for the entire chromosome
        seq = genome_dict.get_sequence(chrom, 0, None)
        chrom_len = len(seq)

        # finds all instances of forward and reverse complement motifs on the chromosome
        forward_motif_pos = np.array([m.start() + 1 for m in re.finditer(f'(?={forward_motif})', seq)])
        reverse_complement_motif_pos = np.array([m.start() + 1 for m in re.finditer(f'(?={reverse_complement_motif})', seq)])

        forward_motif_pos = check_chromosome_placement(chrom_len, forward_motif_pos, min_distance_from_chromosome_end)
        reverse_complement_motif_pos = check_chromosome_placement(chrom_len, reverse_complement_motif_pos, min_distance_from_chromosome_end)
        
        all_motifs += [(
                chrom, motif_pos, motif_pos + motif_len, 
                0, 1000, '+', 
                seq[motif_pos-flank_width:motif_pos+5+flank_width]) 
            for motif_pos in forward_motif_pos
        ]
        all_motifs += [(
                chrom, motif_pos, motif_pos + motif_len, 
                0, 1000, '-', 
                str(Seq.Seq(seq[motif_pos-flank_width:motif_pos+5+flank_width]).reverse_complement())) 
            for motif_pos in reverse_complement_motif_pos
        ]
    
    
    all_motifs = pd.DataFrame(
        data=all_motifs, 
        columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'seq']
    ).sort_values(by=['chrom', 'start']).reset_index(drop=True)   
    all_motifs['name'] = all_motifs.index
    seqs = all_motifs[['seq']]
    all_motifs.pop('seq')
    
    all_motifs = BedTool.from_dataframe(all_motifs)
    return all_motifs, seqs

class GenomeDict(object):
    def __init__(self, path_to_genome_fasta):
        # input: path_to_genome_fasta (path the the genome fasta file (either .fa or .fa.gz))
        #
        # creates a dictionary of Bio.SeqRecord objects for each chromosome
        #
        # Warning: loads the entire genome sequence into memory at once
        #          takes a bit of time initially
        #          loads individual sequences quickly
        fasta_open = None

        if path_to_genome_fasta.endswith('.gz'):
            fasta_open = gzip.open
        else:
            fasta_open = open

        with fasta_open(path_to_genome_fasta, "rt") as handle:
            genome_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))

        self.genome_dict = genome_dict
    
    def get_sequence(self, chrom, start, end):
        return str(self.genome_dict[chrom].seq[start:end].upper())

genome_dict = GenomeDict(genome_path)

# load blacklisted regions
blacklisted_regions = [BedTool(blacklisted_file) for blacklisted_file in blacklisted_region_files] 
blacklisted_regions = BedTool.cat(*blacklisted_regions)

all_motifs, seqs = all_genomic_motifs(
    forward_motif=forward_motif,
    reverse_complement_motif=reverse_complement_motif,
    flank_width=flank_width,
)

genome_dict = None

all_motifs = all_motifs.subtract(blacklisted_regions, A=True)
seqs = seqs.loc[all_motifs.to_dataframe()['name']]

all_motifs = all_motifs.to_dataframe()

all_motifs.to_csv(
    os.path.join(
        DATA_DIR, 'motifs', f'{factor}-{MOTIF_BED}'
    ), compression='gzip', index=None, sep='\t', header=None
)

seqs.to_csv(
    os.path.join(
        DATA_DIR, 'motifs', f'{factor}-{SEQ_FILE}'
    ), index=None
)
