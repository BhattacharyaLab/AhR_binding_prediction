import os
from pybedtools import BedTool


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE_GENOME = 'hg19' #change to the desired reference genome (e.g., hg38)
REFERENCE_GENOME_URL = f'https://hgdownload.soe.ucsc.edu/goldenPath/{REFERENCE_GENOME}/bigZips/{REFERENCE_GENOME}.fa.gz'

DATA_DIR = os.path.join(BASE_DIR, '../data')
MOTIF_BED = 'motifs.bed.gz'
SEQ_FILE = 'sequences.csv.gz'


def set_name(peaks, name):
    peaks = peaks.cat(peaks).to_dataframe()
    
    peaks['name'] = [f'{name}_{x}' for x in + peaks.index]
    peaks = BedTool.from_dataframe(peaks)
    return peaks
    

def load_mcf7_peaks(name='MCF-7'):
    mcf7_1 = BedTool(os.path.join(DATA_DIR, 'chipseq/ahr_ahrr_mcf7/GSE90550_AHR-only_bound_peaks.bed'))
    mcf7_2 = BedTool(os.path.join(DATA_DIR, 'chipseq/ahr_ahrr_mcf7/GSE90550_AHR_AHRR_cobound_peaks.bed'))
    mcf7 = mcf7_1.cat(mcf7_2)
    
    mcf7 = set_name(mcf7, name)
    return mcf7