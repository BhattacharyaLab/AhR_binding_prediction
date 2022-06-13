import os

from utils import download_file
from settings import REFERENCE_GENOME_URL

genome_name = os.path.basename(REFERENCE_GENOME_URL)
print(f'Downloading reference genome "{genome_name}"')
download_file(REFERENCE_GENOME_URL, dest=f'../data/genome_reference/{genome_name}')
