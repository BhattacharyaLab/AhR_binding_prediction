import requests
import os

import pandas as pd
import numpy as np

from utils import download_file

import argparse

download_with_python = False
download_with_wget = True

parser = argparse.ArgumentParser(description="download input files")
parser.add_argument('--download-with', default='python')

script_args = parser.parse_args()
download_with_python = script_args.download_with == 'python'
download_with_wget = script_args.download_with == 'wget'

def get_encode_filesize(encode_file):
    url = f'https://www.encodeproject.org/files/{encode_file}/'
    r = requests.get(url, headers={'accept': 'application/json'})
    response = r.json()
    return response['file_size']

if download_with_python:
    print('Downloading all input files through python.\n\n')
else:
    #print('Downloading all input files through wget. Will only print wget commands.\n\n')
    pass

filetypes = ['DNase', 'HMs', 'TFs']
for filetype in filetypes:
    filelist = pd.read_csv(f'../data/encode_filelists/{filetype}-ENCODE_file-accessions.csv', index_col=0).replace(float('nan'), '')
    cells = [x.split('__')[0] for x in filelist.columns]
    
    for cell in cells:
        all_encode_files = []
        path = f'../data/download/{cell}/{filetype}_files'
        try:
            os.makedirs(path)
        except:
            pass

        for encode_file_desc in filelist[f'{cell}__files'].items():
            factor, encode_files  = encode_file_desc
            if not(encode_files):
                continue
            all_encode_files.extend([(
                    encode_file,
                    f'https://www.encodeproject.org/files/{encode_file}/@@download/{encode_file}.bigWig', 
                    f'{path}/{factor}_{encode_file}.bigWig',
                ) for encode_file in encode_files.split(',')
            ])
        if download_with_python:
            for encode_filename, url, dest in all_encode_files:
                print(f'Downloading {dest} input file.')
                download_file(url, dest, file_size=get_encode_filesize(encode_filename))
        elif download_with_wget:
            for encode_filename, url, dest in all_encode_files:
                print(f'wget {url} -O {dest}')
