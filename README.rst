Dependencies 

---------------------------------------------------------------- 

The code for this project has been written in Python 3 and requires the Python dependencies listed below. The dependencies are also available in the *requirements.txt* file and can be installed with the *pip install -r requirements.txt* 

 

* biopython 

* pandas 

* numpy

* pybedtools

* numpy 

 

 

Download the reference genome 

---------------------------------------------------------------- 

 

In order to create a bed file list of viable motifs, we first must download the hg19 reference genome. 

 

In case you would like to perform the pipeline with a different reference genome, change the following settings in src/settings.py 

.. code-block:: python 

    REFERENCE_GENOME = 'hg19'.  # change to the desired reference genome 

    REFERENCE_GENOME_URL = ' https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz' # URL to the appropriate "*.fa.gz" file 

 

To download the reference genome, issue the following command 

.. code-block:: 

    python src/download_reference.py 

 

 

Download the input data 

---------------------------------------------------------------- 

 

The repository comes preloaded with the AhR ChIP-seq bed files, as well as DNase-seq bed files for each investigated cell line/type. These files are contained in the *data/chipseq/* directory. 

All other input data including DNase-seq, transcription factor (TF) and histone modification (HM) *bigWig* signal tracks need to be downloaded.  The lists of all tracks used are located in *data/encode_filelists/* **Note:** these files will require several hundred GB of space. 

 

To download the input data, issue the following command 

 

.. code-block:: 

    python src/download_input_data.py 

 


Create the AHR motif (DRE) bedfile 

---------------------------------------------------------------- 

Creates a bed file containing all AhR motifs – Dioxin Response Elements (DREs) in the reference genome. 

 

To create a motif file, issue the following command 

 

.. code-block:: 

        python src/find_genomic_motifs.py --factor AhR --forward-motif GCGTG --reverse-complement-motif CACGC --flank_width 7 

 

 

Create a list of training/validation AHR motifs  

---------------------------------------------------------------- 

Creates a list of singleton bound (only one DRE under the AhR peak) and isolated unbound DREs (at least 100bps away from any other DRE) in open chromatin, as well as DREs under multi-DRE bound AhR peaks, and five dummy DREs for each 0-DRE AhR peak. 

 

 

Create the files containing training/validation AhR motifs for each cell line/type 

------------------------------------------------------------------------------- 

To create the list of training/validation AhR motifs to be used in the construction of input features used in the subsequent machine learning steps, issue the following commands 

 

.. code-block:: 

    python src/prepare_ML_motifs.py --primary_factor MCF-7 --primary_file data/chipseq/MCF-7_ahr_ahrr/{name_of_file} 

    python src/prepare_ML_motifs.py --primary_factor primary_hepatocytes --primary_file data/chipseq/hepatocytes_ahr/{name_of_file} 

    python src/prepare_ML_motifs.py --primary_factor HepG2 --primary_file data/chipseq/HepG2_ahr/{name_of_file} 

    python src/prepare_ML_motifs.py --primary_factor GM17212 --primary_file data/chipseq/GM17212_ahr/{name_of_file} 

 

 

Prepare input features 

------------------------------ 

To extract DNase-seq, HM, and TF features for each cell line/type, issue the following commands 

 

. code-block:: 

    python src/prepare_ML_features.py --cells MCF-7 

    python src/prepare_ML_motifs.py --cells primary_hepatocytes 

    python src/prepare_ML_motifs.py --cells HepG2 

    python src/prepare_ML_motifs.py --cells GM17212 

 

 

 

 

 
