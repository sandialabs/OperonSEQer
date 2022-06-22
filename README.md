(c) 2021: National Technology & Engineering Solutions of Sandia, LLC (NTESS)
# OperonSEQer

OperonSEQer takes RNA-sequencing data and transforms to output a Kruskal-Wallis statistic and p-value for similarity in expression of contiguous gene pairs and the intergentic region. Following this, the data is run through a suite of 6 machine learning models which determine whether each gene pair is an operon pair or not. Finally, a voting system is implemented to establish confidence in the operon pair calls. 


## Requirements

OperonSEQer is designed to be run on Linux or Mac platforms and requires Python 3.7 or greater. 
See setup.py file for specific requirements (they will be automatically installed if the protocol below is followed, but can also be installed manually)

## Install

It is highly recommended to install OperonSEQer in a virtual environment such as conda:

**conda create -n OperonSEQer_env python=3.7**

**conda activate OperonSEQer_env**

**cd \path\to\cloned\repository**

**pip install -r requirements.txt**

**pip install .**

Note: If you are having trouble downloading dependencies through pip install (due to firewall or proxy), manually download the modules listed in the requirements file

## Quick Start

**python OperonSEQer -h**

~~~
usage: OperonSEQer [-h] -c COVFILE -g GFF -o OUT [-p PREDS] [-k] [-t PATH] [-d DELIMITER]

OperonSEQer

optional arguments:
  -h, --help                          show this help message and exit
  -c COVFILE, --coverage COVFILE      scaled coverage file from RNA-seq
  -g GFF, --gff GFF                   modified gff file for organism (chr archive type start stop . strand . geneName)
  -o OUT, --output_name OUT           output prefix
  -t THRESH, --threshold THRESH       threshold for number of calls to become an operon (default is 3)
  -p PREDS, --preds PREDS             prediction file (optional)
  -k, --pickleit                      if true, saves the processed data as a pickle file
  -t PATH, --path PATH                optional folder path where input files are and output files go (if not specified, current directory is used)
  -d DELIMITER, --delimiter DELIMITER optional delimiter for gff file (if not specified, tab is assumed)
~~~

## Example

In the folder ExampleFiles, the following input files are provided as an example:

**SRR10775302_sort.cov** - Coverage file generated from SRR10775302 (Culviner et al, 2020, DOI: https://doi.org/10.1128/mBio.00010-20)

**MicrobesOnline_Operon_Ecoli_MG1655.txt** - Operon prediciton downloaded from Microbes Online (http://www.microbesonline.org/, Dehal et al, 2009, DOI: 10.1093/nar/gkp919)

**Escherichia_coli_mg1655_lite.txt** - GFF file modified from original downloaded at https://bacteria.ensembl.org/index.html

Use the following to obtain output files from the example:

python OperonSEQer -c SRR10775302_sort.cov -g Escherichia_coli_mg1655_lite.txt -p MicrobesOnline_Operon_Ecoli_MG1655.txt -k -d space [-t if you want to specify an output path please do so]

## Expected output files

Data_CoverageExtracted.txt - intermediate file with extracted coverages

Data_predsMerged.txt - intermediate file with predictions merged if a prediction was provided

Data_Formatted.txt - intermediate file with data formatted (and a Data_Formatted.p pickle file if requested)

[Model]_result.csv - specs for each of the six models (XGB, MLP, SVM, RF, LR and GNB) if a prediction was provided

[Model]_predictions.csv - operon pair predictions for each of the six models

[Model] ROC curve - ROC curve for each of the six models

OperonSEQer_vote_result.csv - final vote tally of operon pair calls

## For advanced users:

A)
Individual scripts are provided for more advanced users that want to run individual parts of the pipeline. 

1) extractCoverageInts.py
2) mergepreds.py
3) formatData.py
4) OperonSEQer.py
5) voting.py

B)
A folder called TrainingModels contains the code to train all six ML models with new data
