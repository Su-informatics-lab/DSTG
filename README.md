# Deconvoluting Spatial Transcriptomics data through Graph-based convolutional networks (DSTG)

This is a TensorFlow implementation of DSTG for decomposing spatial transcriptomics data, which is described in our paper: 

## Installation

```bash
python setup.py install
```

## Requirements
* tensorflow (>0.12)
* networkx

## Run the demo

load the example data using the convert_data.R script
In the example data, we provide two synthetic spatial transcriptomics data generated from scRNA-seq data (GSE72056). Each synthetic data consists of 1,000 spots, which can be found in folder synthetic_data.
```bash
cd DSTG
Rscript convert_data.R # load example data 
python train.py # run DSTG
```
Predicted compositions within each spot are saved in will be shown in the DSTG_Result folder.

Performance of JSD score will be shown if you run
```
Rscript evaluation.R
```
If you want to use your own scRNA-seq data to deconvolute your spatail transcriptomcis data, provide you data to script below:

## Run your own data
When using your own scRNA-seq data to deconvolute your spatail transcriptomcis data, you have to provide 
* the raw scRNA-seq data matrix and label, which are saved as .RDS format (e.g. 'scRNAseq_data.RDS' & 'scRNAseq_label.RDS')
* the raw spatial transcriptomics data matrix saved as .RDS format (e.g. 'spatial_data.RDS')

```
cd DSTG
Rscript  convert_data.R  'scRNAseq_data.RDS'  'spatial_data.RDS'  'scRNAseq_label.RDS'
python train.py # run DSTG
```
Then you will get your results in the DSTG_Result folder.


## Cite

Please cite our paper if you use this code in your own work:

```
J. Su, Q. Song, DSTG: Deconvoluting Spatial Transcriptomics Data through Graph-based Artificial Intelligence, bioRxiv, (2020) 
```
