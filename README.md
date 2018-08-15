# ExPecto
This repository contains code for predicting expression effects of human genome variants ab initio from sequence with ExPecto models and training new sequence-based expression model with any expression profile.


## Install
Clone the repository first:
```bash
git clone https://github.com/FunctionLab/ExPecto.git;
cd ExPecto;
```
Use `pip install -r requirements.txt` to install the dependencies. If you have a CUDA-enabled GPU, we highly recommend installing pytorch with GPU support following instructions from http://pytorch.org/.

Run `sh download_resources.sh; tar xf resources.tar.gz` to download and extract necessary model files and chromatin representations for training new ExPecto models. 
## Usage

### Example:
```bash
python chromatin.py ./example/example.vcf
python predict.py --coorFile ./example/example.vcf --geneFile ./example/example.vcf.bed.sorted.bed.closestgene --snpEffectFilePattern ./example/example.vcf.shift_SHIFT.diff.h5 --modelList ./resources/modellist --output output.csv
```

The output will be saved to output.csv. The first few columns of the csv file will be the same as the vcf files. The additional columns include predicted expression effect (log fold change) for each of the input models in the order given by the modelList file. 

If `./resources/hg19.fa.flat` cannot be loaded, try removing it and it will be generated next time you run the code. 

#### Details:


`chromatin.py` computes the chromatin effects of the variants, using trained convolutional neural network model. The input vcf files need to contain only variant with a single alternative allele. Variants such as `A T,AT` is not recognized and you can split it into biallelic variants to run it.


`predict.py` computes predicted tissue-specific expression effects which takes predicted chromatin effects as input.

`--coorFile ./example/example.vcf` specifies the variants of interest in vcf format (coordinates should be in hg19, support for other genomes can be obtained via switching genome fasta file, but note that the current models are trained on hg19 genome).

`--closestGeneFile ./example/example.vcf.bed.sorted.bed.closestgene` specifies the gene association file which decides for each variant the associated gene for which the expression effect is predicted. The content of the gene association file has to include the following information: the first column and third columns are chromosome names and positions, and the last three columns are the strand of the associated gene, the ENSEMBL gene id (matched with the gene annotation file `./resources/geneanno.csv`) and distance to the representative TSS of that gene. The row order gene association file does not need to be in the same as the vcf file. The distance should be signed and calculated as '' TSS position - variant position" regardless of on which strand the gene is transcribed. The representive TSSes can be found in the provided geneanno.csv file. The associated gene can be specified by finding the closest representative TSS. When is known of gene is of interest, such as for eQTL predictions, the know gene association can be used. This can be done for example using closest-features from [BEDOPS](https://bedops.readthedocs.io/en/latest/) and the representation TSS of protein coding genes that we included, for example:
```bash
closest-features --delim '\t' --closest --dist <(awk '{printf $1"\t"$2-1"\t"$2"\n"}' ./example/example.vcf|sort-bed - ) ./resources/geneanno.pc.sorted.bed > ./example/example.vcf.bed.sorted.bed.closestgene
```

`--snpEffectFilePattern ./example/example.vcf_shiftSHIFT_outdir/infile.vcf.wt2100.fasta.ref.h5.diff.h5` specifies the name pattern of the input epigenomic effect prediction files. Note these files are the output from `chromatin.py`. `SHIFT` string is a placeholder that is substituted automatically to the shift positions (e.g. 0, -200, -400, ...). 


Optional:  For very large input files use the split functionality to distribute the prediction into multiple runs. For `predict.py` you can use for example `--splitFlag --splitIndex 0 --splitFold 10` to divide the input into 10 chunks and process only the first chunk.

### Training example:
```bash
python ./train.py --expFile ./resources/geneanno.exp.csv --targetIndex 1 --output model.adipose
```

This trains an ExPecto model using the Adipose gene expression profile in the first column of the `geneanno.exp.csv` file and the default precomputed epigenomic features. For training new ExPecto model for your custom (differential) expression profile, replace geneanno.exp.csv with your expression profile. The gene order has to be the same as the geneanno.csv. The generated model can be used by `predict.py` by adding the path of the xgboost model file to the `modelList` file.


### Citation

The ExPecto framework is described in the following manuscript: Jian Zhou, Chandra L. Theesfeld, Kevin Yao, Kathleen M. Chen, Aaron K. Wong,  and Olga G. Troyanskaya, Deep learning sequence-based ab initio prediction of variant effects on expression and disease risk

### Contact me
Jian Zhou [jzhoup@gmail.com](mailto:jzhoup@gmail.com)
