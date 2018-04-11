This repository contains code for predicting expression effects for human genome variants with ExPecto models and training new sequence-based expression model with any expression profile. 

The ExPecto framework is described in the following manuscript: Jian Zhou, Chandra L. Theesfeld, Kevin Yao, Kathleen M. Chen, Aaron K. Wong,  and Olga G. Troyanskaya, Deep learning sequence-based ab initio expression prediction and disease-risk identification

## Install
Use `pip install -r requirements.txt` to install the dependencies.
Run `sh download_resources.sh; tar xf resources.tar.gz` to download and extract necessary model files and chromatin representations for training new ExPecto models.

## Usage

##### Example run :
```bash
python chromatin.py ./example/example.vcf
python predict.py --coorFile ./example/example.vcf --geneFile ./example/example.vcf.bed.sorted.bed.closestgene --snpEffectFilePattern ./example/example.vcf.shift_SHIFT.diff.h5 --modelList ./resources/modellist --output output.csv
```

The output will be saved to output.csv. The first few columns of the csv file will be the same as the vcf files. The additional columns include predicted expression effect (log fold change) for each of the input models in the order given by the modelList file. `chromatin.py` computes the chromatin effects of the variant both on-site and to nearby regions, using trained convolutional neural network model. `predict.py` computes predicted tissue-specific expression effects from the predicted chromatin effects.

#####  Explanations for the arguments:

`--coorFile ./example/example.vcf` includes the variants of interest in vcf format.
`--closestGeneFile ./example/example.vcf.bed.sorted.bed.closestgene` is the gene association file which specify the associated gene for each variant and the expression effect is predicted wrt that gene. The content of the gene association file has to have to following information: the first column and third columns are chromosome names and positions, and the last three columns are the strand of the associated gene, the ENSEMBL gene id (matched with the provided geneanno.csv file) and distance to the representative TSS of that gene. The gene association file does not need to be in the same order as the vcf file. The distance is signed and calculated as '' TSS position - variant position" regardless of on which strand the gene is transcribed. The representive TSSes can be found in the provided geneanno.csv file. The associated gene can be specified by finding the closest representative TSS. When is known of gene is of interest, such as for eQTL predictions, the know gene association can be used. This can be done for example using closest-features from [BEDOPS](https://bedops.readthedocs.io/en/latest/) and the representation TSS of protein coding genes that we included, for example:
```
closest-features --delim '\t' --closest --dist <(awk '{printf $1"\t"$2-1"\t"$2"\n"}' ./example/example.vcf|sort-bed - ) ./resources/geneanno.pc.sorted.bed > ./example/example.vcf.bed.sorted.bed.closestgene
```

`--snpEffectFilePattern ./example/example.vcf_shiftSHIFT_outdir/infile.vcf.wt2100.fasta.ref.h5.diff.h5` specify the name pattern of the input epigenomic effect prediction files. `SHIFT` string is used as a placeholder that is substituted automatically to the shift positions (e.g. 0, -200, -400, ...). For generating the epigenomic effect predictions, use the scripts under ./convnet directory and see instructions. Note that we will soon release a pure python version of convnet.

Optional:  For very large input files use the split functionality to distribute the prediction into multiple runs. For example, `--splitFlag --splitIndex 0 --splitFold 10` will divide the input into 10 chunks and process only the first chunk.

##### For training new models:

Example run:
```bash
python ./train.py --expFile ./resources/geneanno.exp.csv --targetIndex 1 --output model.adipose
```

This trains an ExPecto model using the Adipose gene expression profile in the first column of the `geneanno.exp.csv` file. `./Xreducedall.2002.npy` is the default precomputed epigenomic features. For training new ExPecto model(s) for your custom (differential) expression profile, replace geneanno.exp.csv with the expression profile. Each row in the file has to be the expression value of gene. The gene order has to be the same as the geneanno.csv. The generated model can be used by `predict.py` by adding the path of the model file ('.save') to the `modelList` file.


##### Contact me: 
Jian Zhou [jzhoup@gmail.com](mailto:jzhoup@gmail.com)
