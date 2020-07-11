# ExPecto
This repository contains code for predicting expression effects of human genome variants ab initio from sequence with ExPecto models and training new sequence-based expression model with any expression profile.

The ExPecto framework is described in the following manuscript: Jian Zhou, Chandra L. Theesfeld, Kevin Yao, Kathleen M. Chen, Aaron K. Wong,  and Olga G. Troyanskaya, [Deep learning sequence-based ab initio prediction of variant effects on expression and disease risk](https://www.nature.com/articles/s41588-018-0160-6), Nature Genetics (2018).


## Install
#### Important Note: that for training new models with `train.py`, the default hyperparameters are only compatible with xgboost version 0.7.post4, because in newer xgboost versions the interpretation of eta parameter have been substantially different. Please make sure to install the version 0.7.post4 `pip install xgboost=0.7.post4`. The default hyperparameters are not compatible with newer xgboost versions. 

Clone the repository first:
```bash
git clone https://github.com/FunctionLab/ExPecto.git
cd ExPecto
sh download_resources.sh; tar xf resources_20190807.tar.gz
```
Install PyTorch following instructions from https://pytorch.org/.  Use `pip install -r requirements.txt` to install the other dependencies.


 
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
closest-features --delim '\t' --closest --dist <(awk '{printf $1"\t"$2-1"\t"$2"\n"}' ./example/example.vcf|sed s/chr//g|sed s/^/chr/g|sort-bed - ) ./resources/geneanno.pc.sorted.bed > ./example/example.vcf.bed.sorted.bed.closestgene
```

`--snpEffectFilePattern ./example/example.vcf_shiftSHIFT_outdir/infile.vcf.wt2100.fasta.ref.h5.diff.h5` specifies the name pattern of the input epigenomic effect prediction files. Note these files are the output from `chromatin.py`. `SHIFT` string is a placeholder that is substituted automatically to the shift positions (e.g. 0, -200, -400, ...). 


Optional:  For very large input files use the split functionality to distribute the prediction into multiple runs. For `predict.py` you can use for example `--splitFlag --splitIndex 0 --splitFold 10` to divide the input into 10 chunks and process only the first chunk.

### Training example:
```bash
python ./train.py --expFile ./resources/geneanno.exp.csv --targetIndex 1 --output model.adipose
```

This trains an ExPecto model using the Adipose gene expression profile in the first column of the `geneanno.exp.csv` file and the default precomputed epigenomic features. For training new ExPecto model for your custom (differential) expression profile, replace geneanno.exp.csv with your expression profile. The gene order has to be the same as the geneanno.csv. 

The new trained model(s) can be used by `predict.py` by adding the path of the xgboost model file to a new `modelList` file. The new models should be put in a separate modellist file not mixed with provided models, because the provided models were in an old legacy format incompatible with new trained models.

*** Note that for training new models with train.py, the default hyperparameters are only compatible with xgboost version 0.7.post4. Make sure to install the correct version e.g. pip install xgboost=0.7.post4. The default hyperparameters are not compatible with newer xgboost versions. ***


### Contact me
Jian Zhou [jzhoup@gmail.com](mailto:jzhoup@gmail.com)



### Use agreement

If you are interested in obtaining the software for commercial use, please contact Office of Technology Licensing, Princeton University (Laurie Tzodikov 609-258-7256, tzodikov@princeton.edu, or Linda Jan, 609-258-3653,  ljan@princeton.edu). For academic use, downloading or using the software means you agree with the following Academic Use SOFTWARE Agreement.

```
PRINCETON Academic Use SOFTWARE Agreement

The Trustees of Princeton University, a non-profit educational corporation organized and existing under the
laws of the State of New Jersey with its Office of Technology Licensing at 87 Prospect Avenue, Princeton NJ
08544 (“PRINCETON”) is willing to make the ExPecto Software (software for predicting expression effects of 
genome variants ab initio from sequence ) (“SOFTWARE”) available to you
(“RECIPIENT”) under the following terms and conditions: 

1. The above SOFTWARE is the property of PRINCETON and is made available as a service to the research
community. The SOFTWARE is a research tool still in the development stage and is being provided “as is”
without any support, services or improvements. PRINCETON makes no representations and extends no
warranties of any kind, and expressly disclaims any representations or warranties of any kind
(including, but not limited to, any warranties of merchantability, fitness for a particular purpose, or
non-infringement).
2. The SOFTWARE will be used for teaching or not-for-profit research purposes only. The RECIPIENT agrees
not to use the Software for commercial purposes, or for diagnosis, treatment, cure, prevention or mitigation of
disease or any other conditions in man. The RECIPIENT agrees that the Software is not intended to substitute
for care by a licensed healthcare professional.
3. The SOFTWARE will not be further distributed to others without PRINCETON’S prior written consent.
The RECIPIENT shall promptly refer any request for the SOFTWARE to PRINCETON.
4. RECIPIENT
acknowledges that any programs or software created based on the SOFTWARE will be considered a
derivative of SOFTWARE and owned by PRINCETON.
5. The RECIPIENT agrees to acknowledge the source of the SOFTWARE in any publications.
6. RECIPIENT will use the SOFTWARE in compliance with all applicable laws, policies and regulations
including, but not limited to, any approvals, informed consent, patient confidentiality principles and US
government and local export control regulations. RECIPIENT acknowledges that the SOFTWARE may
not be exported to Cuba, Iran, North Korea, or Syria.
7. In no event shall PRINCETON be responsible or liable to RECIPIENT or any third party for
RECIPIENT’s activities under or related to this Agreement, including RECIPIENT’S use of the
SOFTWARE. RECIPIENT agrees to indemnify, defend, and hold harmless PRINCETON (including its
trustees, officers, faculty, employees, students and agents), from and against any and all claims,
liabilities, or damages based on, arising out of, or relating to this Agreement, including RECIPIENT’s
use of the SOFTWARE or RECIPIENT’S breach of this Agreement.
8. All disputes regarding the construction, interpretation and the parties’ obligations under this Agreement
shall be governed by the laws of the State of New Jersey, notwithstanding any of that state’s laws to the
contrary. The venue and jurisdiction for the resolution of any such disputes shall be Mercer County,
New Jersey
```
