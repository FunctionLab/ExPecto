# -*- coding: utf-8 -*-
"""Predict variant expression effects

This script takes the predicted chromatin effects computed by chromatin.py and
expression model file list, and predicts expression effects in all models provided
in the model list.

Example:
        $ python predict.py --coorFile ./example/example.vcf --geneFile ./example/example.vcf.bed.sorted.bed.closestgene --snpEffectFilePattern ./example/example.vcf.shift_SHIFT.diff.h5 --modelList ./resources/modellist --output output.csv
For very large input files use the split functionality to distribute the
prediction into multiple runs. For example, `--splitFlag --splitIndex 0 --splitFold 10`
will divide the input into 10 chunks and process only the first chunk.

"""
import argparse
import xgboost as xgb
import pandas as pd
import numpy as np
import h5py
from six.moves import reduce

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--coorFile', action="store", dest="coorFile")
parser.add_argument('--geneFile', action="store",
                    dest="geneFile")
parser.add_argument('--snpEffectFilePattern', action="store", dest="snpEffectFilePattern",
                    help="SNP effect hdf5 filename pattern. Use SHIFT as placeholder for shifts.")
parser.add_argument('--modelList', action="store", dest="modelList",
                    help="A list of paths of binary xgboost model files (if end with .list) or a combined model file (if ends with .csv).")
parser.add_argument('--nfeatures', action="store",
                    dest="nfeatures", type=int, default=2002)
parser.add_argument('--output', action="store", dest="output")
parser.add_argument('--fixeddist', action="store",
                    dest="fixeddist", default=0, type=int)
parser.add_argument('--maxshift', action="store",
                    dest="maxshift", type=int, default=800)
parser.add_argument('--batchSize', action="store",
                    dest="batchSize", type=int, default=500)
parser.add_argument('--splitFlag', action="store_true", default=False)
parser.add_argument('--splitIndex', action="store",
                    dest="splitIndex", type=int, default=0)
parser.add_argument('--splitFold', action="store",
                    dest="splitFold", type=int, default=10)
parser.add_argument('--threads', action="store", dest="threads",
                    type=int, default=16, help="Number of threads.")
args = parser.parse_args()


def compute_effects(snpeffects, snpdists, snpstrands, all_models, maxshift=800, nfeatures=2002, batchSize=500,old_format=False):
    """Compute expression effects (log fold-change).

    Args:
        snpeffects: list of chromatin effect numpy arrays
        snpdists:  integer array or pandas Series representing distances to TSS
        snpstrands: string array or pandas Series containing only '+' and '-'s
                     representing the strand of the TSS for each variant
        all_models: list of ExPecto model files.
        maxshift:  maximum shift distance for chromatin effects.
        nfeatures: number of chromatin/epigenomic features.
        batchSize: batch size when computing ExPecto predictions.

    Returns:
        numpy array of size num_variants x num_models. Each value represents
        predicted log fold-change
    """
    snpdists = snpdists * ((snpstrands == '+') * 2 - 1)
    Xreducedall_diffs = [np.vstack([
    np.exp(-0.01 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) <= 0),
    np.exp(-0.02 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) <= 0),
    np.exp(-0.05 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) <= 0),
    np.exp(-0.1 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) <= 0),
    np.exp(-0.2 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) <= 0),
    np.exp(-0.01 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) >= 0),
    np.exp(-0.02 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) >= 0),
    np.exp(-0.05 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) >= 0),
    np.exp(-0.1 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) >= 0),
    np.exp(-0.2 * np.floor(np.abs((snpdists + dist * ((snpstrands == '+') * 2 - 1)
           ) / 200.0))) * ((snpdists + dist * ((snpstrands == '+') * 2 - 1)) >= 0)
     ]).T for dist in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200))]
    n_snps = len(snpdists)
    effect = np.zeros((n_snps, len(all_models)))


    for i in range(int( (n_snps - 1) / batchSize) + 1):
        print("Processing " + str(i) + "th batch of "+str(batchSize))
        # compute gene expression change with models
        diff = reduce(lambda x, y: x + y, [np.tile(np.asarray(snpeffects[j][i * batchSize:(i + 1) * batchSize, :]), 10)
                                 * np.repeat(Xreducedall_diffs[j][i * batchSize:(i + 1) * batchSize, :], nfeatures, axis=1) for j in range(len(Xreducedall_diffs))])
        if old_format:
            # backward compatibility
            diff = np.concatenate([np.zeros((diff.shape[0], 10, 1)), diff.reshape(
                (-1, 10, 2002))], axis=2).reshape((-1, 20030))
        dtest_ref = xgb.DMatrix(diff * 0)
        dtest_alt = xgb.DMatrix(diff)

        for j in range(len(all_models)):
            effect[i * batchSize:(i + 1) * batchSize, j] = all_models[j].predict(dtest_alt) - \
                            all_models[j].predict(dtest_ref)

    return effect

#load resources
modelList = pd.read_csv(args.modelList,sep='\t',header=0)
models = []
for file in modelList['ModelName']:
        bst = xgb.Booster({'nthread': args.threads})
        bst.load_model(file.strip())
        models.append(bst)

# backward compatibility with earlier model format
if len(bst.get_dump()[0].split('\n')) == 20034:
    old_format = True
else:
    old_format = False


#load input data
maxshift = int(args.maxshift)
snpEffects = []
for shift in [str(n) for n in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200))]:
    h5f = h5py.File(args.snpEffectFilePattern.replace(
        'SHIFT', shift), 'r')['/pred']

    if args.splitFlag:
        index_start = int((args.splitIndex - 1) *
                          np.ceil(float(h5f.shape[0] / 2) / args.splitFold))
        index_end = int(np.minimum(
            (args.splitIndex) * np.ceil(float(h5f.shape[0] / 2) / args.splitFold), (h5f.shape[0] / 2)))
    else:
        index_start = 0
        index_end = int(h5f.shape[0] / 2)

    snp_temp = (np.asarray(h5f[index_start:index_end,:])+ np.asarray(h5f[index_start+int(h5f.shape[0]/2):index_end+int(h5f.shape[0]/2),:]))/2.0
    snpEffects.append(snp_temp)


coor = pd.read_csv(args.coorFile,sep='\t',header=None,comment='#')
coor = coor.iloc[index_start:index_end,:]

#Fetch the distance to TSS information
gene = pd.read_csv(args.geneFile,sep='\t',header=None,comment='#')
geneinds = pd.match(coor.iloc[:,0].map(str).str.replace('chr','')+' '+coor.iloc[:,1].map(str),
            gene.iloc[:,0].map(str).str.replace('chr','')+' '+gene.iloc[:,2].map(str))
if np.any(geneinds==-1):
    raise ValueError("Gene association file does not match the vcf file.")
if args.fixeddist == 0:
    dist = - np.asarray(gene.iloc[geneinds,-1])
else:
    dist = args.fixeddist
genename = np.asarray(gene.iloc[geneinds,-2])
strand= np.asarray(gene.iloc[geneinds,-3])

#comptue expression effects
snpExpEffects = compute_effects(snpEffects, \
                                dist, strand,\
                                models, maxshift=maxshift, nfeatures=args.nfeatures,
                                batchSize = args.batchSize, old_format = old_format)
#write output
snpExpEffects_df = coor
snpExpEffects_df['dist'] = dist
snpExpEffects_df['gene'] = genename
snpExpEffects_df['strand'] = strand
snpExpEffects_df=pd.concat([snpExpEffects_df.reset_index(),pd.DataFrame(snpExpEffects, columns = modelList.iloc[:,1])],axis=1,ignore_index =False)
snpExpEffects_df.to_csv(args.output, header = True)
