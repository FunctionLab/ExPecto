# -*- coding: utf-8 -*-
"""Compute chromatin representation of variants (required by predict.py).

This script takes a vcf file, compute the effect of the variant both at the
variant position and at nearby positions, and output the effects as the
representation that can be used by predict.py.

Example:
        $ python chromatin.py ./example/example.vcf

"""
import argparse
import math
import pyfasta
import torch
from torch.utils.serialization import load_lua
import numpy as np
import pandas as pd
import h5py

parser = argparse.ArgumentParser(description='Predict variant chromatin effects')
parser.add_argument('inputfile', type=str, help='Input file in vcf format')
parser.add_argument('--maxshift', action="store",
                    dest="maxshift", type=int, default=800,
                    help='Maximum shift distance for computing nearby effects')
parser.add_argument('--inputsize', action="store", dest="inputsize", type=int,
                    default=2000, help="The input sequence window size for neural network")
parser.add_argument('--batchsize', action="store", dest="batchsize",
                    type=int, default=32, help="Batch size for neural network predictions.")
parser.add_argument('--cuda', action='store_true')
args = parser.parse_args()

genome = pyfasta.Fasta('./resources/hg19.fa')
model = load_lua('./resources/deepsea.beluga.2002.cpu')
model.evaluate()
if args.cuda:
    model.cuda()

CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']


inputfile = args.inputfile
maxshift = args.maxshift
inputsize = args.inputsize
batchSize = args.batchsize
windowsize = inputsize + 100


def encodeSeqs(seqs, inputsize=2000):
    """Convert sequences to 0-1 encoding and truncate to the input size.
    The output concatenates the forward and reverse complement sequence
    encodings.

    Args:
        seqs: list of sequences (e.g. produced by fetchSeqs)
        inputsize: the number of basepairs to encode in the output

    Returns:
        numpy array of dimension: (2 x number of sequence) x 4 x inputsize

    2 x number of sequence because of the concatenation of forward and reverse
    complement sequences.
    """
    seqsnp = np.zeros((len(seqs), 4, inputsize), np.bool_)

    mydict = {'A': np.asarray([1, 0, 0, 0]), 'G': np.asarray([0, 1, 0, 0]),
            'C': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
            'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
            'a': np.asarray([1, 0, 0, 0]), 'g': np.asarray([0, 1, 0, 0]),
            'c': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
            'n': np.asarray([0, 0, 0, 0]), '-': np.asarray([0, 0, 0, 0])}

    n = 0
    for line in seqs:
        cline = line[int(math.floor(((len(line) - inputsize) / 2.0))):int(math.floor(len(line) - (len(line) - inputsize) / 2.0))]
        for i, c in enumerate(cline):
            seqsnp[n, :, i] = mydict[c]
        n = n + 1

    # get the complementary sequences
    dataflip = seqsnp[:, ::-1, ::-1]
    seqsnp = np.concatenate([seqsnp, dataflip], axis=0)
    return seqsnp


def fetchSeqs(chr, pos, ref, alt, shift=0, inputsize=2000):
    """Fetches sequences from the genome.

    Retrieves sequences centered at the given position with the given inputsize.
    Returns both reference and alternative allele sequences . An additional 100bp
    is retrived to accommodate indels.

    Args:
        chr: the chromosome name that must matches one of the names in CHRS.
        pos: chromosome coordinate (1-based).
        ref: the reference allele.
        alt: the alternative allele.
        shift: retrived sequence center position - variant position.
        inputsize: the targeted sequence length (inputsize+100bp is retrived for
                reference allele).

    Returns:
        A string that contains sequence with the reference allele,
        A string that contains sequence with the alternative allele,
        A boolean variable that tells whether the reference allele matches the
        reference genome

        The third variable is returned for diagnostic purpose. Generally it is
        good practice to check whether the proportion of reference allele
        matches is as expected.

    """
    windowsize = inputsize + 100
    mutpos = int(windowsize / 2 - 1 - shift)
    # return string: ref sequence, string: alt sequence, Bool: whether ref allele matches with reference genome
    seq = genome.sequence({'chr': chr, 'start': pos + shift -
                           int(windowsize / 2 - 1), 'stop': pos + shift + int(windowsize / 2)})
    return seq[:mutpos] + ref + seq[(mutpos + len(ref)):], seq[:mutpos] + alt + seq[(mutpos + len(ref)):], seq[mutpos:(mutpos + len(ref))].upper() == ref.upper()


vcf = pd.read_csv(inputfile, sep='\t', header=None, comment='#')

# standardize
vcf.iloc[:, 0] = 'chr' + vcf.iloc[:, 0].map(str).str.replace('chr', '')
vcf = vcf[vcf.iloc[:, 0].isin(CHRS)]

for shift in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200)):
    refseqs = []
    altseqs = []
    ref_matched_bools = []
    for i in range(vcf.shape[0]):
        refseq, altseq, ref_matched_bool = fetchSeqs(
            vcf[0][i], vcf[1][i], vcf[3][i], vcf[4][i], shift=shift)
        refseqs.append(refseq)
        altseqs.append(altseq)
        ref_matched_bools.append(ref_matched_bool)

    if shift == 0:
        # only need to be checked once
        print("Number of variants with reference allele matched with reference genome:")
        print(np.sum(ref_matched_bools))
        print("Number of input variants:")
        print(len(ref_matched_bools))

    ref_encoded = encodeSeqs(refseqs).astype(np.float32)
    alt_encoded = encodeSeqs(altseqs).astype(np.float32)
    #print(ref_encoded.shape)
    #print(alt_encoded.shape)

    ref_preds = []
    for i in range(int(1 + ref_encoded.shape[0] / batchSize)):
        input = torch.from_numpy(ref_encoded[int(i*batchSize):int((i+1)*batchSize),:,:]).unsqueeze(3)
        if args.cuda:
            input = input.cuda()
        ref_preds.append(model.forward(input).cpu().numpy().copy())
    ref_preds = np.vstack(ref_preds)

    alt_preds = []
    for i in range(int(1 + alt_encoded.shape[0] / batchSize)):
        input = torch.from_numpy(alt_encoded[int(i*batchSize):int((i+1)*batchSize),:,:]).unsqueeze(3)
        if args.cuda:
            input = input.cuda()
        alt_preds.append(model.forward(input).cpu().numpy().copy())
    alt_preds = np.vstack(alt_preds)

    #ref_preds = np.vstack(map(lambda x: model.forward(x).numpy(),
    #        np.array_split(ref_encoded, int(1 + len(refseqs) / batchSize))))
    #alt_encoded = torch.from_numpy(encodeSeqs(altseqs).astype(np.float32)).unsqueeze(3)
    #alt_preds = np.vstack(map(lambda x: model.forward(x).numpy(),
    #        np.array_split(alt_encoded, int(1 + len(altseqs) / batchSize))))
    diff = alt_preds - ref_preds
    f = h5py.File(inputfile + '.shift_' + str(shift) + '.diff.h5', 'w')
    f.create_dataset('pred', data=diff)
    f.close()
