import pysam,sys
import numpy as np
import pandas as pd
from src.utils.clust import readClusterFile
from src.utils.motif import getCounts,clusterStats

DEFAULTNAME='./motif.summary.csv'
DEFAULTLF  ='totalBp'

def main(parser):
    args = parser.parse_args()

    clusters = readClusterFile(args.clusterFile)
    if np.any(clusters.start.isnull()) or args.full:
        extract = False
        if not args.full:
            print('WARNING: using full read sequence for motif counting')
    else:
        extract = True
    seqs        = seqGen(args.inBAM,clusters,extract)
    motifs      = args.motifs.split(',')
    motifCounts = getCounts(seqs,motifs,lengthField=DEFAULTLF)
    clusterIdx  = motifCounts.index.map(clusters.cluster.to_dict())
    results     = clusterStats(motifCounts,clusterIdx,
                               motifs+[DEFAULTLF],
                               randomSeed=args.seed)
    
    results.to_csv(args.outName, float_format='%.1f')

    return results

def seqGen(bamfile,clusterDf,extract=False):
    bam = pysam.AlignmentFile(bamfile)
    parse = (lambda seq,s,e:seq[s:e]) if extract else (lambda seq,s,e:seq)
    for rec in bam:
        if rec.query_name in clusterDf.index:
            s,e = clusterDf.loc[rec.query_name,['start','stop']]
            yield rec.query_name,parse(rec.query_sequence,s,e)


class MotifCounter_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='motifCounter.py', description='summarize motif counts in clustered reads')
    parser.add_argument('clusterFile', metavar='clusterFile', type=str,
                    help='cluster.txt output from clustering script.')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='aligned BAM containing reads id\'d in clusterFile')
    parser.add_argument('-m,--motifs', dest='motifs', type=str, default=None,required=True,
                    help='comma-separated list of motifs to count')
    parser.add_argument('-o,--outName', dest='outName', type=str, default=DEFAULTNAME,
                    help=f'Output filename. Default {DEFAULTNAME}')
    parser.add_argument('-s,--seed', dest='seed', type=int, default=42,
                    help='Seed for resampling confidence interval.  Default 42')
    parser.add_argument('-f,--full', dest='full', action='store_true', default=False,
                    help='Count motifs for full read.  Default False (try to extract subsequence if defined in cluster.txt)')

    try:
        main(parser)
    except MotifCounter_Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1)
