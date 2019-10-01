#! ~/anaconda3/bin/python

import os,sys

from src.model.kmer   import *
from src.model.models import MODELS
from src.utils.bam import addHPtag

DEFAULTMODEL    = 'dbscan'
DEFAULTKMER     = 11
DEFAULTNORM     = 'l1'
DEFAULTMINREADS = 5
DEFAULTCOMP     = 2
DEFAULTEPS      = 0.0
DEFAULTSIMP     = 0.0
DEFAULTPREFIX   = './clustered'
PALETTE         = 'husl'
CLUSTNAME       = 'cluster{0}_numreads{1}'.format #args: int,int


def main(parser):
    args = parser.parse_args()

    #load dataframe with samples(row) by kmer counts (cols)
    data = loadKmers(args.inBAM,args.minQV,args.kmer,
                     collapse   =args.hpCollapse,
                     region     =args.region,
                     minimizer  =args.minimizer,
                     ignoreEnds =args.ignoreEnds,
                     whitelist  =args.whitelist,
                     flanks     =args.flanks,
                     simpson    =args.simpson,
                     norm       =args.normalize,
                     components =args.components,
                     agg        =args.agg)

    if args.testPlot:
        #plot kdist and exit
        print("Plotting distance to k-neighbors")
        f = plotEPS(data,args.minReads)
        f.savefig('%s.eps_estimator.png' % args.prefix)
        return data
    
    #Clustering
    cluster     = MODELS[args.model](args)
    print("Clustering with {m}\n{p}".format(m=args.model,
                                            p=printParams(cluster)))
    result      = cluster.fit(data)
    clusterIdx  = result.labels_

    #TODO
    #cluster size and warning if too much noise as frac of total

    #write cluster file
    print("Writing clusters")
    with open('%s.clusters.txt' % args.prefix, 'w') as namefile:
        for cluster,reads in data.groupby(clusterIdx):
            nreads = len(reads)
            name = 'Noise_numreads%i'%nreads if cluster==-1 else CLUSTNAME(cluster,nreads)
            namefile.write('>%s\n' % name)
            namefile.write('\n'.join(reads.index) + '\n')
    
    if not args.noBam:
        print("Adding HP tag to bam")
        clusterMap = dict(zip(data.index,clusterIdx))
        outBam     = '%s.hptagged.bam' % args.prefix
        addHPtag(args.inBAM,outBam,clusterMap,dropNoClust=args.drop,splitBam=args.splitBam)

    if args.plotReads:
        fig = plotReads(data,clusterIdx)
        fig.tight_layout()
        fig.savefig('%s.clusters.png' % args.prefix)
    
    return None

def printParams(model):
    return '\n'.join(['\t' + '='.join(map(str,v)) for v in model.defaults.items()])

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='ClusterAmplicons.py', description='Clustering by kmer counts')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='input BAM of CCS alignments')
    parser.add_argument('-j,--njobs', dest='njobs', type=int, default=None,
                    help='j parallel jobs (only for OPTICS). Default 1')
    kmer = parser.add_argument_group('kmers')
    kmer.add_argument('-k,--kmer', dest='kmer', type=int, default=DEFAULTKMER,
                    help='kmer size for clustering. Default %i'%DEFAULTKMER)
    kmer.add_argument('-z,--minimizer', dest='minimizer', type=int, default=0,
                    help='group kmers by minimizer of length z. Default 0 (no minimizer)')
    kmer.add_argument('-H,--noHPcollapse', dest='hpCollapse', action='store_false', default=True,
                    help='do not compress homopolymers.  Default collapse HP')
    clust = parser.add_argument_group('cluster')
    clust.add_argument('-M,--model', dest='model', type=str, choices=MODELS.keys(), default=DEFAULTMODEL,
                    help='clustering model. Default %s'%DEFAULTMODEL)
    clust.add_argument('-a,--agg', dest='agg', type=str, choices=['pca','featagg'],default='pca',
                    help='Feature reduction method. Default pca')
    clust.add_argument('-c,--components', dest='components', type=int, default=DEFAULTCOMP,
                    help='Use first c components of PCA/FeatAgg for clustering. Set to 0 for no reduction. Default %i'%DEFAULTCOMP)
    clust.add_argument('-e,--eps', dest='eps', type=float, default=None,
                    help='eps cluster tolerance. Default None')
    clust.add_argument('-m,--minReads', dest='minReads', type=int, default=DEFAULTMINREADS,
                    help='Minimum reads to be a cluster. Default %i'%DEFAULTMINREADS)
    clust.add_argument('-n,--normalize', dest='normalize', type=str, choices=['l1','l2','none'], default=DEFAULTNORM,
                    help='normalization of kmer counts.  Default %s'%DEFAULTNORM)
    clust.add_argument('-i,--ignoreEnds', dest='ignoreEnds', type=int, default=0,
                    help='ignore i bases at ends of amplicons for clustering.  Default 0')
    filt = parser.add_argument_group('filter')
    filt.add_argument('-r,--region', dest='region', type=str, default=None,
                    help='Target region for selection of reads, format \'[chr]:[start]-[stop]\'.  Example \'4:3076604-3076660\'. \nDefault all reads (no region)')
    filt.add_argument('-q,--minQV', dest='minQV', type=float, default=0.99,
                    help='Minimum quality [0-1] to use for clustering. Default 0.99')
    filt.add_argument('-s,--simpsonDominance', dest='simpson', type=float, default=DEFAULTSIMP,
                    help='Dominance filter for kmers.  Remove kmers with > s (dominance). Default %.2f (no filter)'%DEFAULTSIMP)
    filt.add_argument('-w,--whitelist', dest='whitelist', type=str, default=None,
                    help='whitelist of reads to cluster. Default None')
    filt.add_argument('-f,--flanks', dest='flanks', type=str, default=None,
                    help='fasta of flanking/primer sequence. Reads not mapping to both will be filtered. Default None')
    out = parser.add_argument_group('output')
    out.add_argument('-p,--prefix', dest='prefix', type=str, default=DEFAULTPREFIX,
                    help='Output prefix. Default %s'%DEFAULTPREFIX)
    out.add_argument('-S,--splitBam', dest='splitBam', action='store_true',
                    help='split clusters into separate bams (noise and no-cluster dropped). Default one bam')
    out.add_argument('-x,--noBam', dest='noBam', action='store_true',
                    help='Do not export HP-tagged bam of clustered reads')
    out.add_argument('-d,--drop', dest='drop', action='store_true',
                    help='Drop reads with no cluster in output bam.  Default keep all reads.')
    out.add_argument('-t,--testPlot', dest='testPlot', action='store_true',
                    help='Plot reads vs dist to nearest m-neighbors without clustering')
    out.add_argument('-g,--plotReads', dest='plotReads', action='store_true',
                    help='Plot first 2 axes of PCA for each read.  Default no plot generated')

    try:
        main(parser)
    except Clustering_Exception as e:
        print('ERROR: %s' % e)
        sys.exit(1)
