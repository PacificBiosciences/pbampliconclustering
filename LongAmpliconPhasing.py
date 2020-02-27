from src.phase.variant import VariantGrouper,VariantPartition_Error
from src.utils.bam import addHPtag

DEFAULTPREFIX  = './longamp'
DEFAULTMINREADS= 10
DEFAULTMINFRAC = 0.1
DEFAULTMINSIG  = 0.05
DEFAULTMINSPAN = 0.95
DEFAULTMINCOV  = 10
DEFAULTMINBQ   = 0

def main(parser):
    args = parser.parse_args()

    readGrouper = VariantGrouper(args.inBAM,args.reference,args.sampleName,
                                 region=args.region,
                                 truncate=args.truncate,
                                 minCov=args.minCov,
                                 minFrac=args.minFrac,
                                 minReads=args.minReads,
                                 minSpan=args.minSpan,
                                 minSignal=args.minSignal,
                                 minBaseQ=args.minBaseQ,
                                 flagFilter=args.flagFilter)
    readGrouper.run()

    #export stuff
    s = '' if args.prefix.endswith('/') else '.'
    #bam
    if not args.noBam:
        print('Writing BAM(s)')
        outBAM = f'{args.prefix}{s}hptagged.bam' 
        addHPtag(args,outBAM,readGrouper.clusterMap)
    #simple decis:ion tree
    print('Writing Outputs')
    with open(f'{args.prefix}{s}clusterSplits.txt','w') as ofile:
        ofile.write(f'{str(readGrouper)}\n')
    #cluster file/readnames
    with open(f'{args.prefix}{s}clusters.txt','w') as ofile:
        for node,clust in readGrouper.node2cluster.items():
            vclust = readGrouper.vTree[node]
            ofile.write(f'>{vclust.clusterName()}\n')
            ofile.write('\n'.join(vclust.reads) + '\n')
    #summary table
    name    = f'{args.prefix}{s}summary.csv'
    readGrouper.summary().to_csv(name)

    return readGrouper
                                 
class LongAmpliconPhasing_Error(Exception):
    pass


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='LongAmpliconPhasing.py', description='Recursively separate amplicons by shared variants')
    parser.set_defaults(prog=parser.prog)
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='aligned BAM containing CCS reads')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='reference fasta used for aligning')
    parser.add_argument('sampleName', metavar='sampleName', type=str,
                    help='sample name')
    parser.add_argument('--region', dest='region', type=str,default=None,
                    help='reference subregion. Format [chr]:[start]-[stop].  Default None.')
    parser.add_argument('-t','--truncate', dest='truncate', action='store_true', default=False,
                    help='Truncate to columns only within region (pysam pileup).  Default False')
    parser.add_argument('-p','--prefix', dest='prefix', type=str, default='',
                    help=f'Output prefix. Default {DEFAULTPREFIX}')
    parser.add_argument('-d','--drop', dest='drop', action='store_true', default=False,
                    help='Drop noise/other reads.  Default False')
    parser.add_argument('-s','--splitBam', dest='splitBam', action='store_true', default=False,
                    help='Generate one BAM per cluster.  Default False (HP tag single bam)')
    parser.add_argument('-x','--noBam', dest='noBam', action='store_true', default=False,
                    help='Generate BAM outputs.  Default export BAM')
    parser.add_argument('-r','--minReads', dest='minReads', type=int, default=DEFAULTMINREADS,
                    help=f'Minimum abs reads sharing variant to split.  Default {DEFAULTMINREADS}')
    parser.add_argument('-f','--minFrac', dest='minFrac', type=float, default=DEFAULTMINFRAC,
                    help=f'Minimum fraction of reads sharing variant to split.  Default {DEFAULTMINFRAC}')
    parser.add_argument('-g','--minSignal', dest='minSignal', type=float, default=DEFAULTMINSIG,
                    help=f'Minimum fraction of reads with non-ref calls for consideration.  Default {DEFAULTMINSIG}')
    parser.add_argument('-c','--minCov', dest='minCov', type=int, default=DEFAULTMINCOV,
                    help=f'Minimum read coverage to include position.  Default {DEFAULTMINCOV}')
    parser.add_argument('-S','--minSpan', dest='minSpan', type=float, default=DEFAULTMINSPAN,
                    help=f'NOT IMPLEMENTED. Minimum fraction of positions covered to include read.  Default {DEFAULTMINSPAN}')
    parser.add_argument('-q','--minBaseQ', dest='minBaseQ', type=int, default=DEFAULTMINBQ,
                    help=f'Minimum base quality (phred) to include. (Recommend 0) Default {DEFAULTMINBQ}')
    parser.add_argument('-F','--flagFilter', dest='flagFilter', type=hex, default=0x900,
                    help='Filter to pass to pysam.pileup. Default 0x900')

    try:
        main(parser)
    except (VariantPartition_Error,LongAmpliconPhasing_Error) as e:
        print(f'ERROR: {e}')
        sys.exit(1)
