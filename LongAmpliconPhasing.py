import pysam
from src.phase.phaser import Phaser,Phaser_Error
from src.phase.split import Splitter_Error
from src.phase.utils import getFileType,PhaseUtils_Error
from src.utils.bam import addHPtag

DEFAULTPREFIX   = './longamp'
DEFAULTMINREADS = 5
DEFAULTMINLENGTH= 50
DEFAULTMAXLENGTH= 50000
DEFAULTMINFRAC  = 0.1
DEFAULTMINSIG   = 0.05
DEFAULTKMER     = 55

MAXCOMPRESSSIZE = 500000
#DEFAULTMINSPAN  = 0.95
#DEFAULTMINCOV   = 0

def main(parser):
    args = parser.parse_args()

    #Make sure the inputs are logically consistent with outputs
    ftype = getFileType(args.inFile)
    if args.reference is None:
        if args.variants:
            raise LongAmpliconPhasing_Error('Cannot export variant info without reference')
        if args.maxHP != 0 and args.method == 'align':
            print(f'WARNING: compressed reads will be realigned to {args.template} input read')
    else:
        if args.maxHP !=0 and args.method=='align':
            refLens = [len(rec.sequence) for rec in pysam.FastxFile(args.reference)]
            if sum(refLens) > MAXCOMPRESSSIZE:
                raise LongAmpliconPhasing_Error('Reference size too large for compressed realignment.  Try subsetting reference')
            if len(refLens) > 1:
                raise LongAmpliconPhasing_Error('Compression of multi-record references is not yet supported')
        if ftype != 'bam':
            if args.region is not None:
                raise LongAmpliconPhasing_Error('Region option can only be used with BAM input')
            if args.method == 'debruijn':
                print('WARNING: Reads should be oriented before running')
        else:
            if args.maxHP !=0 and args.method=='align':
                print('WARNING: Reads will be realigned after compression')
         
    #########

    if (args.method == 'align'  and args.maxHP != 0) \
      or (ftype == 'fastq' and (args.method == 'align' or args.variants)):
        #need to remap
        from src.phase.utils import PilerUpper
        refseq = next(pysam.FastxFile(args.reference)).sequence if args.reference else None
        varDf  = PilerUpper(args.inFile,
                            region=args.region,
                            refSeq=refseq,
                            method=args.template,
                            minLength=args.minLength,
                            maxLength=args.maxLength,
                            maxHP=args.maxHP).varDf
    else:
        varDf = None

    #load splitter class
    if args.method == 'align':
        from src.phase.split import VariantGrouper
        splitter = VariantGrouper(args.inFile,
                                  args.reference,
                                  region=args.region,
                                  truncate=args.truncate,
                                  minFrac=args.minFrac,
                                  minReads=args.minReads,
                                  minSignal=args.minSignal,
                                  flagFilter=args.flagFilter,
                                  aggressive=args.aggressive,
                                  vTable=varDf)
    elif args.method == 'debruijn':
        from src.phase.split import sparseDBG
        splitter = sparseDBG(args.k,
                             collapse=args.maxHP,
                             ignoreEnds=args.ignore,
                             minReads=args.minReads,
                             minFrac=args.minFrac)
        splitter.loadReads(args.inFile,
                           region=args.region,
                           minLength=args.minLength,
                           maxLength=args.maxLength)
    else:
        raise LongAmpliconPhasing_Error('unknown method')

    phaser = Phaser(splitter,
                    args.sampleName,
                    aggressive=args.aggressive)
    phaser.run()
    
    #export stuff
    s = '' if args.prefix.endswith('/') else '.'
    #bam
    if not args.noBam and ftype == 'bam':
        print('Writing BAM(s)')
        #hacky
        args.inBAM = args.inFile
        outBAM = f'{args.prefix}{s}hptagged.bam' 
        addHPtag(args,outBAM,phaser.clusterMap)
    #simple decision tree
    print('Writing Outputs')
    with open(f'{args.prefix}{s}clusterSplits.txt','w') as ofile:
        ofile.write(f'{str(phaser)}\n')
    #cluster file/readnames
    with open(f'{args.prefix}{s}clusters.txt','w') as ofile:
        for node,clust in phaser.node2cluster.items():
            vclust = phaser.vTree[node]
            ofile.write(f'>{vclust.clusterName()}\n')
            ofile.write('\n'.join(vclust.reads) + '\n')
    #check if anything was produced, else exit
    hasResults = True
    if len(phaser.node2cluster) == 0:
        print(f'WARNING: No results produced for sample {phaser.sampleName}.  Nreads: {phaser.splitter.nReads}')
        hasResults = False
    #variants/consensus
    if args.variants:
        from src.phase.caller import VariantCaller,FIGFORMAT 
        if args.method == 'debruijn' or args.maxHP != 0:
            #need to build the variant splitter
            from src.phase.split import VariantGrouper
            vsplitter = VariantGrouper(args.inFile,
                                       args.reference,
                                       region=args.region,
                                       truncate=args.truncate,
                                       minFrac=args.minFrac,
                                       minReads=args.minReads,
                                       minSignal=args.minSignal,
                                       flagFilter=args.flagFilter,
                                       aggressive=args.aggressive,
                                       vTable=None)
        else:
            vsplitter = splitter

        endpoints = vsplitter.vTable.groupby(phaser.clusterMap).apply(vsplitter.getEndpoints)
        caller    = VariantCaller(vsplitter.sigVar,
                                  vsplitter.minCount,
                                  args.reference,
                                  phaser.clusterMap,
                                  endpoints)
        #summary table
        name = f'{args.prefix}{s}alleleClusterSummary.csv'
        caller.alleleFrequency(caller.clusterMap).to_csv(name)
        name = f'{args.prefix}{s}sampleVariantSummary.csv'
        caller.alleleFrequency(caller.variantGroupMap)\
              .reindex(list(set(caller.variantGroupMap.values())))\
              .to_csv(name)
        
        #total allele fractions (draft)
        name = f'{args.prefix}{s}variantFraction.csv'
        caller.variantFractions.to_csv(name)
        #draft consensus
        name = f'{args.prefix}{s}draft.consensus.fasta'
        with open(name,'w') as ofasta:
            for clust,seq in caller.draftConsensus().items():
                clust = clust if clust!=-1 else 'Noise'
                ofasta.write(f'>cluster{clust}.plurality\n')
                ofasta.write(seq + '\n')
        #entropy
        name = f'{args.prefix}{s}entropy.csv'
        caller.entropy.to_csv(name)
        if args.entropyPlot:
            name = f'{args.prefix}{s}entropy.{FIGFORMAT}'
            caller.entropyPlot(name)

    return phaser
                                 
class LongAmpliconPhasing_Error(Exception):
    pass


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='LongAmpliconPhasing.py', description='Recursively separate amplicons by shared variants')
    parser.set_defaults(prog=parser.prog)
    parser.add_argument('inFile', metavar='inFile', type=str,
                    help='Fasta,Fastq, or BAM containing CCS reads')
    parser.add_argument('sampleName', metavar='sampleName', type=str,
                    help='sample name')
    parser.add_argument('--reference', dest='reference', type=str, default=None,
                    help='reference fasta used for alignment.  Required for draft consensus and variant counts')
    parser.add_argument('--region', dest='region', type=str,default=None,
                    help='reference subregion. Format [chr]:[start]-[stop].  Default None.')
    parser.add_argument('-t','--truncate', dest='truncate', action='store_true', default=False,
                    help='Truncate to columns only within region (pysam pileup).  Default False')
    parser.add_argument('-k','--kmer', dest='k', type=int, default=DEFAULTKMER,
                    help=f'Kmer size for debruijn method.  Ignored for align method. Default {DEFAULTKMER}')
    parser.add_argument('-i','--ignore', dest='ignore', type=int, default=0,
                    help=f'Ignore first and last N bases when building debruijn graph. Not used for align method. Default 0')
    parser.add_argument('-p','--prefix', dest='prefix', type=str, default='',
                    help=f'Output prefix. Default {DEFAULTPREFIX}')
    parser.add_argument('-d','--drop', dest='drop', action='store_true', default=False,
                    help='Drop noise/other reads.  Default False')
    parser.add_argument('-s','--splitBam', dest='splitBam', action='store_true', default=False,
                    help='Generate one BAM per cluster.  Default False (HP tag single bam)')
    parser.add_argument('-x','--noBam', dest='noBam', action='store_true', default=False,
                    help='Generate BAM outputs from BAM inputs.  Default export BAM')
    parser.add_argument('-v','--variants', dest='variants', action='store_true', default=False,
                    help='Generate variant and draft consensus outputs. Only valid with reference. Default do not export variants/consensus')
    parser.add_argument('-r','--minReads', dest='minReads', type=int, default=DEFAULTMINREADS,
                    help=f'Minimum abs reads sharing variant to split.  Default {DEFAULTMINREADS}')
    parser.add_argument('-l','--minLength', dest='minLength', type=int, default=DEFAULTMINLENGTH,
                    help=f'Minimum length read.  Default {DEFAULTMINLENGTH}')
    parser.add_argument('-L','--maxLength', dest='maxLength', type=int, default=DEFAULTMAXLENGTH,
                    help=f'Maximum length read.  Default {DEFAULTMAXLENGTH}')
    parser.add_argument('-f','--minFrac', dest='minFrac', type=float, default=DEFAULTMINFRAC,
                    help=f'Minimum fraction of reads sharing variant to split.  Default {DEFAULTMINFRAC}')
    parser.add_argument('-g','--minSignal', dest='minSignal', type=float, default=DEFAULTMINSIG,
                    help=f'Minimum fraction of reads with non-ref calls for consideration.  Default {DEFAULTMINSIG}')
    parser.add_argument('-H','--maxHP', dest='maxHP', type=int, default=0,
                    help='Max length of homopolymer seqs.  Set to 0 to turn of hp collapsing.  Default 0')
    parser.add_argument('-a','--aggressive', dest='aggressive', action='store_true', default=False,
                    help='More aggressively split groups (experimental).  Default False')
    parser.add_argument('-e','--entropyPlot', dest='entropyPlot', action='store_true', default=False,
                    help='Generate heatmap of entropy for clusters.  Default False')
    parser.add_argument('-F','--flagFilter', dest='flagFilter', type=hex, default=0x900,
                    help='Filter to pass to pysam.pileup. Default 0x900')
    parser.add_argument('--template', dest='template', choices=['first','median'], default='median',
                    help='Method for choosing reference template from inputs (if no reference passed). Default median')
    parser.add_argument('-m','--method', dest='method', choices=['align','debruijn'], default='align',
                    help='Splitting method.  If align and maxHP != 0, reads will be realigned after compression for clustering; Output variants are from non-compressed pileup. Default align')
#    parser.add_argument('-c','--minCov', dest='minCov', type=int, default=DEFAULTMINCOV,
#                    help=f'Minimum read coverage to include position.  Default {DEFAULTMINCOV}')
#    parser.add_argument('-S','--minSpan', dest='minSpan', type=float, default=DEFAULTMINSPAN,
#                    help=f'NOT IMPLEMENTED. Minimum fraction of positions covered to include read.  Default {DEFAULTMINSPAN}')

    try:
        main(parser)
    except (Phaser_Error,Splitter_Error,
            PhaseUtils_Error,LongAmpliconPhasing_Error) as e:
        print(f'ERROR: {e}')
        sys.exit(1)
