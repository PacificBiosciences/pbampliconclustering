import pysam,sys,traceback
import pandas as pd
from src.phase.phaser import Phaser,Phaser_Error
from src.phase.split import Splitter_Error
from src.phase.utils import getFileType,writeSimpleBED,writeRegionBam,summary,PhaseUtils_Error,hpCollapse
from src.phase.caller import Caller_Error
from src.utils.bam import addHPtag,getCoordinates,exportFastq
from src.utils.logging import getLogger

#DEFAULTPREFIX   = './longamp'
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
    args   = parser.parse_args()
    #exporting prefix
    pref   = args.prefix if args.prefix else args.sampleName
    s      = '' if pref.endswith('/') else '.'
    prefix = f'{pref}{s}'

    #logging
    log = getLogger('lap',f'{prefix}laphase.log',stdout=args.verbose)
    log.debug(f'Command: {" ".join(sys.argv)}')

    #Make sure the inputs are logically consistent with outputs
    ftype = getFileType(args.inFile)
    if args.minSignal > args.minFrac:
        log.warning(f'minSignal ({args.minSignal}) > minFraction ({args.minFrac}).  Some groups may be missed.')
    if args.reference is None:
        #if args.method == 'align':
        #    raise LongAmpliconPhasing_Error('Need reference for this method')
        if args.variants:
            raise LongAmpliconPhasing_Error('Cannot export variant info without reference')
        if args.maxHP != 0 and args.method == 'align':
            log.warning(f'compressed reads will be realigned to {args.template} input read')
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
                log.warning('Reads should be oriented before running')
        else:
            if args.maxHP !=0 and args.method=='align':
                log.warning('Reads will be realigned after compression')
         
    #########

    if ((args.method == 'align' or args.method == 'cluster')  and args.maxHP != 0) \
      or (ftype == 'fastq' and (args.method == 'align' or args.variants or args.method == 'cluster')):
        #need to remap
        from src.phase.utils import PilerUpper
        refseq = next(pysam.FastxFile(args.reference)).sequence if args.reference else None
        pileup  = PilerUpper(args.inFile,
                             region=args.region,
                             refSeq=refseq,
                             method=args.template,
                             minLength=args.minLength,
                             maxLength=args.maxLength,
                             maxHP=args.maxHP,
                             log=log,
                             nproc=args.nproc)
                             #multifunc=getRow)
        varDf   = pileup.varDf
        counts = pileup.recGen.counter
        secNsup = sum(counts.get(k,0) for k in ['secondary','supplementary'])
        filtlen = sum(counts.get(k,0) for k in counts.keys() if k.startswith('long') or k.startswith('short'))
        prim    = pileup.recGen.counter['pass']
        stats   = {'total alignments'  : prim + secNsup + filtlen,
                   'primary alignments': prim,
                   'pileup alignments' : len(varDf)}
    else:
        varDf = None
        stats = {}

    #load splitter class
    if args.method != 'debruijn':
        if args.method == 'align':
            from src.phase.split import VariantGrouper
            sclass = VariantGrouper
        elif args.method == 'cluster':
            from src.phase.split import VariantSubCluster
            sclass = VariantSubCluster
        else: 
            raise LongAmpliconPhasing_Error(f'unknown method: {args.method}')
        splitter = sclass(args.inFile,
                           args.reference,
                           region=args.region,
                           truncate=args.truncate,
                           minFrac=args.minFrac,
                           minReads=args.minReads,
                           minSignal=args.minSignal,
                           aggressive=args.aggressive,
                           indels=args.indels,
                           vTable=varDf,
                           prefix=prefix,
                           nproc=args.nproc, 
                           makeDf=makeDf,
                           log=log,stats=stats)
    elif args.method == 'debruijn':
        from src.phase.split import sparseDBG
        splitter = sparseDBG(args.k,
                             collapse=args.maxHP,
                             ignoreEnds=args.ignore,
                             minReads=args.minReads,
                             minFrac=args.minFrac,log=log)
        splitter.loadReads(args.inFile,
                           region=args.region,
                           minLength=args.minLength,
                           maxLength=args.maxLength)
    else:
        raise LongAmpliconPhasing_Error(f'unknown method: {args.method}')

    phaser = Phaser(splitter,
                    args.sampleName,
                    aggressive=args.aggressive,
                    log=log)
    phaser.run()
    
    #export stuff
    #bam
    if not args.noBam and ftype == 'bam':
        log.info('Writing clustered BAM(s)')
        #hacky
        args.inBAM = args.inFile
        outBAM = f'{prefix}hptagged.bam' 
        addHPtag(args,outBAM,phaser.clusterMap)
        if args.region:
            log.info('Writing region BAM')
            outBED = f'{prefix}region.bed'
            writeSimpleBED(*getCoordinates(args.region),
                            args.sampleName,
                            phaser.splitter.stats['clustered reads'],
                            outBED)
            outBAM = f'{prefix}region.bam'
            writeRegionBam(args.inFile,outBAM,args.region)
    #fastq
    if args.exportFq:
        if args.inFile.endswith('a'):
            log.warning('Cannot export fasta from fasta')
        else:
            ft = 'bam' if args.inFile.endswith('bam') else 'fastq'   
            log.info("Exporting Fastq files")
            exportFastq(args.inFile,ft,prefix,
                        phaser.clusterMap,region=args.region)
    #simple decision tree
    log.info('Writing Splits')
    with open(f'{prefix}clusterSplits.txt','w') as ofile:
        ofile.write(f'{str(phaser)}\n')
    #cluster file/readnames
    log.info('Writing Cluster file')
    with open(f'{prefix}clusters.txt','w') as ofile:
        for node,clust in phaser.node2cluster.items():
            vclust = phaser.vTree[node]
            ofile.write(f'>{vclust.clusterName()}\n')
            ofile.write('\n'.join(vclust.reads) + '\n')
    #check if anything was produced, else exit
    hasResults = True
    if len(phaser.node2cluster) == 0:
        nreads = phaser.splitter.stats['clustered reads']
        log.warning(f'No results produced for sample {phaser.sampleName}.  Nreads: {nreads}')
        hasResults = False
    #variants/consensus
    if args.variants:
        from src.phase.caller import VariantCaller,FIGFORMAT 
        if args.method == 'debruijn' or args.maxHP != 0:
            #need to build the variant splitter
            log.info('Building uncompressed variant table from bam')
            from src.phase.split import VariantGrouper
            vsplitter = VariantGrouper(args.inFile,
                                       args.reference,
                                       region=args.region,
                                       truncate=args.truncate,
                                       minFrac=args.minFrac,
                                       minReads=args.minReads,
                                       minSignal=args.minSignal,
                                       aggressive=args.aggressive,
                                       indels=args.indels,
                                       vTable=None,log=log,
                                       prefix=prefix,
                                       makeDf=makeDf,
                                       nproc=args.nproc)
        else:
            vsplitter = splitter

        endpoints = vsplitter.vTable.groupby(phaser.clusterMap).apply(vsplitter.getEndpoints)
        caller    = VariantCaller(vsplitter.sigVar,
                                  vsplitter.minCount,
                                  args.reference,
                                  phaser.clusterMap,
                                  endpoints,
                                  log=log)
        #reads counts
        log.info('Writing Run Summary')
        name = f'{prefix}summary.txt'
        with open(name,'w') as ofile:
            ofile.write(f'{summary(splitter,caller)}\n')

        #summary table
        log.info('Writing Cluster Summary')
        name = f'{prefix}alleleClusterSummary.csv'
        caller.alleleFrequency(caller.clusterMap).to_csv(name)
        log.info('Writing Sample Variant Summary')
        name = f'{prefix}sampleVariantSummary.csv'
        caller.alleleFrequency(caller.variantGroupMap)\
              .reindex(list(set(caller.variantGroupMap.values())))\
              .to_csv(name)
        
        #total allele fractions (draft)
        log.info('Writing Variant Fraction File')
        name = f'{prefix}variantFraction.csv'
        caller.variantFractions.to_csv(name)
        #draft consensus
        log.info('Generating draft consensus')
        name = f'{prefix}draft.consensus.fasta'
        with open(name,'w') as ofasta:
            for clust,seq in caller.draftConsensus().items():
                clust = clust if clust!=-1 else 'Noise'
                ofasta.write(f'>cluster{clust}.plurality\n')
                ofasta.write(seq + '\n')
        #entropy
        log.info('Writing entropy table')
        name = f'{prefix}entropy.csv'
        caller.entropy.to_csv(name)
        if args.entropyPlot:
            log.info('Writing entropy heatmap')
            name = f'{prefix}entropy.{FIGFORMAT}'
            caller.entropyPlot(name)
    
    log.info('Done')

    for handler in log.handlers:
        handler.close()
        log.removeFilter(handler)

    return phaser

###multiproc definitions for pickling

def makeDf(bamfile,reference,region=None,truncate=False):
    bam = pysam.AlignmentFile(bamfile,'r')
    ref = pysam.FastaFile(reference)
    df  = pd.DataFrame({(column.reference_name,
                         column.reference_pos)  : dict(zip(column.get_query_names(),
                                                        column.get_query_sequences(mark_matches=True,
                                                                                   add_indels=True)))
                       for column in bam.pileup(flag_filter=0x900,
                                                 fastafile=ref,region=region,
                                                 truncate=truncate,
                                                 min_base_quality=0,
                                                 compute_baq=False)})
    df.columns.names = ('contig','pos')
    return df


########

class LongAmpliconPhasing_Error(Exception):
    pass

if __name__ == '__main__':
    import argparse,sys

    parser = argparse.ArgumentParser(prog='LongAmpliconPhasing.py', description='Recursively separate amplicons by shared variants')
    parser.set_defaults(prog=parser.prog)
    parser.add_argument('inFile', metavar='inFile', type=str,
                    help='Fasta,Fastq, or BAM containing CCS reads')
    parser.add_argument('sampleName', metavar='sampleName', type=str,
                    help='sample name')
    parser.add_argument('-j', dest='nproc', type=int, default=1,
                    help=f'Number of procs to use for loading data')
    parser.add_argument('--reference', dest='reference', type=str, default=None,
                    help='reference fasta used for alignment.  Required for draft consensus and variant counts')
    parser.add_argument('--region', dest='region', type=str,default=None,
                    help='reference subregion. Format [chr]:[start]-[stop].  Default None.')
    parser.add_argument('-t','--truncate', dest='truncate', action='store_true', default=False,
                    help='Truncate to columns only within region (pysam pileup).  Default False')
    parser.add_argument('-k','--kmer', dest='k', type=int, default=DEFAULTKMER,
                    help=f'Kmer size for debruijn method.  Ignored for align method. Default {DEFAULTKMER}')
    parser.add_argument('--no-indels', dest='indels', action='store_false', default=True,
                    help='Do not split on indels (only for alignm method).  Default use indels')
    parser.add_argument('-i','--ignore', dest='ignore', type=int, default=0,
                    help=f'Ignore first and last N bases when building debruijn graph. Not used for align method. Default 0')
    parser.add_argument('-p','--prefix', dest='prefix', type=str, default=None,
                    help=f'Output prefix. Default cwd')
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
    parser.add_argument('-F','--exportFastq', dest='exportFq', action='store_true', default=False,
                    help='Export Fastq per phase. Default False')
    parser.add_argument('--template', dest='template', choices=['first','median'], default='median',
                    help='Method for choosing reference template from inputs (if no reference passed). Default median')
    parser.add_argument('-m','--method', dest='method', choices=['align','debruijn','cluster'], default='align',
                    help='Splitting method.  If align and maxHP != 0, reads will be realigned after compression for clustering; Output variants are from non-compressed pileup. Default align')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=False,
                    help='Print progress messages to stdout. Default False')
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
