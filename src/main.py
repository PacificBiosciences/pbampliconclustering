from src.model.kmer   import *
from src.model.models import MODELS
from src.utils.bam import addHPtag,exportFastq
from src.utils.clust import clusterName
from src.utils.extract import Extract_Exception

def main(args):
    if args.normalize == 'none':
        args.normalize = None

    #load dataframe with samples(row) by kmer counts (cols)
    data = loadKmers(args.inBAM,args.minQV,args.kmer,
                     nproc      =args.njobs,
                     collapse   =args.hpCollapse,
                     region     =args.region,
                     minimizer  =args.minimizer,
                     ignoreEnds =args.ignoreEnds,
                     whitelist  =args.whitelist,
                     flanks     =args.flanks,
                     #simpson    =args.simpson,
                     trim       =args.trim,
                     norm       =args.normalize,
                     components =args.components,
                     agg        =args.agg,
                     extractRef =args.reference)

    #Plot k-nearest neighbors
    if args.testPlot:
        from src.figures.kdist import plotEPS
        print(f"Plotting distance to {args.minReads} neighbors")
        f = plotEPS(data,args.minReads,args.normalize)
        f.savefig(f'{args.prefix}.eps_estimator.png')
        return data

    #Clustering
    cluster     = MODELS[args.model](args)
    print(f'Clustering {len(data)} reads with {args.model}\n{printParams(cluster)}')
    result      = cluster.fit(data)
    clusterIdx  = result.labels_

    #TODO
    #cluster size and warning if too much noise as frac of total

    #write cluster file
    with open(f'{args.prefix}.clusters.txt', 'w') as namefile:
        grouped = data.groupby(clusterIdx)
        cnts    = sorted([len(idx) for c,idx in grouped if c!=-1],reverse=True)
        print(f'Writing {len(cnts)} clusters with nreads {",".join(map(str,cnts))}')
        if -1 in grouped.groups:
            print(f'{len(grouped.groups[-1])} reads identified as noise')
        for clust,reads in grouped:
            nreads = len(reads)
            name = f'Noise_numreads{nreads}' if clust==-1 else clusterName((clust,nreads))
            namefile.write(f'>{name}\n')
            namefile.write('\n'.join(reads.index) + '\n')

    #tag BAM
    if not args.noBam:
        print("Adding HP tag to bam")
        names      = data.index.map(stripReadname)
        clusterMap = dict(zip(names,clusterIdx))
        outBam     = f'{args.prefix}.hptagged.bam'
        addHPtag(args.inBAM,outBam,clusterMap,region=args.region,dropNoClust=args.drop,splitBam=args.splitBam)

    #export fastq
    if args.fastq:
        print("Exporting fastq")
        exportFastq(args.inBAM,args.prefix,clusterMap,region=args.region)

    #plot samples
    if args.plotReads:
        from src.figures.cluster import plotReads
        fig = plotReads(data,clusterIdx)
        fig.savefig(f'{args.prefix}.clusters.png')

    return data,cluster,result

def stripReadname(readname):
    '''strip off all but <movie>/<zmw>/ccs'''
    return '/'.join(readname.split('/')[:3])

def printParams(model):
    return '\n'.join(['\t' + '='.join(map(str,v)) for v in model.defaults.items()])
