from src.model.kmer   import *
from src.model.models import MODELS
from src.utils.bam import addHPtag
from src.utils.extract import Extract_Exception

CLUSTNAME = 'cluster{0}_numreads{1}'.format #args: int,int

def main(args):
    if args.normalize == 'none':
        args.normalize = None

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
                     agg        =args.agg,
                     extractRef =args.reference)

    if args.testPlot:
        #plot kdist and exit
        from src.figures.kdist import plotEPS
        print("Plotting distance to m-neighbors")
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
        cnts    = sorted([len(idx) for c,idx in grouped.groups.items() if c!=-1],reverse=True)
        print(f'Writing clusters with nreads {",".join(map(str,cnts))}')
        if -1 in grouped.groups:
            print(f'{len(grouped.groups[-1])} reads identified as noise')
        for cluster,reads in grouped:
            nreads = len(reads)
            name = f'Noise_numreads{nreads}' if cluster==-1 else CLUSTNAME(cluster,nreads)
            namefile.write(f'>{name}\n')
            namefile.write('\n'.join(reads.index) + '\n')

    if not args.noBam:
        print("Adding HP tag to bam")
        clusterMap = dict(zip(data.index,clusterIdx))
        outBam     = f'{args.prefix}.hptagged.bam'
        addHPtag(args.inBAM,outBam,clusterMap,dropNoClust=args.drop,splitBam=args.splitBam)

    if args.plotReads:
        from src.figures.cluster import plotReads
        fig = plotReads(data,clusterIdx)
        fig.savefig(f'{args.prefix}.clusters.png')

    return None

def printParams(model):
    return '\n'.join(['\t' + '='.join(map(str,v)) for v in model.defaults.items()])
