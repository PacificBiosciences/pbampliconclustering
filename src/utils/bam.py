import pysam
import seaborn as sns
import numpy as np

PALETTE='husl'
NOCLUST=999

def addHPtag(inBAM,outBAM,clusterMap,splitBam=False,noCluster=NOCLUST,dropNoClust=False):
    '''clusterMap is map of {readname:cluster::int}'''
    cvals   = set(clusterMap.values())
    ncolors = len(cvals)
    rgb     = np.array(sns.color_palette(PALETTE,ncolors))
    colors  = {clust:','.join(map(str,(col*255 +1).astype(int)))
               for clust,col in zip(cvals,rgb)}
    colors[noCluster] = '255,255,255'
    #noise
    colors[-1] = '200,200,200'

    with pysam.AlignmentFile(inBAM) as inbam:
        if splitBam:
            names  = [re.sub('.bam$',
                             '.%s.bam'%('Noise' if c==-1 else str(c)),
                             outBAM)
                      for c in cvals ]
            outMap = {c : pysam.AlignmentFile(n,'wb',template=inbam)
                      for c,n in zip(cvals,names)}
            dropNoClust = True #no place to put these
        else:
            names  = [outBAM]
            outbam = pysam.AlignmentFile(outBAM,'wb',template=inbam)
            outMap = {c:outbam for c in cvals.union([noCluster])}
        getbam = (lambda c: outMap[c])
        for rec in inbam:
            if dropNoClust:
                if rec.query_name not in clusterMap: #filtered
                    continue
                elif clusterMap[rec.query_name] == -1: #noise
                    continue
            clust = int(clusterMap.get(rec.query_name,noCluster))
            rec.set_tag('YC',colors[clust])
            rec.set_tag('HP',clust)
            getbam(clust).write(rec)
    for b,n in zip(outMap.values(),names):
        b.close()
        pysam.index(n)

    return None
