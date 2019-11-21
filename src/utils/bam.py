import pysam,re
import seaborn as sns
import numpy as np
from .extract import getCoordinates

PALETTE     ='husl'
NOCLUST     =999
NOCLUSTCOLOR='255,255,255'
NOISECOLOR  ='200,200,200'

def addHPtag(inBAM,outBAM,clusterMap,region=None,splitBam=False,noCluster=NOCLUST,dropNoClust=False):
    '''clusterMap is map of {readname:cluster::int}'''
    cvals   = set(clusterMap.values())
    ncolors = len(cvals)
    rgb     = np.array(sns.color_palette(PALETTE,ncolors))
    colors  = {clust:','.join(map(str,(col*255 +1).astype(int)))
               for clust,col in zip(cvals,rgb)}
    colors[noCluster] = NOCLUSTCOLOR
    #noise
    colors[-1] = NOISECOLOR

    with pysam.AlignmentFile(inBAM) as inbam:
        if splitBam:
            names  = [re.sub('.bam$',
                             f'.{"Noise" if c==-1 else str(c)}.bam',
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
        recGen = inbam.fetch(*getCoordinates(region)) if region else inbam
        for rec in recGen:
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

def fqRec(name,seq,qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'

def exportFastq(inBAM,outPrefix,clusterMap,region=None):
    cvals  = set(clusterMap.values())
    ofiles ={c : open(f'{outPrefix}.cluster{c}.fastq','w') 
             for c in cvals if c!=-1} 
    with pysam.AlignmentFile(inBAM) as inbam:
        recGen = inbam.fetch(*getCoordinates(region)) if region else inbam
        for rec in recGen:
            if rec.flag & 0x900:
                continue
            if rec.query_name in clusterMap:
                cluster = clusterMap[rec.query_name]
                if cluster == -1: #noise
                    continue
                #TODO?? export in native orientation (revcomp if neg strand)
                qual = ''.join([chr(q+33) for q in rec.query_qualities])
                ofiles[cluster].write(fqRec(rec.query_name,
                                            rec.query_sequence,
                                            qual))
    for c,f in ofiles.items():
        f.close()

    return None

