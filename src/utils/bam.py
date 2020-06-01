import pysam,re
import seaborn as sns
import numpy as np
from .extract import getCoordinates, \
                     fastqReader

PALETTE     ='husl'
NOCLUST     =999
NOCLUSTCOLOR='255,255,255'
NOISECOLOR  ='200,200,200'
DUMMYVERSION=0.1

def addHPtag(args,outBAM,clusterMap,noCluster=NOCLUST):
    '''clusterMap is map of {readname:cluster::int}'''
    cvals   = set(clusterMap.values())
    ncolors = len(cvals)
    rgb     = np.array(sns.color_palette(PALETTE,ncolors))
    colors  = {clust:','.join(map(str,(col*255 +1).astype(int)))
               for clust,col in zip(cvals,rgb)}
    colors[noCluster] = NOCLUSTCOLOR
    #noise
    colors[-1] = NOISECOLOR
    dropFilt   = args.drop

    with pysam.AlignmentFile(args.inBAM) as inbam:
        header = inbam.header.to_dict()
        header['PG'].append(getCL(args))
        if args.splitBam:
            names  = [re.sub('.bam$',
                             f'.{"Noise" if c==-1 else str(c)}.bam',
                             outBAM)
                      for c in cvals ]
            outMap = {c : pysam.AlignmentFile(n,'wb',header=header)
                      for c,n in zip(cvals,names)}
            dropFilt = True #no place to put these
        else:
            names  = [outBAM]
            outbam = pysam.AlignmentFile(outBAM,'wb',header=header)
            outMap = {c:outbam for c in cvals.union([noCluster])}
        getbam = (lambda c: outMap[c])
        recGen = inbam.fetch(*getCoordinates(args.region)) if args.region else inbam
        for rec in recGen:
            if dropFilt:
                if rec.flag & 0x900:
                    continue
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

def getCL(args):
    baseProg = args.prog.split(' ')[0]
    cl = f'{args.prog} ' + ' '.join(f'--{a} {getattr(args,a)}' 
                                    for a in vars(args)
                                    if getattr(args,a) and not a=='func') 
    return {'ID': baseProg,
            'PN': baseProg,
            'VN': DUMMYVERSION,
            'CL': cl}

def fqRec(name,seq,qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'

def stripReadname(readname):
    '''strip off all but <movie>/<zmw>/ccs'''
    return '/'.join(readname.split('/')[:3])

def exportFastq(inFile,fileType,outPrefix,clusterMap,region=None):
    cvals  = set(clusterMap.values())
    ofiles ={c : open(f'{outPrefix}.cluster{c}.fastq','w') 
             for c in cvals if c!=-1}
    if fileType == 'bam':
        inBam  = pysam.AlignmentFile(inFile)
        recGen = inBam.fetch(*getCoordinates(region)) if region else inBam 
    elif fileType == 'fastq':
        recGen = fastqReader(inFile)
    else:
        raise Bam_Exception(f'Invalid input type {fileType}')

    for rec in recGen:
        if rec.flag & 0x900:
            continue
        stripped = stripReadname(rec.query_name)
        if stripped in clusterMap:
            cluster = clusterMap[stripped]
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

class Bam_Exception(Exception):
    pass
