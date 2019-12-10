import re,pysam
from collections import defaultdict
from multiprocessing import Pool,Manager,cpu_count
import pandas as pd
import numpy as np
import mappy as mp
#import skbio.diversity.alpha as ad
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.cluster import FeatureAgglomeration
from ..utils.extract import getCoordinates, \
                            extractRegion

FLANKSIZE=500

_PATT = re.compile(r'([ATGC])\1+')
def hpCollapse(seq):
    return _PATT.sub(r'\1',seq)

def countKmers(k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
    transform =  hpCollapse if collapseHP else ident
    minim     =  getMinimizer(minimizer) if minimizer>0 else ident
    end       =  -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
    def getKmers(seq):
        s = transform(seq[ignoreEnds:end])
        return Counter(minim(s[i:i+k]) for i in range(len(s)-k))
    return getKmers

def qualityCrit(minqv):
    def crit(rec):
        return rec.get_tag('rq') >= minqv and not (rec.flag & 0x900)
    return crit

def getFlankCrit(flankFa):
    flanks = [seq for _,seq,_ in mp.fastx_read(flankFa)]
    def includesFlanks(rec):
        a = mp.Aligner(seq=rec.query_sequence,preset='sr')
        return np.all([len(list(a.map(f)))==1 for f in flanks])
    return includesFlanks

def getMinimizer(m=6):
    def minimizer(seq):
        return sorted(seq[i:i+m] for i in range(0,len(seq)-m))[0]
    return minimizer

#def simpsonFilter(data,maxDom):
#    dom = data.apply(lambda feat: ad.dominance(np.bincount(feat)))
#    return data.loc[:,dom<=maxDom]

def noFilter(x):
    return True

def ident(x):
    return x

def loadKmers(inBAM,qual,k,nproc=0,
              collapse=True,region=None,
              minimizer=0,ignoreEnds=0,
              whitelist=None,flanks=None,
              trim=0,norm=None,
              components=3,agg='pca',
              extractRef=None):
    '''
    kmer loader
    '''
    bam         = pysam.AlignmentFile(inBAM,'rb')
    if region:
        if extractRef:
            recGen = extractRegion(inBAM,extractRef,region,flanksize=FLANKSIZE)
        else:
            recGen = bam.fetch(*getCoordinates(region))
    else:
        recGen = bam

    if whitelist:
        wl      = open(whitelist).read().split()
        useRead = (lambda read: read in wl)
    else:
        useRead = noFilter
    passQuality = qualityCrit(qual)
    flankCrit   = getFlankCrit(flanks) if flanks else noFilter

    print("Reading Sequence")
    sequences =  pd.DataFrame([{'qname':rec.query_name,'seq':rec.query_sequence}
                               for rec in recGen
                               if useRead(rec.query_name) and passQuality(rec) and flankCrit(rec)])    
    if len(sequences) == 0:
        raise Kmer_Exception('No sequences returned for clustering!')

    counts    = defaultdict(lambda: np.zeros(len(sequences),dtype=np.int16))
    parser    = seqParser(k,collapseHP=collapse,
                          minimizer=minimizer,
                          ignoreEnds=ignoreEnds)
    for i,seq in enumerate(sequences.seq):
        for kmer in parser(seq):
            counts[kmer][i] +=1

    data = pd.DataFrame(np.array(list(counts.values())).T,
                        index=sequences.qname)
#    if simpson>0:
#        print('Filtering by dominance')
#        data = simpsonFilter(data,simpson)

    if trim:
        print("Trimming low-freq kmers")
        freqs = data.sum()/len(data)
        data  = data.loc[:,freqs>=trim] 

    if norm:
        print('Normalizing data')
        data = pd.DataFrame(normalize(data,norm=norm),
                            index=data.index,
                            columns=data.columns)
    if components:
        print(f'Reducing Features with {agg}')
        if agg == 'pca':
            tool = PCA(n_components=components)
        elif agg == 'featagg':
            tool = FeatureAgglomeration(n_clusters=components)
        else:
            raise Kmer_Error(f'{agg} is not a valid reduction tool')
        data = pd.DataFrame(tool.fit_transform(data),
                            index=data.index)

    return data

class seqParser:
    def __init__(self,k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
        self.k         = k
        self.transform = hpCollapse if collapseHP else ident
        self.minim     = getMinimizer(minimizer) if minimizer>0 else ident
        self.start     = ignoreEnds
        self.end       = -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
    def __call__(self, seq):
        s = self.transform(seq[self.start:self.end])
        for i in range(len(s)-self.k):
            yield self.minim(s[i:i+self.k])

#class kmerCounter:
#    def __init__(self,k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
#        self.k         = k
#        self.transform = hpCollapse if collapseHP else ident
#        self.minim     = getMinimizer(minimizer) if minimizer>0 else ident
#        self.start     = ignoreEnds
#        self.end       = -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
#    def __call__(self, seq):
#        s = self.transform(seq[self.start:self.end])
#        return Counter(self.minim(s[i:i+self.k]) for i in range(len(s)-self.k))

class kmerCounter2:
    def __init__(self,size,k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
        self.size      = size
        self.k         = k
        self.transform = hpCollapse if collapseHP else ident
        self.minim     = getMinimizer(minimizer) if minimizer>0 else ident
        self.start     = ignoreEnds
        self.end       = -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
    def __call__(self, i, seq, multiDict):
        s = self.transform(seq[self.start:self.end])
        for j in range(len(s)-self.k):
            kmer = self.minim(s[j:j+self.k])
            if kmer not in multiDict:
                multiDict[kmer] = np.zeros(self.size,int)
            arr = multiDict[kmer]
            arr[i] +=1
            multiDict[kmer] = arr
            print(f'adding {i}:{arr[:i+1]}  :  {kmer} {multiDict[kmer][:i+1]}')

class recordMaker:
    def __init__(self, counter):
        self.counter = counter
    def __call__(self, rec):
        return pd.Series(self.counter(rec[1]),
                         name=rec[0])

class Kmer_Exception(Exception):
    pass

