import re,pysam
from collections import Counter
from multiprocessing import Pool,cpu_count
import pandas as pd
import numpy as np
import mappy as mp
import skbio.diversity.alpha as ad
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.cluster import FeatureAgglomeration
from .models import Clustering_Exception
from ..utils.extract import getCoordinates, \
                            extractRegion

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

def simpsonFilter(data,maxDom):
    dom = data.apply(lambda feat: ad.dominance(np.bincount(feat)))
    return data.loc[:,dom<=maxDom]

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
            recGen = extractRegion(inBAM,extractRef,region)
        else:
            recGen = bam.fetch(*getCoordinates(region))
    else:
        recGen = bam

    counter = kmerCounter(k,collapseHP=collapse,
                          minimizer=minimizer,
                          ignoreEnds=ignoreEnds)
    if whitelist:
        wl      = open(whitelist).read().split()
        useRead = (lambda read: read in wl)
    else:
        useRead = noFilter
    passQuality = qualityCrit(qual)
    flankCrit   = getFlankCrit(flanks) if flanks else noFilter

    print("Reading Sequence")
    proc = nproc if nproc else cpu_count()
    pool = Pool(proc)
    records = pool.map(recordMaker(counter),
                       [(rec.query_name,rec.query_sequence) for rec in recGen 
                        if useRead(rec.query_name) and passQuality(rec) and flankCrit(rec)])
    #records = [recordMaker(counter)((rec.query_name,rec.query_sequence))
    #           for rec in recGen 
    #           if useRead(rec.query_name) and passQuality(rec) and flankCrit(rec)]
    data = pd.concat(records,axis=1,sort=False)\
             .fillna(0).astype(int).T
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

class kmerCounter:
    def __init__(self,k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
        self.k         = k
        self.transform = hpCollapse if collapseHP else ident
        self.minim     = getMinimizer(minimizer) if minimizer>0 else ident
        self.start     = ignoreEnds
        self.end       = -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
    def __call__(self, seq):
        s = self.transform(seq[self.start:self.end])
        return Counter(self.minim(s[i:i+self.k]) for i in range(len(s)-self.k))
    

class recordMaker:
    def __init__(self, counter):
        self.counter = counter
    def __call__(self, rec):
        return pd.Series(self.counter(rec[1]),
                         name=rec[0])

class Kmer_Exception(Exception):
    pass

