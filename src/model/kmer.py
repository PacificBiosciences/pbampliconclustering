import re,pysam
from collections import Counter
import pandas as pd
import numpy as np
import mappy as mp
import skbio.diversity.alpha as ad
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.cluster import FeatureAgglomeration
from .models import Clustering_Exception

_PATT = re.compile(r'([ATGC])\1+')
def hpCollapse(seq):
    return _PATT.sub(r'\1',seq)

def countKmers(k=11,collapseHP=True,minimizer=0,ignoreEnds=0):
    transform =  hpCollapse if collapseHP else (lambda x:x)
    minim     =  getMinimizer(minimizer) if minimizer>0 else (lambda x:x)
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

def getCoordinates(regionString):
    try:
        ctg,start,stop = re.search('(.*):(\d+)-(\d+)',regionString).groups()
    except AttributeError as e:
            #catch when the region format doesn't match
            raise Clustering_Exception('Invalid region format %s. Correct \'[chr]:[start]-[stop]\'' % regionString)
    return ctg.strip(),int(start),int(stop)

def simpsonFilter(data,maxDom):
    print('Filtering by dominance')
    useMers = data.apply(pd.Series.value_counts)\
                  .fillna(0).astype(int)\
                  .apply(ad.dominance) < maxDom
    return data.loc[:,data.columns[useMers]]

def noFilter(x):
    return True

def loadKmers(inBAM,qual,k,
              collapse=True,region=None,
              minimizer=0,ignoreEnds=0,
              whitelist=None,flanks=None,
              simpson=None,norm=None,
              components=3,agg='pca',
              extractRef=None):
    '''
    kmer loader
    '''
    bam         = pysam.AlignmentFile(inBAM)
    if region:
        if extractRef:
            from ..utils.extract import extractRegion
            recGen = extractRegion(inBAM,extractRef,region)
        else:
            recGen = bam.fetch(*getCoordinates(region))
    else:
        recGen = bam

    kmerCounter = countKmers(k,collapseHP=collapse,
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
    data = pd.concat([pd.Series(kmerCounter(rec.query_sequence),
                                name=rec.query_name)
                      for rec in recGen
                      if useRead(rec.query_name)
                         and passQuality(rec)
                         and flankCrit(rec)],
                     axis=1,sort=True)\
             .fillna(0).astype(int).T
    if norm:
        print('Normalizing data')
        data = pd.DataFrame(normalize(data,norm=norm),
                            index=data.index,
                            columns=data.columns)
    if components:
        print('Reducing Features with %s' % agg)
        if agg == 'pca':
            tool = PCA(n_components=components)
        elif agg == 'featagg':
            tool = FeatureAgglomeration(n_clusters=components)
        else:
            raise Kmer_Error('%s is not a valid reduction tool' % agg)
        data = pd.DataFrame(tool.fit_transform(data),
                            index=data.index)

    return data

class Kmer_Error(Exception):
    pass

