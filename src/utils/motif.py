import re
import numpy as np
import pandas as pd
from collections import Counter
from src.utils.clust import clusterName

DEFAULTCI= [0.025,0.975] #95%

def countMotifs(motifs,lengthField=False):
    ''' takes a list of motifs
        returns a function that takes a sequence and returns a dict of counts
        optionally includes length of seq in output
    '''
    patt = re.compile('(' + ')|('.join(motifs) + ')')
    def getCounts(seq):
        counts = Counter(m.group() for m in patt.finditer(seq))
        if lengthField:
            counts[lengthField] = len(seq)
        return counts
    return getCounts

def getCounts(seqGen,motifs,lengthField=None):
    motifCounter = countMotifs(motifs,lengthField)
    data         = {name:motifCounter(seq)
                    for name,seq in seqGen}
    columns      = motifs+[lengthField] if lengthField else motifs
    return pd.DataFrame(data.values(),index=data.keys())\
             .reindex(columns,axis=1)\
             .fillna(0)

def resampleCI(data,nboot=10000,ci=DEFAULTCI):
    n = max(nboot,len(data))
    resamp = np.random.choice(data,size=n,replace=True)
    return '({} - {})'.format(*map(int,np.quantile(resamp,ci)))

def clusterStats(motifCounts,clusterIdx,outColumns,
                 aggFuncs=[np.median,np.mean],
                 randomSeed=None,ci=DEFAULTCI):
    '''
    motifCounts: df with cols =  motif (+length), index = readnames
    clusterIdx: vector of cluster indices, same order as motifCounts.index
    outColumns: list of column names to describe in output
    aggFuncs: list of functions to apply to each column
    randomSeed: random seed for resampling ci
    '''
    clusters    = motifCounts.groupby(clusterIdx)
    clusterSize = clusters.size().rename(('Read','count'))

    #set random seed
    np.random.seed(randomSeed)
    results = clusters[outColumns].agg(aggFuncs+[resampleCI])\
                                  .join(clusterSize)
    #rename ci column
    name = f'ci{int(100*(ci[1]-ci[0]))}'
    results.rename(columns={'resampleCI':name},level=1,inplace=True)

    #rename clusters
    names = clusterSize.reset_index().apply(clusterName,axis=1)
    results.index = names.values

    return results
