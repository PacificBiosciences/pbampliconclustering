import re,pysam
from collections import defaultdict
from multiprocessing import Pool,Manager,cpu_count
import pandas as pd
import numpy as np
import mappy as mp
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.cluster import FeatureAgglomeration
from ..utils.extract import getCoordinates, \
                            extractRegion

FLANKSIZE=500
MINLEN   =500
MAXLEN   =25000
MAXPALOVR=250 #max rev comp overlap bwtn primary/supp alignments
RANDSEED =17

_PATT = re.compile(r'([ATGC])\1+')
def hpCollapse(seq):
    return _PATT.sub(r'\1',seq)

def qualityCrit(minqv):
    def crit(rec):
        return rec.get_tag('rq') >= minqv #and not (rec.flag & 0x900)
    return crit

def getLengthCrit(minLen,maxLen):
    def crit(rec):
        rlen = rec.query_length
        return rlen >= minLen and rlen <= maxLen
    return crit

def getFlankCrit(flankFa):
    flanks = [seq for _,seq,_ in mp.fastx_read(flankFa)]
    def includesFlanks(rec):
        a = mp.Aligner(seq=rec.query_sequence,preset='sr')
        return np.all([len(list(a.map(f)))==1 for f in flanks])
    return includesFlanks

def getMinimizer(m=6):
    def minimizer(seq):
        return sorted(seq[i:i+m] for i in range(0,len(seq)-m+1))[0]
    return minimizer

def noFilter(x):
    return True

def ident(x):
    return x

def overlap(a,b):
    return min(a.r_en,b.r_en) - max(a.r_st,b.r_st)

def isArtifact(alns,maxOvr=MAXPALOVR):
    for i,aln1 in alns.iterrows():
        for j,aln2 in alns.loc[i:].iterrows():
            if i==j: continue
            if overlap(aln1,aln2) > maxOvr and aln1.isrev!=aln2.isrev:
                return True
    return False

def loadKmers(inBAM,qual,k,nproc=0,
              collapse=True,region=None,
              minLength=MINLEN,maxLength=MAXLEN,
              minimizer=0,ignoreEnds=0,
              whitelist=None,flanks=None,
              trim=0,norm=None,
              components=3,agg='pca',
              extractRef=None,palfilter=True,
              exportKmers=None,subsample=0,
              randseed=RANDSEED):
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
    lengthCrit  = getLengthCrit(minLength,maxLength)
    flankCrit   = getFlankCrit(flanks) if flanks else noFilter

    print("Reading Sequence")
    seqDB = pd.DataFrame([{'qname'   :rec.query_name,
                           'seq'     :rec.query_sequence,
                           'r_st'    :rec.reference_start,
                           'r_en'    :rec.reference_end,
                           'isSecond':bool(rec.flag & 0x900),
                           'isrev'   :bool(rec.flag & 0x10)}
                          for rec in recGen
                          if useRead(rec.query_name) 
                            and passQuality(rec)
                            and lengthCrit(rec) 
                            and flankCrit(rec)])    

    if palfilter:
        #filter out palindromic sequences
        multiMapped = seqDB[seqDB.isSecond].qname
        palindromes = seqDB.query('qname in @multiMapped')\
                           .groupby('qname')\
                           .filter(isArtifact)\
                           .qname.unique()
        sequences = seqDB.query('qname not in @palindromes and not isSecond')
    else:
        #just remove secondary/supplemental alignments
        sequences = seqDB.query('not isSecond')

    if len(sequences) == 0:
        raise Kmer_Exception('No sequences returned for clustering!')

    if subsample and len(sequences) > subsample:
        print(f"Downsampling to {subsample} reads from {len(sequences)}")
        sequences = sequences.sample(subsample,replace=False,random_state=randseed)

    counts    = defaultdict(lambda: np.zeros(len(sequences),dtype=np.int16))
    parser    = seqParser(k,collapseHP=collapse,
                          minimizer=minimizer,
                          ignoreEnds=ignoreEnds)
    for i,seq in enumerate(sequences.seq):
        for kmer in parser(seq):
            counts[kmer][i] +=1

    data = pd.DataFrame(np.array(list(counts.values())).T,
                        index=sequences.qname)

    if trim:
        #print("Trimming low-freq kmers")
        print("Trimming high- and low-freq kmers")
        freqs = data.sum()/len(data)
        #data  = data.loc[:,freqs>=trim] 
        data  = data.loc[:,(freqs>=trim) & (freqs<=(1-trim))] 

    if exportKmers:
        print('Exporting kmer counts')
        data.rename(columns=dict(enumerate(counts.keys()))).to_csv(exportKmers)

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
        for i in range(len(s)-self.k+1):
            yield self.minim(s[i:i+self.k])

class Kmer_Exception(Exception):
    pass

