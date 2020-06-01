import pysam,os,re
from math import ceil
import pandas as pd
import numpy as np
from itertools import chain
from collections import Counter
from scipy.stats import entropy
from operator import itemgetter
from collections import defaultdict
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils.validation import check_symmetric
from .utils import hpCollapse,RecordGenerator
from ..utils.extract import getCoordinates

DIAGNOSTICS=False

class VariantGrouper:
    def __init__(self,inFile,refFasta,
                 region=None,truncate=False,
                 minCov=1,minFrac=0.1,
                 minReads=10,minSpan=0.9,
                 minSignal=0.05,flagFilter=0x900,
                 aggressive=False,indels=True,
                 hpmask=0,hptol=0,
                 vTable=None,nproc=1,prefix=None,
                 makeDf=None,log=None,stats={},
                 diagnostics=DIAGNOSTICS):
        self.bamfile    = inFile             #bam file
        self.refFasta   = refFasta           #ref
        self.region     = region             #region eg X:12345-678900
        self.truncate   = truncate           #truncate pileup to region (pysam pileup kwarg)
        self.minCov     = minCov             #min pos coverage for inclusion
        self._minFrac   = minFrac            #min fraction to split off
        self._minReads  = minReads           #abs min reads to split off
        self.minSpan    = minSpan            #not implem yet
        self.minSignal  = minSignal          #min frac of reads that are diff from ref to consider
        self.flagFilter = flagFilter         #hex filters, passed to pysam.pileup engine
        self.aggressive = aggressive
        self.indels     = indels             #use indels for variants
        self.hpmask     = hpmask             #mask variants in homopolymers larger than this size; incl ins at beginning
        self.hptol      = hptol              #mask tolerance (added to hpmask on either side of each hp)
        self.nproc      = nproc
        self.prefix     = prefix
        self.makeDf     = makeDf             #pickle-able function for parallel processing
        self.log        = log
        self.stats      = stats
        self.vTable     = vTable if vTable is not None else self._makeDf()
        self.sigVar     = self._makeSigVar()
        self.minCount   = self._getMinCount()
        self.readnames  = self.sigVar.index
        if diagnostics:
            self._runDiagnostics()
        
    def __repr__(self):
        return f'VariantGrouper: {self.bamfile}'
        
    def _getSignalPos(self,vdf):
        #identify positions with signal
        matchOrNa   = ((vdf == '.') | (vdf == ',') | vdf.isna())
        if not self.indels:
            #don't count indels
            matchOrNa |= vdf.apply(lambda pos:pos.str.contains('\+|-')).fillna(False)
        fracVar     = (~matchOrNa).sum(axis=0)/len(vdf)
        if self.log:
            self.log.info(f'Reducing feature space: using {sum(fracVar >= self.minSignal)} positions')
        cols = vdf.columns[fracVar >= self.minSignal]
        if self.hpmask > 0:
            patt = re.compile(rf'([ATGC])\1{{{self.hpmask},}}')
            ref  = pysam.FastaFile(self.refFasta)
            if self.region:
                region      = self.region
                ctg,start,_ = getCoordinates(self.region)
                start -= 1
            else:
                region    = ref.references[0] 
                ctg,start = region,0
            sequence = ref.fetch(region=region)
            mask = list(chain(*([(ctg,p+start) for p in range(m.start()-1,m.end())] 
                                for m in patt.finditer(sequence))))
        else:
            mask = []
        if self.truncate and self.region:
            ctg,start,stop = getCoordinates(self.region)
            rgn  = [(ctg,p) for p in range(start,stop+1)]
            cols = cols[cols.isin(rgn)]
        return cols[~cols.isin(mask)]
    
    def _getMinCount(self):
        return max(self._minReads,ceil(self._minFrac*self.stats['clustered reads']))    
    
    def _checkBam(self,inFile):
        #weak check
        assert inFile.endswith('.bam')
        return inFile

    def _sanitizeValues(self,vdf):
        '''homogenize across strands and drop reads not covering all selected pos'''
        return vdf.apply(self._homogenizeStrands).dropna(axis='index',how='any')
    
    def _makeDf(self):
        if self.log:
            self.log.info(f'Reading alignments from input BAM using {self.nproc} procs')
        if self.nproc == 1:
            df  = self.makeDf(self.bamfile,self.refFasta,self.region,self.truncate)
            bam = pysam.AlignmentFile(self.bamfile) 
            self.stats['total alignments'] = bam.count()
            bam.reset()
            self.stats['primary alignments'] = sum(1 for r in bam if not r.flag & 0x900)
            self.stats['pileup alignments'] = len(df)
            return df            
        else:
            from multiprocessing import Pool
            pool     = Pool(self.nproc)
            result   = []
            for chunk in self.chunkBam():
                callback = self._chunkCallback(chunk,result)
                pool.apply_async(self.makeDf,
                                 args=(chunk,self.refFasta),
                                 callback=callback) 
            pool.close()
            pool.join()
            res = pd.concat(result)
            self.stats['pileup alignments'] = len(res)
            return res

    def _chunkCallback(self,chunk,out):
        def f(res):
            if self.log:
                self.log.debug(f'Finished reading {chunk}')
            out.append(res)
            os.remove(chunk)
            os.remove(f'{chunk}.bai')
        return f

    def chunkBam(self):
        outFmt = f'{self.prefix}chunk{{}}.bam'
        chunk  = 0
        oname  = outFmt.format(chunk)
        chunks = [oname]
        secNsup = 0
        with pysam.AlignmentFile(self.bamfile) as inBam:
            nreads = inBam.count(region=self.region)
            breaks = np.linspace(0,nreads,self.nproc+1,dtype=int)[1:-1]
            outBam = pysam.AlignmentFile(oname,'wb',template=inBam)
            for i,rec in enumerate(inBam.fetch(region=self.region)):
                if rec.flag & 0x900:
                    secNsup += 1
                else:
                    outBam.write(rec)
                if i in breaks:
                    outBam.close()
                    pysam.index(oname)
                    yield oname
                    chunk += 1
                    oname = outFmt.format(chunk)
                    outBam = pysam.AlignmentFile(oname,'wb',template=inBam)
                    chunks.append(oname)
            outBam.close()
            pysam.index(oname)
            yield oname
            self.stats['total alignments'] = i+1
            self.stats['primary alignments'] = i+1 - secNsup
            if self.log:
                self.log.info(f'{i+1} Alignments chunked into {self.nproc} tmp bams; {secNsup} filtered secondary/supplemental alignments')

    def _runDiagnostics(self):
        if self.log:
            self.log.info('Running diagnostics')
        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        for name,df in {'all':self.vTable,'sigvar':self.vTable[self.sigVar.columns]}.items():
            fig,ax = plt.subplots(ncols=2)
            fig.set_figheight(10)
            fig.set_figwidth(20)
            df.isnull().sum().reset_index('contig',drop=True).plot(ax=ax[0])
            ax[0].set_ylabel('n Null Values')
            sns.heatmap(df.isnull(),ax=ax[1],
                        cmap="cubehelix",cbar=False,
                        xticklabels=False,yticklabels=False)
            fig.savefig(f'{self.prefix}{name}_cov.png',format='png',dpi=800)
        return None

    def _makeSigVar(self):
        useCols = self._getSignalPos(self.vTable)
        outdf = self._sanitizeValues(self.vTable[useCols])
        self.stats['clustered reads'] = len(outdf)
        if self.log:
            self.log.debug(f'Removing reads not covering all positions. Input: {len(self.vTable)} Passing: {len(outdf)}')
        return outdf

    def _homogenizeStrands(self,s):
        '''variants are coded by strand (.+2CC vs ,+2cc). This makes them the same'''
        return s.str.upper().str.replace(',','.')

    def getEndpoints(self,vtbl,minCov=1):
        usePos = vtbl.notnull().sum() >= minCov
        return vtbl.columns[usePos][::sum(usePos)-1].get_level_values('pos')

    def split(self,reads):
        '''Returns largest group of reads from position with highest entropy'''
        if self.aggressive:
            maxReads = len(reads) - 1
        else:
            maxReads = len(reads) - self.minCount
        #vTable = self.sigVar.loc[reads]
        vTable = self.sigVar.reindex(reads)
        counts = vTable.apply(pd.Series.value_counts).fillna(0)
        for ch in '.*':
            if ch in counts.index:
                counts.drop(ch,inplace=True)
        counts.drop(columns=counts.columns[counts.sum()==0],inplace=True)
        ent = counts.apply(lambda p: p.sum()*entropy(p.dropna()))\
                    .sort_values(ascending=False)
        for pos in ent.index:
            try:
                vnt = counts[pos].dropna().idxmax()
            except ValueError as e: #empty return
                break
            subset = vTable[vTable[pos]==vnt].index
            if len(subset) >= self.minCount and len(subset) <= maxReads: 
                return subset,pos,vnt
        return None,None,None

class VariantSubCluster(VariantGrouper):
    '''subclass of VariantGrouper with different split method'''
    maxFeatures = 3
    def _rankEntropy(self,counts):
        return counts.apply(lambda p: pd.Series({'score'  :p.sum()*entropy(p.dropna()),
                                                 'entropy':entropy(p.dropna())}))\
                     .T.sort_values('score',ascending=False)

    def _similarity(self,features,entropy):
        sim = pd.DataFrame(0,index=features.index,columns=features.index,dtype=float)
        for pos in features.columns:
            for vnt,subset in features.groupby(pos):
                #might need to change incr function to just 1
                if len(subset) < self.minCount:
                    continue
                sim.loc[subset.index,subset.index] += entropy[pos]
        return sim

    def split(self,reads):
        if self.aggressive:
            maxReads = len(reads) - 1
        else:
            maxReads = len(reads) - self.minCount
        #vTable = self.sigVar.loc[reads]
        vTable = self.sigVar.reindex(reads)
        counts = vTable.apply(pd.Series.value_counts).fillna(0)
        for ch in '.*':
            if ch in counts.index:
                counts.drop(ch,inplace=True)
        counts.drop(columns=counts.columns[counts.sum()==0],inplace=True)
        ent     = self._rankEntropy(counts)
        useCols = ent.index[:self.maxFeatures]
        #ent = counts.apply(lambda p: p.sum()*entropy(p.dropna()))\
        #            .sort_values(ascending=False)    
        #useCols    = ent[ent>=np.percentile(ent,80)].index[:self.maxFeatures]
        if self.log:
            self.log.debug(f'Checking for groups using pos {tuple(useCols)}')
        features   = self.sigVar.reindex(reads)[useCols]     
        similarity = self._similarity(features,ent.entropy)
        spectral   = SpectralClustering(n_clusters=2,affinity='precomputed')
        scaled     = check_symmetric(MinMaxScaler().fit_transform(similarity),raise_warning=False)
        clustv     = spectral.fit_predict(scaled)
        #use group with most non-ref calls
        useClust   = features.groupby(clustv).apply(lambda d:((d!='.').sum()/len(d)).mean()).idxmax()
        size       = sum(clustv==useClust)
        if size >= self.minCount and size <= maxReads: 
            subset = features.index[clustv == useClust]
            var    = features.reindex(subset).apply(pd.Series.value_counts).idxmax().values
            return subset,tuple(useCols),tuple(var)
        else:
            return None,None,None



class sparseDBG:
    def __init__(self,k,collapse=1,ignoreEnds=0,minReads=5,minFrac=0.1,log=None,stats={}):
        self.parser   = seqParser(k,collapse,ignoreEnds=ignoreEnds)
        self.k        = k
        self.collapse = collapse
        self.minReads = minReads
        self.minFrac  = minFrac
        self.log      = log
        self.nodes    = {}
        self.stats    = stats

    def __repr__(self):
        return f'SparseDBG <k:{self.k} r:{self.minReads} f:{self.minFrac}>'
        
    def __getitem__(self,key):
        return self.nodes[key]
    
    def __setitem__(self,key,item):
        self.nodes[key] = item
        
    def loadReads(self,inFile,region=None,minLength=50,maxLength=50000):
        recGen = RecordGenerator(inFile,minLength=minLength,maxLength=maxLength)
        self.name2idx  = recGen.getNameIdx()
        self.readnames = list(self.name2idx.keys())
        nReads         = len(self.readnames)
        self.minCount  = max(ceil(self.minFrac*nReads),self.minReads)
        allNodes       = {}
        if self.log:
            self.log.info('Building debruijn graph')
        for i,rec in enumerate(recGen):
            for kmer in self.parser(rec.sequence):
                nseq = kmer[:-1]
                if nseq in allNodes:
                    node = allNodes[nseq]
                    if node.add(i,kmer):
                        self.nodes[nseq] = node
                else:
                    allNodes[nseq] = Node(readID=i,
                                          kmer=kmer,
                                          minCount=self.minCount)
        #record counts
        counts  = recGen.counter
        secNsup = sum(counts.get(k,0) for k in ['secondary','supplementary'])
        filtlen = sum(counts.get(k,0) for k in counts.keys() if k.startswith('long') or k.startswith('short'))
        prim    = counts['pass']
        self.stats.update({'total alignments'  : prim + secNsup + filtlen,
                           'primary alignments': prim,
                           'pileup alignments' : nReads,
                           'clustered reads'   : nReads})
        if self.log:
            self.log.debug(recGen.report())
                
    def split(self,readnames):
        rIdxs = [self.name2idx[name] for name in readnames] 
        try:
            kmer,node = sorted(filter(itemgetter(1),
                                      ((k,n.subset(rIdxs)) 
                                       for k,n in self.nodes.items())),
                               key=lambda t: t[1].score())[-1]
            e,idx     = node.maxEdge
            reads     = [self.readnames[i] for i in idx]
            return reads,kmer,e
        except IndexError as e:
            return None,None,None
            
class Node:
    def __init__(self,edges=None,readID=None,kmer=None,minCount=0):
        self.isLoop   = False
        self.minCount = minCount
        if edges is None:
            assert readID is not None and kmer is not None,'need input, either (edges,) or (readId,kmer)'
            self.outEdges = defaultdict(set)
            self.add(readID,kmer)
        else:
            self.outEdges = edges
            
    def score(self):
        splits = self.splits()
        return entropy(splits)*sum(splits)
            
    def add(self,readID,kmer):
        wasUsed = self.isUsable
        if sum(readID in edge for edge in self.outEdges.values()):
            self.isLoop = True
        self.outEdges[kmer[-1]].add(readID)
        return not self.isLoop and self.isUsable and not wasUsed
    
    def splits(self,readIDs=None):
        if readIDs:
            edges = self.intersect(readIDs)
        else:
            edges = self.outEdges
        return list(map(len,edges.values()))
    
    def intersect(self,readIDs):
        edges = {i:s.intersection(readIDs) for i,s in self.outEdges.items()}
        return edges
    
    def subset(self,readIDs):
        subnode = Node(edges=self.intersect(readIDs),
                       minCount=self.minCount)
        return subnode if subnode.isUsable else None
    
    @property
    def maxEdge(self):
        return sorted(self.outEdges.items(),key=lambda e:len(e[1]))[-1]
    
    @property
    def isUsable(self):
        return sum(s>self.minCount for s in self.splits()) >= 2
    
    def __repr__(self):
        return f'Node: {self.splits()}'

def ident(x): return x

class seqParser:
    def __init__(self,k=11,collapseHP=1,minimizer=0,ignoreEnds=0):
        self.k         = k
        self.transform = hpCollapse(collapseHP) if collapseHP >= 1 else ident
        self.minim     = self.getMinimizer(minimizer) if minimizer>0 else ident
        self.start     = ignoreEnds
        self.end       = -ignoreEnds if ignoreEnds else 1000000 #really big to get everything
    def __call__(self, seq):
        s = self.transform(seq[self.start:self.end])
        for i in range(len(s)-self.k+1):
            yield self.minim(s[i:i+self.k])
    def getMinimizer(self,m=6):
        def minimizer(seq):
            return sorted(seq[i:i+m] for i in range(0,len(seq)-m+1))[0]
        return minimizer

class Splitter_Error(Exception):
    pass

