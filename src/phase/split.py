import pysam
from math import ceil
import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import entropy
from operator import itemgetter
from collections import defaultdict
from .utils import hpCollapse,RecordGenerator

class VariantGrouper:
    def __init__(self,inFile,refFasta,
                 region=None,truncate=False,
                 minCov=1,minFrac=0.1,
                 minReads=10,minSpan=0.9,
                 minSignal=0.05,flagFilter=0x900,
                 aggressive=False,vTable=None):
        self.bamfile    = self._checkBam(inFile)   #bam file
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
        self.vTable     = vTable if vTable is not None else self._makeDf()
        self.sigVar     = self._makeSigVar()
        self.nReads     = len(self.sigVar)
        self.minCount   = self._getMinCount()
        self.readnames  = self.vTable.index
        
    def __repr__(self):
        return f'VariantGrouper: {self.bamfile}'
        
    def _getSignalPos(self,vdf):
        #identify positions with signal
        matchOrNa   = ((vdf == '.') | (vdf == ',') | vdf.isna())
        fracVar     = (~matchOrNa).sum(axis=0)/len(vdf)
        return vdf.columns[fracVar >= self.minSignal]
    
    def _getMinCount(self):
        return max(self._minReads,ceil(self._minFrac*self.nReads))    
    
    def _checkBam(self,inFile):
        #weak check
        assert inFile.endswith('.bam')
        return inFile

    def _sanitizeValues(self,vdf):
        '''homogenize across strands and drop reads not covering all selected pos'''
        return vdf.apply(self._homogenizeStrands).dropna(axis='index',how='any')
        
    def _makeDf(self):
        '''DF with read names as index and ref positions as columns and variants as elements'''
        bam = pysam.AlignmentFile(self.bamfile,'r')
        ref = pysam.FastaFile(self.refFasta)
        df  = pd.DataFrame({(column.reference_name,
                             column.reference_pos)  : dict(zip(column.get_query_names(),
                                                            column.get_query_sequences(mark_matches=True,
                                                                                       add_indels=True)))
                            for column in bam.pileup(flag_filter=self.flagFilter,
                                                     fastafile=ref,region=self.region,
                                                     truncate=self.truncate,
                                                     min_base_quality=0)
                            if len(column.get_query_names()) >= self.minCov})
        df.columns.names = ('contig','pos')
        return df

    def _makeSigVar(self):
        useCols = self._getSignalPos(self.vTable)
        return self._sanitizeValues(self.vTable[useCols])

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
        vTable = self.sigVar.reindex(reads)
        counts = vTable.apply(pd.Series.value_counts).fillna(0)
        for ch in '.*':
            if ch in counts.index:
                counts.drop(ch,inplace=True)
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

class sparseDBG:
    def __init__(self,k,collapse=1,ignoreEnds=0,minReads=5,minFrac=0.1):
        self.parser   = seqParser(k,collapse,ignoreEnds=ignoreEnds)
        self.k        = k
        self.collapse = collapse
        self.minReads = minReads
        self.minFrac  = minFrac
        self.nodes    = {}

    def __repr__(self):
        return f'SparseDBG <k:{self.k} r:{self.minReads} f:{self.minFrac}>'
        
    def __getitem__(self,key):
        return self.nodes[key]
    
    def __setitem__(self,key,item):
        self.nodes[key] = item
        
    def loadReads(self,inFile,region=None,minLength=50,maxLength=50000):
        recGen = RecordGenerator(inFile,minLength=minLength,maxLength=maxLength)
        self.name2idx  = {rec.name:i for i,rec in enumerate(recGen)}
        self.nReads    = len(self.readnames)
        self.minCount  = max(ceil(self.minFrac*self.nReads),self.minReads)
        allNodes       = {}
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
                
    @property
    def readnames(self):
        return list(self.name2idx.keys())
    
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
        #self.transform = hpCollapse if collapseHP else ident
        self.transform = hpCollapse(collapseHP) if collapseHP >= 1 else ident
        #self.minim     = getMinimizer(minimizer) if minimizer>0 else ident
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

