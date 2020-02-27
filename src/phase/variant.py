import pysam
from math import ceil
import pandas as pd
import numpy as np
from collections import Counter

class VariantGrouper:
    def __init__(self,bamfile,refFasta,sampleName,
                 region=None,truncate=False,
                 minCov=10,minFrac=0.1,
                 minReads=10,minSpan=0.9,
                 minSignal=0.05,flagFilter=0x900,
                 minBaseQ=0):
        self.bamfile    = bamfile      #bam file
        self.refFasta   = refFasta     #ref
        self.region     = region       #region eg X:12345-678900
        self.truncate   = truncate     #truncate pileup to region (pysam pileup kwarg)
        self.minCov     = minCov       #min pos coverage for inclusion
        self._minFrac   = minFrac      #min fraction to split off
        self._minReads  = minReads     #abs min reads to split off
        self.minSpan    = minSpan      #not implem yet
        self.minSignal  = minSignal    #min frac of reads that are diff from ref to consider
        self.minBaseQ   = minBaseQ     #min baseQV, to pysam.pileup engine
        self.flagFilter = flagFilter   #hex filters, passed to pysam.pileup engine
        self.sampleName = sampleName   #sample name
        self.vTree      = None
        self.pending    = []
        
    def __repr__(self):
        if self.vTree:
            return str(self.vTree)
        else:
            return f'VariantGrouper: {self.sampleName}:{self.bamfile}'
        
    def run(self):
        #read variants in pileup
        print('Reading Alignments')
        self.varDf   = self._makeDf()
        self.inReads = len(self.varDf)
        
        #TODO filter reads by min span?? 
        #(reads not covering selected pos are
        #removed below)
        
        #identify positions with signal
        print('Reducing Features')
        self.usePos = self._getSignalPos(self.varDf)
        
        #homogenize across strands and drop reads not covering all selected pos
        self.sigVar = self._sanitizeValues(self.varDf[self.usePos])
        
        #reads left covering all selected pos
        self.NclustReads = len(self.sigVar) 
        
        #set smallest cluster allowed, det by frac or abs min
        self.minCount    = max(ceil(self._minFrac*self.NclustReads),self._minReads)
        
        #initialize starting conditions
        print('Initializing Tree')
        self._initGroups()
        #split groups
        print('Splitting Groups')
        self._splitGroups()
        # identify residual node and dump in noise bin
        print('Cleaning Up')
        try:
            residualNode = list(filter(self.vTree.isResidual,self.vTree.nodes))[0]
            #self.rescue(residualNode)
        except IndexError as e:
            residualNode = None
            pass
        self.residual = residualNode
        self._dumpResidual()
        #get clusters
        self._getClusters()
        
    def _getSignalPos(self,vdf):
        #identify positions with signal
        matchOrNa   = ((vdf == '.') | (vdf == ',') | vdf.isna())
        fracVar     = (~matchOrNa).sum(axis=0)/len(vdf)
        return vdf.columns[fracVar >= self.minSignal]
    
    def _sanitizeValues(self,vdf):
        '''homogenize across strands and drop reads not covering all selected pos'''
        return vdf.apply(self._homogenizeStrands).dropna(axis='index',how='any')
        
    def _makeDf(self):
        '''DF with read names as index and ref positions as columns and variants as elements'''
        bam = pysam.AlignmentFile(self.bamfile,'r')
        ref = pysam.FastaFile(self.refFasta)
        return pd.DataFrame({column.reference_pos : dict(zip(column.get_query_names(),
                                                             column.get_query_sequences(mark_matches=True,
                                                                                        add_indels=True)))
                             for column in bam.pileup(flag_filter=self.flagFilter,
                                                      fastafile=ref,region=self.region,
                                                      truncate=self.truncate,
                                                      min_base_quality=self.minBaseQ)
                             if len(column.get_query_names()) >= self.minCov})
    def _homogenizeStrands(self,s):
        '''variants are coded by strand (.+2CC vs ,+2cc). This makes them the same'''
        return s.str.upper().str.replace(',','.')
    
    def _initGroups(self):
        '''id recall and start tree'''
        #init labeler
        self.labeler = self.indexer()
        #all reads
        rootLbl = next(self.labeler)
        root = vCluster(self.sigVar.index,
                        rootLbl,
                        split=self.sampleName)
        #group tree
        self.vTree = vTree(root)
        #reads matching reference at all selected pos
        refcall = self.sigVar.index[(self.sigVar == '.').all(axis=1)]
        #if refcall reads < minCount --> add to noise
        if len(refcall) < self.minCount:
            self.vTree.noise.update(refcall)
        #else it is first split group
        else:
            clust1  = vCluster(refcall,
                               next(self.labeler),
                               parent=rootLbl,
                               split=('all','refcall'),
                               pending=False)
            self.vTree.append(clust1)
        #reads with vars in selected pos
        pending = self.sigVar.index.difference(refcall)
        clust2  = vCluster(pending,
                           next(self.labeler),
                           parent=rootLbl,
                           split=(-1,'variable'))
        self.vTree.append(clust2)
        self.vTree[rootLbl].pending = False
        self._updatePending()
        return
    
    def _updatePending(self):
        self.pending = sorted([node.label 
                               for node in self.vTree.nodes.values() 
                               if node.pending],
                              reverse=True)
    
    def _splitGroups(self,vTable=None):
        if vTable is None:
            vTable = self.sigVar
        #func to det if big enough to split
        canSplit = (lambda t: len(t) >= 2*self.minCount)
        while len(self.pending):
            label   = self.pending.pop()
            self.vTree[label].pending = False
            testSet = self.vTree[label].reads
            subset,pos,vnt = self.vGen(vTable.loc[testSet])
            if subset is None:
                continue
                #self.vTree[label].pending = False
            else:
                #joined by variant
                newGroup = vCluster(subset,
                                    next(self.labeler),
                                    parent =label,
                                    split  =(pos,vnt),
                                    pending=canSplit(subset))
                self.vTree.append(newGroup)
                #leftovers
                complement = testSet.difference(subset)
                remainder  = vCluster(complement,
                                      next(self.labeler),
                                      parent =label,
                                      split  =(pos,'.'),
                                      pending=canSplit(complement))
                self.vTree.append(remainder)
            self._updatePending()
        return

    def _dumpResidual(self):
        '''put reads from "residual" node in noise bin'''
        node = self.vTree[self.residual]
        self.vTree.noise.update(node.reads)
        node.setNoise()
    
    def _getClusters(self):
        '''map of node label -> cluster'''
        #clusters sorted by size, descending
        sortedLeaves = filter(lambda n: n!= self.residual,
                              sorted(self.vTree.leaves,
                                     key=lambda node:-len(self.vTree[node])))
        self.node2cluster = {node:idx for idx,node in enumerate(sortedLeaves)}
        #update nodes
        for node,clust in self.node2cluster.items():
            self.vTree[node].set_cluster(clust)
        return
    
    def _getPlurality(self,vTable):
        return vTable.apply(pd.Series.value_counts)\
                     .fillna(0).idxmax()
    
    def summary(self):
        '''return summary of results in dataframe format'''
        #most common call by cluster/group
        plurality = self.sigVar.groupby(self.clusterMap).apply(self._getPlurality)
        #remove columns with all refcall (.) or del-continue (*)
        varCols   = plurality.columns[~plurality.isin(['.','*']).all(axis=0)]
        #new table
        summary              = plurality[varCols].copy()
        countMap             = Counter(self.clusterMap.values())
        summary['nReads']    = summary.index.map(countMap) 
        summary.index        = np.where(summary.index >= 0,
                                        'cluster' + summary.index.astype(str),
                                        'ResidualNoise')
        summary.index.name   = 'Group'
        summary.columns.name = 'RefPos'
        return summary

    def rescue(self,label):
        '''go back to original and try to identify any significant pos that will group'''
        node   = self.vTree[label]
        print(f'Trying to rescue/split {node}')
        vTable = self.varDf.loc[node.reads]
        newPos = self._getSignalPos(vTable)
        if len(newPos) == 0:
            return None
        newSigTable  = self._sanitizeValues(vTable[newPos])
        node.pending = True
        self._updatePending()
        self._splitGroups(newSigTable)
        return
    
    def vGen(self,vTable):
        '''Returns smallest group of reads with shared variant given minCount'''
        maxReads = len(vTable) - self.minCount
        counts = vTable.apply(pd.Series.value_counts).fillna(0)
        for ch in '.*':
            if ch in counts.index:
                counts.drop(ch,inplace=True)
        shared = pd.concat([counts.idxmax().rename('var'),
                            counts.max().rename('cnt')],
                           axis=1).sort_values('cnt')['var']
        for pos,vnt in shared.items():
            subset = vTable[vTable[pos]==vnt].index
            if len(subset) >= self.minCount and len(subset) <= maxReads: 
                return subset,pos,vnt
        return None,None,None
    
    @property
    def clusterMap(self):
        '''dict of {readname:cluster}'''
        clustMap = {name:-1 for name in self.vTree.noise}
        for node,clust in self.node2cluster.items():
            clustMap.update({name:clust for name in self.vTree[node].reads})
        return clustMap
        
    def indexer(self,start=0):
        i=start
        while True:
            yield i
            i+=1

class vCluster:
    def __init__(self,reads,label,parent=None,split=None,pending=True):
        self.label    = label
        self.reads    = reads
        self.parent   = parent
        self.split    = split
        self.children = []
        self.pending  = pending
        #private
        self._isLeaf   = False
        self._cluster  = None
        self._isNoise  = False
        #pending = potential to split
    def set_cluster(self,cluster):
        self._isLeaf  = True
        self._cluster = cluster
    def setNoise(self):
        self._isNoise = True
    def clusterName(self):
        assert self._cluster is not None
        return f'cluster{self._cluster}_numreads{len(self.reads)}'
    def __len__(self):
        return len(self.reads)
    def __repr__(self):
        outstr = f'Node {self.label} ({len(self.reads)} reads)'
        if self._isLeaf and self._cluster is not None:
            outstr += f' Cluster {self._cluster}'
        elif self._isNoise:
            outstr += ' ResidualNoise'
        return outstr 
        
class vTree:
    def __init__(self,root):
        if not root.parent is None:
            raise VariantPartition_Error('root parent must = None')
        self.root   = root
        self.label  = root.label
        self.size   = len(root)
        self.leaves = [root.label]
        self.nodes  = {root.label:root}  
        self.noise  = set()
    def append(self,other):
        if other.parent is None or other.parent not in self.nodes:
            raise VariantPartition_Error('invalid parent')
        if other.split is None:
            raise VariantPartition_Error('please set split value')
        self.nodes[other.label] = other
        parent = self.nodes[other.parent]
        if len(parent.children) == 0:
            #remove parent as leaf
            self.leaves.pop(self.leaves.index(parent.label))
        #add this leaf
        self.leaves.append(other.label)
        parent.children.append(other.label)
    def getSplits(self,label):
        if label not in self.nodes:
            raise VariantPartition_Error('invalid label')
        splits = []
        leaf   = label
        while leaf:
            vclust = self.nodes[leaf]
            splits.insert(0,vclust.split)
            leaf = vclust.parent
        return splits
    def isResidual(self,label):
        '''identify residual group -- leftovers at all branches'''
        spl = self.getSplits(label)[1:]
        if len(spl) == 0:
            return False
        else:
            return np.all([s=='.' for p,s in spl]) and len(self[label].children)==0
    def _splitStr(self,label):
        '''join splits for output'''
        splits = self.getSplits(label)
        splits.insert(0,'all')
        return ' -> '.join(map(str,splits))
    def __getitem__(self,node):
        return self.nodes[node]
    def __repr__(self):
        leaves = [f'{str(self[lbl])} : {self._splitStr(lbl)}' for lbl in sorted(self.leaves)]
        return '\n'.join(leaves)        

class VariantPartition_Error(Exception):
    pass
