import pysam
from math import ceil
import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import entropy

class Phaser:
    def __init__(self,splitter,sampleName,
                 aggressive=False,log=None):
        self.sampleName = sampleName   #sample name
        self.aggressive = aggressive   #aggressively separate groups
        self.splitter   = splitter     #split generator
        self.log        = log
        #################
        self.vTree      = None
        self.pending    = []
        
    def __repr__(self):
        if self.vTree:
            return str(self.vTree)
        else:
            return f'Phaser: {self.sampleName}'

    def run(self):
        #initialize starting conditions
        msg = 'Initializing Tree'
        if self.log:
            self.log.info('Initializing Tree')
        self._initGroups()
        #split groups
        if self.log:
            self.log.info('Splitting Groups')
        self._splitGroups()
        # identify residual node and dump in noise bin
        if self.log:
            self.log.info('Cleaning Up')
        #self._dumpResidual()
        #get clusters
        self._getClusters()
        
    def _initGroups(self):
        '''start splitting tree'''
        #init labeler
        self.labeler = self.indexer()
        #all reads
        rootLbl = next(self.labeler)
        rootGrp = self.splitter.readnames
        root = vCluster(rootGrp,
                        rootLbl,
                        split=self.sampleName,
                        pending=-True)
        #group tree
        if self.log:
            self.log.debug(f'Starting vTree with {root}')
        self.vTree = vTree(root)
        #reads matching reference at all selected pos
        #refcall = self.sigVar.index[(self.sigVar == '.').all(axis=1)]
        self._updatePending()
        return
    
    def _updatePending(self):
        self.pending = sorted([node.label 
                               for node in self.vTree.nodes.values() 
                               if node.pending],
                              reverse=True)
        if self.log:
            self.log.debug(f'Pending nodes: {self.pending}')
    
    def _splitGroups(self):
        #func to det if big enough to split
        if self.aggressive:
            canSplit = (lambda t: len(t) > self.splitter.minCount)
        else:
            canSplit = (lambda t: len(t) >= 2*self.splitter.minCount)
        while len(self.pending):
            label   = self.pending.pop()
            self.vTree[label].pending = False
            testSet = set(self.vTree[label].reads)
            subset,pos,vnt = self.splitter.split(testSet)
            if self.log and subset is not None:
                self.log.debug(f'Attempting to split {self.vTree[label]}; split: {len(subset)},{pos},{vnt}')
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

    def _isRefCall(self,lbl):
        plrty = self.splitter.sigVar.reindex(self.vTree[lbl].reads)\
                             .apply(pd.Series.value_counts).fillna(0).idxmax()
        return (plrty.astype(str)== '.').all()

    def _getClusters(self):
        '''map of node label -> cluster'''
        #check for refcall -- all '.' variants for alignment method.  passes through for dbg meth
        offset = 0
        self.node2cluster = {}
        isRefcall = self._isRefCall if hasattr(self.splitter,'refFasta') else (lambda lbl: False)
        for lbl in self.vTree.leaves:
            leaf = self.vTree[lbl]
            if len(leaf) < self.splitter.minCount:
                #dump in noise bin
                self.vTree.noise.update(leaf.reads)
                leaf.setNoise()
            elif isRefcall(lbl):
                #set as first
                leaf.setRefCall()
                self.node2cluster[lbl] = offset
                offset += 1
        #check if refcall was found, else increase offset if not found but there is one
        if offset == 0 and hasattr(self.splitter,'refFasta') and not self.splitter.refFasta is None:
            offset += 1
        #clusters sorted by size, descending
        sortedLeaves = sorted(filter(lambda n: not (self.vTree[n]._isNoise or self.vTree[n]._isRefCall),
                                     self.vTree.leaves),
                              key=lambda n:-len(self.vTree[n]))
        self.node2cluster.update({node : idx + offset 
                                  for idx,node in enumerate(sortedLeaves)
                                  if not (self.vTree[node]._isNoise or self.vTree[node]._isRefCall)})
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
        self._isLeaf    = False
        self._cluster   = None
        self._isNoise   = False
        self._isRefCall = False
    def set_cluster(self,cluster):
        self._isLeaf  = True
        self._cluster = cluster
    def setNoise(self):
        self._isNoise = True
    def setRefCall(self):
        self._isRefCall = True
        self._cluster   = 0
    def clusterName(self):
        assert self._cluster is not None
        return f'cluster{self._cluster}_numreads{len(self.reads)}'
    def __len__(self):
        return len(self.reads)
    def __repr__(self):
        outstr = f'Node {self.label} ({len(self.reads)} reads)'
        if self._isNoise:
            outstr += ' Noise'
        elif self._isRefCall:
            outstr += f' Cluster {self._cluster} (refcall)'
        elif self._isLeaf and self._cluster is not None:
            outstr += f' Cluster {self._cluster}'
        return outstr 
        
class vTree:
    def __init__(self,root):
        if not root.parent is None:
            raise Phaser_Error('root parent must = None')
        self.root   = root
        self.label  = root.label
        self.size   = len(root)
        self.leaves = [root.label]
        self.nodes  = {root.label:root}  
        self.noise  = set()
    def append(self,other):
        if other.parent is None or other.parent not in self.nodes:
            raise Phaser_Error('invalid parent')
        if other.split is None:
            raise Phaser_Error('please set split value')
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
            raise Phaser_Error('invalid label')
        splits = []
        leaf   = label
        while leaf:
            vclust = self.nodes[leaf]
            splits.insert(0,vclust.split)
            leaf = vclust.parent
        return splits
    def isRefcall(self,label):
        '''identify refcall group -- reads that did not split out by any variants'''
        spl = self.getSplits(label)
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
    def __iter__(self):
        for lbl,node in self.nodes.items():
            yield node
    def __repr__(self):
        leaves = [f'{str(self[lbl])} : {self._splitStr(lbl)}' for lbl in sorted(self.leaves)]
        return '\n'.join(leaves)        

class Phaser_Error(Exception):
    pass
