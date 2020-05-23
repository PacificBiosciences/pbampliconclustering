import pysam
import pandas as pd
import numpy as np
from scipy.stats import entropy
from collections import Counter

FIGFORMAT= 'pdf'
MINCOUNT = 3 #absolute minimum shared min variants

class VariantCaller:
    def __init__(self,sigVar,minCount,reference,clusterMap,endPoints,log=None):
        self.vTable     = sigVar
        self.minCount   = minCount
        self.reference  = pysam.FastaFile(reference)
        self.clusterMap = clusterMap
        self.endPoints  = endPoints
        self.readCounts = Counter(self.clusterMap.values())
        self.totalReads = len(self.vTable)
        self.log        = log
    
    def _getVarCounts(self,vtbl,fillna=True):
        res = vtbl.apply(pd.Series.value_counts)
        return res.fillna(0) if fillna else res
    
    def _getPlurality(self,vtbl):
        return self._getVarCounts(vtbl).idxmax()
    
    def plurality(self,grpMap=None):
        gmap = self.clusterMap if grpMap is None else grpMap
        plrty = self.vTable.groupby(gmap).apply(self._getPlurality)
        vcols = plrty.columns[~plrty.isin(['.','*']).all(axis=0)]
        plrty.index.name = 'cluster'
        return plrty[vcols]
    
    @property
    def variantGroupMap(self):
        varPos    = self._sigCounts().index
        #catch the case where all reads are exactly the reference in the sig positions
        if len(varPos) == 0:
            return {r:0 for r in self.vTable.index}
        asTuple   = self.vTable[varPos].apply(tuple,axis=1)
        refCall   = (('.'),)*len(varPos)
        varGroups = self.vTable.groupby(asTuple)
        vkeys     = sorted(varGroups.groups.keys(),key=lambda k:-len(varGroups.get_group(k)))
        clust     = 0
        if refCall in vkeys:
            varGrpMap = {r:clust for r in varGroups.get_group(refCall).index}
            vkeys.pop(vkeys.index(refCall))
        else:
            varGrpMap = {}
        clust += 1
        for vkey in vkeys:
            size = len(varGroups.get_group(vkey))
            if size >= MINCOUNT:
                cnumber = clust
                clust  += 1
            else:
                if self.log:
                    self.log.debug(f'Assigning variant tuple {vkey} as noise (size={size})')
                cnumber = -1 #noise
            varGrpMap.update({r:cnumber for r in varGroups.get_group(vkey).index})
        return varGrpMap

    @property
    def byCluster(self):
        return self.vTable.groupby(self.clusterMap)
    
    @property
    def entropy(self):
        ent = self.byCluster.apply(lambda d: self._getVarCounts(d).apply(entropy))
        ent.index.name = 'cluster'
        return ent.assign(meanEntropy=ent.mean(axis=1))
    
    @property
    def variants(self):
        grp = self.plurality().unstack().groupby(level=(0,1))
        vtbl = grp.apply(lambda d:d.drop_duplicates()\
                                   .sort_values()\
                                   .reset_index(drop=True)\
                                   .rename('variant'))
        if len(vtbl):
            vtbl.index.rename('idx',level=-1,inplace=True)
        return vtbl
   
    def _sigCounts(self):
        allCounts = self._getVarCounts(self.vTable,fillna=False)
        hasDel = '*' in allCounts.index
        drop = ['.','*'] if hasDel else '.' 
        try:
            usePos    = allCounts.drop(drop).max() >= self.minCount
            useVars   = allCounts.drop('*' if hasDel else []).max(axis=1) >= self.minCount
        except KeyError as e:
            usePos,useVars = [False],[False]
        #check for null output
        if sum(usePos) == 0 or sum(useVars) == 0:
            midx = pd.MultiIndex(levels=[[],[]],codes=[[],[]],names=allCounts.columns.names)
            return pd.DataFrame(columns=['-1'], index=midx)
        else:
            return allCounts.reindex(index=useVars[useVars].index,columns=usePos[usePos].index).T
    
    @property
    def variantFractions(self):
        sigCounts = self._sigCounts()
        sigCounts['noise'] = self.totalReads - sigCounts.sum(axis=1)
        return (sigCounts / self.totalReads).fillna(0)
    
    def alleleFrequency(self,grpMap=None):
        gmap    = self.clusterMap if grpMap is None else grpMap
        readcnt = Counter(gmap.values())
        summary = self.plurality(gmap)
        summary['nReads']    = summary.index.map(readcnt)
        summary['frequency'] = summary.nReads / self.totalReads
        return summary

    def _updateSequence(self,row,start,end):
        bases    = 'ATGC'
        prevPos  = start
        ctg      = self.variants.index.get_level_values(0)[0]
        for (ctg,pos),vnt in row[row!='.'].items():
            refseq = self.reference.fetch(reference=ctg,start=prevPos,end=pos+1)
            if vnt in bases: #snp
                prevPos = pos + 1
                yield refseq[:-1] + vnt
            elif vnt.startswith('.+'): #insert
                prevPos = pos + 1
                ins = ''.join(b for b in vnt if b in bases)
                yield refseq + ins
            elif vnt.startswith('.-'): #deletion
                delsize = sum(1 for b in vnt if b in bases)
                prevPos = pos + delsize + 1
                yield refseq
            elif vnt == '*': #del continue
                continue
            else:
                raise Caller_Error(f'unknown variant {vnt}')
        yield self.reference.fetch(reference=ctg,start=prevPos,end=end)
        
    def draftConsensus(self):
        try:
            return {clust: ''.join(self._updateSequence(row,*self.endPoints[clust]))
                    for clust,row in self.plurality().iterrows()}
        except IndexError as e:
            return {}

    def entropyPlot(self,savename):
        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        ent = self.entropy
        f,ax = plt.subplots()
        #catch when everything is identical
        if not np.all(ent.isna().values):
            sns.heatmap(self.entropy,ax=ax,cbar_kws={'label': 'Entropy'})
            #plt.tight_layout()
        else:
            pass
        f.savefig(savename,dpi=400,format=FIGFORMAT)
        return f    


class Caller_Error(Exception):
    pass
