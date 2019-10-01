from sklearn.cluster import DBSCAN,\
                            OPTICS, \
                            AgglomerativeClustering, \
                            AffinityPropagation
import numpy as np
from collections import Counter

class ClusterModel:
    defaults = {}
    pmap     = {}
    MODEL    = None

    def __init__(self,args):
        self.defaults.update({mparam:getattr(args,aparam) for aparam,mparam in self.pmap.items() 
                              if getattr(args,aparam)})
        self.model = self.MODEL(**self.defaults)

    def fit(self,X):
        return self.model.fit(X)

class ClusterModel_wNoise(ClusterModel):
    def __init__(self,args):
        self.minCnt = args.minReads
        super().__init__(args)
    def fit(self,X):
        res = self.model.fit(X)
        for val,count in Counter(res.labels_).items():
            if count < self.minCnt:
                res.labels_[res.labels_==val] = -1
        return res

class Dbscan(ClusterModel):
    MODEL    = DBSCAN
    defaults = {'eps'        : 0.01,
                'min_samples': 3}

    pmap     = {'eps'        : 'eps',
                'minReads'  : 'min_samples'}

class Optics(ClusterModel):
    MODEL    = OPTICS
    defaults = {'max_eps'    : np.inf,
                'min_samples': 3,
                'n_jobs'     : 1,
                'metric'     : 'l2',
                'xi'         : 0.1}

    pmap     = {'eps'        : 'max_eps',
                'minReads'   : 'min_samples',
                'njobs'      : 'n_jobs',
                'normalize'  : 'metric'}

class Aggcluster(ClusterModel_wNoise):
    MODEL    = AgglomerativeClustering
    defaults = {'affinity'          : 'euclidean',
                'compute_full_tree' : True,
                'distance_threshold': 0.01,
                'n_clusters'        : None}

    pmap     = {'eps':'distance_threshold'}

class Affprop(ClusterModel_wNoise):
    MODEL    = AffinityPropagation
    defaults = {'damping' : 0.5}
    pmap     = {'eps'     : 'damping'}

    def __init__(self,args):
        if args.eps and (args.eps < 0.5 or args.eps >1):
            raise Clustering_Exception('Damping (-e) must be in [0.5-1] for AffinityPropagation')
        super().__init__(args)

MODELS = {'dbscan'    : Dbscan,
          'optics'    : Optics,
          'aggcluster': Aggcluster,
          'affprop'   : Affprop}

class Clustering_Exception(Exception):
    pass

