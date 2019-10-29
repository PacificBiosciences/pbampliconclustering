from sklearn.cluster import DBSCAN,\
                            OPTICS, \
                            AgglomerativeClustering, \
                            AffinityPropagation, \
                            KMeans, \
                            MeanShift
import numpy as np
from collections import Counter
import json

class ClusterModel:
    '''
    defaults: { modelKwargs     : defaultValue }
    pmap    : { CommandLineArgs : modelKwargs } 
    
    returns : parameterized model instance with fit method 
    '''
    MODEL     = None
    defaults  = {}
    pmap      = {}
    exemplars = None

    def __init__(self,args):
        #load CL args first to override defaults
        self.defaults.update({mparam:getattr(args,aparam) 
                              for aparam,mparam in self.pmap.items() 
                              if getattr(args,aparam)})
        #override with param file
        if args.params:
            with open(args.params) as configFile:
                config = json.load(configFile)
            self.defaults.update(config)
        self.model = self.MODEL(**self.defaults)
    def fit(self,X):
        return self.model.fit(X)
    def __repr__(self):
        name          = self.MODEL.__name__
        dashes        = ''.join(['-']*((40 - len(name)) // 2))
        pmap,defaults = ['\n'.join([f'{arg:>20} => {kw}' 
                                    for arg,kw in getattr(self,var).items()]) 
                         for var in ['pmap','defaults']]
        return f'{dashes}{name}{dashes}\nArgMap:\n{pmap}\nDefaults:\n{defaults}\n'

class ClusterModel_wNoise(ClusterModel):
    '''
    Subclass for models without a min cluster size.
    Re-labels small clusters < minReads to noise (-1)
    '''
    def __init__(self,args):
        self.minCnt = args.minReads
        self._noiseLabels = []
        super().__init__(args)
    def fit(self,X):
        res = self.model.fit(X)
        for val,count in Counter(res.labels_).items():
            if count < self.minCnt:
                self._noiseLabels.extend(np.where(res.labels_==val)[0])
                res.labels_[res.labels_==val] = -1
        return res

class Dbscan(ClusterModel):
    MODEL    = DBSCAN
    defaults = {'eps'        : 0.01,
                'min_samples': 3,
                'metric'     : 'euclidean',
                'n_jobs'     : 2}
    pmap     = {'eps'        : 'eps',
                'minReads'   : 'min_samples',
                'njobs'      : 'n_jobs'}

class Optics(ClusterModel):
    MODEL    = OPTICS
    defaults = {'max_eps'    : np.inf,
                'min_samples': 3,
                'n_jobs'     : 1,
                'metric'     : 'l2',
                'xi'         : 0.05}
    pmap     = {'eps'        : 'max_eps',
                'minReads'   : 'min_samples',
                'njobs'      : 'n_jobs',
                'normalize'  : 'metric'}

class Kmeans(ClusterModel):
    MODEL    = KMeans
    defaults = {'n_clusters'  : 2,
                'max_iter'    : 300,
                'tol'         : 1e-4,
                'random_state': None,
                'n_jobs'      : 1}
    pmap     = {'eps'         : 'tol',
                'njobs'       : 'n_jobs'}                

class Aggcluster(ClusterModel_wNoise):
    MODEL    = AgglomerativeClustering
    defaults = {'affinity'          : 'euclidean',
                'linkage'           : 'ward',
                'compute_full_tree' : True,
                'distance_threshold': 0.01,
                'n_clusters'        : None}
    pmap     = {'eps':'distance_threshold'}

class Affprop(ClusterModel_wNoise):
    MODEL    = AffinityPropagation
    defaults = {'preference' : None, #median of inputs
                'damping'    : 0.5}
    pmap     = {'eps'     : 'preference'}

    #def __init__(self,args):
    #    if args.eps and (args.eps < 0.5 or args.eps >1):
    #        raise Clustering_Exception('Damping (-e) must be in [0.5-1] for AffinityPropagation')
    #    super().__init__(args)

class Meanshift(ClusterModel):
    MODEL    = MeanShift
    defaults = {'bandwidth'   : None, #estimate from data
                'bin_seeding' : True,
                'min_bin_freq': 3,
                'cluster_all' : False,
                'n_jobs'      : 2}
    pmap     = {'eps'         : 'bandwidth',
                'minReads'    : 'min_bin_freq',
                'njobs'       : 'n_jobs'}

MODELS = {'dbscan'    : Dbscan,
          'optics'    : Optics,
          'aggcluster': Aggcluster,
          'affprop'   : Affprop,
          'meanshift' : Meanshift,
          'kmeans'    : Kmeans}

def showModels(args):
    models = [MODELS[args.model]] if args.model else MODELS.values()

    for m in models:
        print(m.__repr__(m))

class Clustering_Exception(Exception):
    pass

