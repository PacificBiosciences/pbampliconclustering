import re
import pandas as pd

def clusterName(vals):
    '''vals: tuple of (cluster,numreads)'''
    return f'cluster{vals[0]}_numreads{vals[1]}'

def getCluster(name):
    return int(re.search('cluster(\d+)',name).group(1))

def readClusterFile(clustfile,nFields=3):
    '''
    nFields: use n /-separated fields from input read names in result.
             nFields = 0 -> all fields
    '''
    res  = {}
    name,cluster = None,None
    with open(clustfile) as f:
        for line in f:
            if line.startswith('>'):
                cluster = getCluster(line[1:])
            else:
                fields = line.split('/') 
                read   = '/'.join(fields[:nFields]) if nFields else line.strip()
                try: #get subset of sequence if fourth field has {start}_{stop}
                    start,stop = map(int,fields[3].split('_'))
                except:
                    start,stop = None,None
                res[read] = {'cluster':cluster,'start':start,'stop':stop}
    return pd.DataFrame(res.values(),index=res.keys())

class Cluster_Exception(Exception):
    pass
