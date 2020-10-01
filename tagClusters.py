from src.utils.bam import addHPtag
import pandas as pd
import sys,argparse

parser = argparse.ArgumentParser(description='Add HP cluster tag and YC color tag to reads given readname-->cluster(int) pairing. Default setting for use with read_info.txt output from hifi amplicon clustering')
parser.add_argument('inbam',help='bam contaning hifi reads to tag',
                    type=str)
parser.add_argument('readInfo',help='read_info.txt output from hifi amplicon cluster',
                    type=str)
parser.add_argument('outbam',help='output bam',type=str)
parser.add_argument('--names',help='zero-indexed column number with read names. default 0',
                    type=int,default=0)
parser.add_argument('--cluster',help='zero-indexed column number with cluster number. default 8',
                    type=int,default=8)
parser.add_argument('--sep',help='column separator. default <space>',
                    type=str,default=' ')
parser.add_argument('--header',action='store_true',help='treat first row as column headers. default no header',
                    default=False)

args = parser.parse_args()

class arghandle:
    def __init__(self,inbam):
        self.inBAM    = inbam
        self.drop     = False
        self.splitBam = False
        self.prog     ='foo'
        self.region   = None

table = pd.read_csv(args.readInfo,sep=args.sep,header=0 if args.header else None)

print('making clustermap')
clustermap = dict(table[[args.names,args.cluster]].values) 

print('tagging bam')
addHPtag(arghandle(args.inbam),args.outbam,clustermap)
print('Done')
