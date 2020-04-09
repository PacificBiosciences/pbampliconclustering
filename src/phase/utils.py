import pysam
import pandas as pd
import mappy as mp
from itertools import chain
from statistics import median

MINLEN=50
MAXLEN=50000

def getFileType(fname):
    ext = fname.rsplit('.',1)[-1]
    if ext == 'bam':
        return 'bam'
    elif ext in ['fastq','fq','fasta','fa']:
        return 'fastq'
    else:
        raise PhaseUtils_Error(f'unknown filetype extension: {ext}')

def hpCollapse(maxLen=1):
    def csgen(sequence):
        last = None
        for char in sequence:
            if char == last:
                n += 1
            else:
                last = char
                n    = 0
            if n < maxLen:
                yield char
    return lambda seq: ''.join(csgen(seq))

class SimpleRecord:
    def __init__(self,name,sequence):
        self.name     = name
        self.sequence = sequence
    def __len__(self):
        return len(self.sequence)

class RecordGenerator:
    def __init__(self,inFile,fileType=None,region=None,minLength=MINLEN,maxLength=MAXLEN):
        self.inFile    = inFile
        self.region    = region
        self.minLen    = minLength
        self.maxLen    = maxLength
        
        ftype          = getFileType(inFile) if fileType is None else fileType
        self.generator = {'bam'  : self._bamIter,
                          'fastx': self._fastxIter}[ftype]

    def _bamIter(self,bamfile,**kwargs):
        bam = pysam.AlignmentFile(bamfile,check_sq=False)
        if kwargs['region']:
            recgen = bam.fetch(region=kwargs['region'])
        else:
            recgen = bam
        for rec in recgen:
            if rec.query_length >= self.minLen and rec.query_length <= self.maxLen:
                yield SimpleRecord(rec.query_name,rec.query_sequence)

    def _fastxIter(self,fastx,**kwargs):
        for rec in pysam.FastxFile(fastx):
            if len(rec.sequence) >= self.minLen and len(rec.sequence) <= self.maxLen:
                yield SimpleRecord(rec.name,rec.sequence)

    def __iter__(self):
        return self.generator(self.inFile,region=self.region)

class PilerUpper:
    def __init__(self,inFile,region=None,refSeq=None,method='median',
                 minLength=50,maxLength=1e6,maxHP=1):
        self.recGen    = RecordGenerator(inFile,region=region,
                                         minLength=minLength,
                                         maxLength=maxLength)
        self.collapse  = hpCollapse(maxHP) if maxHP else (lambda x:x)
        self.refseq    = self._getRef(refSeq,method)
        self.aligner   = self._getAligner()
        self.varDf     = self._fillVDF()
        
    def _getRef(self,reference,method):
        if reference:
            #refSeq must be a string DNA sequence [ATGC]
            return self.collapse(reference)
        elif method == 'first':
            for rec in self.recGen:
                size = len(rec.sequence)
                if size >= self.minLength and size <= self.maxLength:
                    return self.collapse(rec.sequence)
        elif method == 'median':
            seqs   = pd.Series(rec.sequence for rec in self.recGen)
            medIdx = seqs.str.len().sort_values().index[int(len(seqs)/2)]
            return self.collapse(seqs[medIdx])
        #TODO add random selection method
        else:
            raise ValueError(f'No method named {method}')
            
    def _getAligner(self):
        return mp.Aligner(seq=self.refseq,preset='splice',best_n=1)
    
    def map(self,seq):
        assert self.aligner is not None
        return list(self.aligner.map(self.collapse(seq),cs=True))[0]
    
    def _fillVDF(self):
        print('Reading Variants')
        return pd.concat([self.expandCS(self.map(rec.sequence),
                                        name=rec.name)
                          for rec in self.recGen],
                         axis=1).T

    def getOps(self,csString):
        ops = ':*-+~' 
        op  = None
        val = ''
        for s in csString:
            if s in ops:
                if op:
                    yield (op,val)
                    op,val = s,''
                op = s
            else:
                val += s
        yield (op,val)

    def parseOp(self,opGrp):
        op,val = opGrp
        if op == ':':
            return ('.',)*int(val)
        elif op == '-':
            return (op+val,) + ('*',)*(len(val)-1)
        else:
            return (op+val,)

    def expandCS(self,aln,name=None):
        out = pd.Series(chain(*map(self.parseOp,self.getOps(aln.cs))))
        #fix insertion locations
        insIdx        = out[out.str.contains('\+')].index
        out[insIdx-1] = out[insIdx]
        out.drop(insIdx,inplace=True)
        out.index     = range(aln.r_st,len(out)+aln.r_st)
        if name:
            out.name  = name
        return out

class PhaseUtils_Error(Exception):
    pass

