import pysam
import pandas as pd
import mappy as mp
from itertools import chain
from statistics import median
from collections import Counter
from scipy.stats import entropy

MINLEN=50
MAXLEN=50000

def summary(splitter,caller):
    try:
        minsig   = splitter.minSignal
        minfrac  = splitter._minFrac
        minreads = splitter._minReads
        feat     = len(splitter.sigVar)
    except AttributeError:
        minsig,minfrac,minreads = ('NA',)*3
        feat = len(splitter.nodes)
    mincnt = splitter.minCount
    total  = splitter.stats['total alignments']
    prim   = splitter.stats['primary alignments']
    plup   = splitter.stats['pileup alignments']
    clust  = splitter.stats['clustered reads']
    ref    = caller.readCounts[0]
    lrg    = max(caller.readCounts.values())
    noise  = caller.readCounts[-1]
    conly  = {k:v for k,v in caller.readCounts.items() if k!=-1}
    nclust = len(conly)
    entr   = entropy(list(conly.values()))
    tcount = Counter(caller.variantGroupMap.values())
    tref   = tcount[0]
    tlrg   = max(tcount.values())
    tnoise = tcount[-1]
    tconly = {k:v for k,v in tcount.items() if k!=-1}
    ngrps  = len(tconly)
    tentr  = entropy(list(tconly.values()))
    return f'''
Input Params
------------
minSignal:\t\t{minsig}
minFrac:\t\t{minfrac}
minReads:\t\t{minreads}

Read Counts
-----------
Total Alignments:\t{total:,}
Primary Alignments:\t{prim:,}\t({prim/total:.2})
Pileup Reads:\t\t{plup:,}\t({plup/total:.2})
Clustered Reads:\t{clust:,}\t({clust/total:.2}) (covering all variant pos)

Cluster Info
------------
Min Cluster Size:\t{mincnt:,}
Cluster Features:\t{feat:,}

Cluster Stats (out of {clust:,})
-------------
Reference Calls:\t{ref:,}\t({ref/clust:.2})
Largest Fraction:\t{lrg:,}\t({lrg/clust:.2})
Noise Reads:\t\t{noise:,}\t({noise/clust:.2})
N Clusters:\t\t{nclust:,}
Shannon Entropy:\t{entr:.4}

Total Fractions (out of {clust:,})
---------------
Reference Calls:\t{tref:,}\t({tref/clust:.2})
Largest Fraction:\t{tlrg:,}\t({tlrg/clust:.2})
Noise Reads(<3):\t{tnoise:,}\t({tnoise/clust:.2})
Unique Var Comb:\t{ngrps:,}
Shannon Entropy:\t{tentr:.4}
    '''

def getFileType(fname):
    ext = fname.rsplit('.',1)[-1]
    if ext == 'bam':
        return 'bam'
    elif ext in ['fastq','fq','fasta','fa']:
        return 'fastx'
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

def writeSimpleBED(chrm,start,stop,name,cov,filename,mode='w'):
    with open(filename,mode) as ofile:
        ofile.write('\t'.join(map(str,[chrm,start,stop,name,cov])) + '\n')

def writeRegionBam(inBam,outBam,region):
    with pysam.AlignmentFile(inBam) as ibam:
        with pysam.AlignmentFile(outBam,'wb',template=ibam) as obam:
            for rec in ibam.fetch(region=region):
                obam.write(rec)
    pysam.index(outBam)
    return outBam

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
        self.counter   = Counter()

    def getNameIdx(self):
        '''Run through without returning sequence'''
        return {rec.name:i 
                for i,rec in enumerate(self.generator(self.inFile,
                                                      region=self.region,
                                                      track=False))}

    def report(self):
        other = ",".join([f"{n}:{c}" for n,c in self.counter.items() if n!="pass"])
        return f'Alignments loaded: {self.counter["pass"]}; filtered: {other}'

    def _bamIter(self,bamfile,track=True,**kwargs):
        bam = pysam.AlignmentFile(bamfile,check_sq=False)
        if kwargs['region']:
            recgen = bam.fetch(region=kwargs['region'])
        else:
            recgen = bam
        for rec in recgen:
            kind = self._classifyBam(rec)
            if track:
                self.counter[kind] += 1
            if kind == 'pass':
                yield SimpleRecord(rec.query_name,rec.query_sequence)

    def _fastxIter(self,fastx,track=True,**kwargs):
        for rec in pysam.FastxFile(fastx):
            kind = self._classifyFq(rec)
            if track:
                self.counter[kind] += 1
            if kind == 'pass':
                yield SimpleRecord(rec.name,rec.sequence)

    def _classifyBam(self,rec):
        if rec.flag & 0x100:
            return 'secondary'
        if rec.flag & 0x800:
            return 'supplementary'
        if rec.query_length < self.minLen:
            return f'short(<{self.minLength})'
        if rec.query_length > self.maxLen:
            return f'long(>{self.maxLength})'
        else:
            return 'pass'

    def _classifyFq(self,rec):
        if len(rec.sequence) < self.minLen:
            return f'short(<{self.minLength})'
        if len(rec.sequence) > self.maxLen:
            return f'long(>{self.maxLength})'
        else:
            return 'pass'

    def __iter__(self):
        return self.generator(self.inFile,region=self.region)

class RecTracker:
    def __init__(self,):
        self.counter = Counter()
    def __call__(self,rec):
        try:
            if rec.flag & 0x100:
                self.counter['secondary'] +=1
            elif rec.flag & 0x800:
                self.counter['supplementary'] +=1
            else:
                self.counter['primary'] +=1
        except AttributeError: #no flag == fastarec
            self.counter['primary'] +=1
    def __repr__(self):
        other = ",".join([f"{n}:{c}" for n,c in self.counter.items() if n!="primary"])
        return f'Alignments loaded: {self.counter["primary"]}; filtered: {other}'
        

class PilerUpper:
    def __init__(self,inFile,region=None,refSeq=None,method='median',
                 minLength=50,maxLength=1e6,maxHP=1,log=None,
                 multifunc=None,nproc=1):
        self.recGen    = RecordGenerator(inFile,region=region,
                                         minLength=minLength,
                                         maxLength=maxLength)
        self.collapse  = hpCollapse(maxHP) if maxHP else (lambda x:x)
        self.log       = log
        self.nproc     = nproc
        self.refseq    = self._getRef(refSeq,method)
        self.aligner   = self._getAligner()
        self.getRow    = self._getRow if nproc==1 else multifunc
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
        try:
            return list(self.aligner.map(self.collapse(seq),cs=True))[0]
        except IndexError:
            raise PhaseUtils_Error(f'Unable to align sequence: {seq}')

    def _getRow(self,rec):
        return self.expandCS(self.map(rec.sequence),name=rec.name)

    def _fillVDF(self):
        if self.log:
            self.log.info('Aligning compressed reads')
        if self.nproc != 1:
            self.log.warn('Parallel processing of hp compression not implemented yet. Proceeding with 1 proc')
        result = map(self._getRow,self.recGen)
        if self.log:
            self.log.debug(self.recGen.report())  
        return pd.concat(result,axis=1).T

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

