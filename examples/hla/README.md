# Clustering CCS reads from Pooled HLA Alleles

This example starts from a single pooled HLA amplicon dataset containing up to 9 loci: 
    
    A,B,C,DPB1,DQB1,DRB1,DRB3,DRB4,DRB5

The data come from a multiplexed SMRTcell with barcodes.  We will look at just one sample, `lbc90--lbc90`.

## Separate HLA loci by alignment

We use a de-duplicated subset of HLA alleles for separating CCS reads by locus.  This is the same guide fasta used for PacBio [LAA with Guided Clustering](https://github.com/PacificBiosciences/pblaa).  CCS reads from barcode lbc90 are aligned to the guide in `mapped.lbc90--lbc90.consensusalignmentset.bam`.  We can identify the HLA-A alleles by whitelisting CCS reads mapped to HLA-A in the guide and clustering.

    $ samtools view -F 0x900 mapped.lbc90--lbc90.consensusalignmentset.bam | awk '/HLA-A/ {print $1}' > A.whitelist

## Cluster whitelisted CCS Reads

The list of names is used to filter the input dataset and cluster.  

    $ py3 ClusterAmplicons.py cluster \
                              -M meanshift \
                              -S \
                              -g \
                              -q 0.999 \
                              -c 4 \
                              -k 15 \
                              -m 5 \
                              -w A.whitelist \
                              -p clusterA.split \
                              mapped.lbc90--lbc90.consensusalignmentset.bam
    Reading Sequence
    Trimming low-freq kmers
    Normalizing data
    Reducing Features with pca
    Clustering 189 reads with meanshift
            bandwidth=None
            bin_seeding=True
            min_bin_freq=5
            cluster_all=False
            n_jobs=2
    Writing 2 clusters with nreads 81,67
    41 reads identified as noise
    Adding HP tag to bam 

The two alleles are easily separated, and noisy reads are removed.
![HLA-A Clusters](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/clusterA.clusters.png)

One of the alleles, showing filtered reads in white (removed before clustering) and noise reads in gray (passed filter but not part of a cluster).  Filtered and Noise reads are automatically removed when the option `--splitBam` is used or if the `--drop` flag is present.
![HLA-A IGV](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/A_02-06-01.cluster.png) 

## Quick Consensus with [bcftools](http://samtools.github.io/bcftools/bcftools.html)
We can generate a quick consensus for cluster `1` using the first clustered read as a 'reference':

    $ samtools fasta clusterA.split.hptagged.1.bam | awk '/^>/ {n++} n>1 {exit} {print}' > clust1_1.fasta && samtools faidx clust1_1.fasta
    
    $ pbmm2 align --sort --preset CCS clust1_1.fasta clusterA.split.hptagged.1.bam \
        | bcftools mpileup -Q 0 -B -Ou -f clust1_1.fasta - \
        | bcftools call --ploidy 1 -mv -Ob \
        | bcftools view -i 'QUAL>30' -Ob > clust1_polish.bcf && bcftools index clust1_polish.bcf
    
    $ bcftools consensus -f clust1_1.fasta clust1_polish.bcf > clust1_cons.fasta

For this sample, the known allele is HLA-A\*02:06:01.  Aligning our consensus with the expected allele shows zero mismatches:

    $ minimap2 --cs HLA00011_A_02-06-01_3517_bp.fasta clust1_cons.fasta 2>/dev/null
    m54043_190914_194303/22282443/ccs       3227    0       3227    +       HLA:HLA00011_A*02:06:01_3517_bp 3517    89      3316    3227    3227    60      NM:i:0  ms:i:6454       AS:i:6454      nn:i:0   tp:A:P  cm:i:598        s1:i:3218       s2:i:0  de:f:0  cs:Z::3227  


## Clustering >1 Locus
For a second example, we can easily cluster all the DRB alleles in a single pass.  

    $ samtools view -F 0x900 mapped.lbc90--lbc90.consensusalignmentset.bam | awk '/HLA-DRB/ {print $1}' > DRB.whitelist
    
    $ py3 ClusterAmplicons.py cluster \
                              -M meanshift \
                              -S \
                              -g \
                              -q 0.999 \
                              -c 6 \
                              -k 15 \
                              -m 5 \
                              -w DRB.whitelist \
                              -p clusterDRB.split \
                              -F \
                              mapped.lbc90--lbc90.consensusalignmentset.bam
    Reading Sequence
    Trimming low-freq kmers
    Normalizing data
    Reducing Features with pca
    Clustering 516 reads with meanshift
            bandwidth=None
            bin_seeding=True
            min_bin_freq=5
            cluster_all=False
            n_jobs=2
    Writing 4 clusters with nreads 233,158,76,21
    28 reads identified as noise
    Adding HP tag to bam
    Exporting fastq

![DRB split](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/clusterDRB.clusters.png)

## Quick Consensus with [c3s](https://github.com/armintoepfer/c3s)
We can also make consensus calls using `c3s` which uses a POA weighted by the QV scores to generate a _denovo_ consensus.

    $ parallel 'c3s --name Barcode0_ClusterDRB_Phase{}_NumReads$(grep -c ^@ clusterDRB.split.cluster{}.fastq) clusterDRB.split.cluster{}.fastq DRB.cluster{}.consensus.fasta' ::: 0 1 2 3

In this case, the IMGT database does not contain the full-length genomic DNA for these alleles, but we can match perfect exons for all expected alleles:

    $ cat DRB.typing                                                      
    Sequence                               gLen  gType                gPctId nMis Indel cLen cType             cPctId Type
    Barcode0_ClusterDRB_Phase3_NumReads21  4982  HLA-DRB1*04:03:01    99.88  2    4     552  HLA-DRB1*04:06:01 100.0  HLA-DRB1*04:06:01
    Barcode0_ClusterDRB_Phase1_NumReads159 3836  HLA-DRB1*12:01:01    99.61  3    12    578  HLA-DRB1*12:02:01 100.0  HLA-DRB1*12:02:01
    Barcode0_ClusterDRB_Phase0_NumReads233 3833  HLA-DRB3*01:01:02:01 97.79  66   19    577  HLA-DRB3*03:01:03 100.0  HLA-DRB3*03:01:03
    Barcode0_ClusterDRB_Phase2_NumReads76  4313  HLA-DRB4*01:03:01:01 99.93  0    3     578  HLA-DRB4*01:03:01 100.0  HLA-DRB4*01:03:01



