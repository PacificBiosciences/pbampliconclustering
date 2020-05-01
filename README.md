# Tools for clustering and phasing of PacBio CCS reads

This repo contains python3 tools to cluster PB CCS reads using kmer counts and clustering algorithms provided by the [Python Scikit-learn machine learning toolset](https://scikit-learn.org/stable/index.html).  The primary use case is for amplicon data, where reads cover a specific region in a reference dataset. For non-targeted data, options are provided to cluster any mapped sequence data covering a defined region in a reference sequence (i.e. WGS data).

An alternate tool `LongAmpliconPhasing.py` recently added provides a variant splitting and phasing model for phasing PacBio HiFi/CCS targeted reads.  

## Dependencies

Python 3 is used to take full advantage of the scikit-learn library.  The following packages are required:
 * [scikit-learn](https://scikit-learn.org/stable/install.html)
 * [numpy](https://numpy.org/)
 * [pandas](https://pandas.pydata.org/)
 * [pysam](https://pysam.readthedocs.io/en/latest/api.html)
 * [mappy](https://pypi.org/project/mappy/)
 * [matplotlib](https://matplotlib.org/)
 * [seaborn](https://seaborn.pydata.org/)



# Long Amplicon Phasing

A program for splitting PacBio HiFi reads using shared variants.  This tool will split groups of reads sharing a single variant, including SNV and indels of >20bp in size.  Large insertions, sometimes called structural variants (SV), within amplicon target regions of >1kb can be separated with this tool.

Two models are provided for splitting reads: `align` and `debruijn`.  The align model identifies variants using alignments to a reference *or* by selecting an exemplar read to which all other reads are aligned.  The debruijn model generates a debruijn graph from kmers and splits along graph edges one node at a time.

For general help, see `LongAmpliconPhasing.py -h`.

Example:
    [coming soon]

## Quickstart

The most important parameters for all models are the settings for minimum output group sizes, as well as initial data reduction.

The minimum cluster size is determined by the parameters `-r,--minReads`, `-f, --minFrac`, and the input read count after filtering.  Minimum size is defined as `max( ceil( minFrac * nReads ), minReads )`.

For the `align` method, the parameter `-g,--minSignal` determines how reference positions are filtered prior to splitting.  Only positions for which at least `minSignal` fraction of reads are different from the reference will be considered as candidates for splitting groups. In general, `-g` should be set <= `-f`.  

Note that any read that does not cover _all candidate positions_ will be filtered from the phasing process.

Align example using defaults

    $ python3 LongAmpliconPhasing.py -m align -p outdir/example --reference myref.fasta aligned.bam mySampleName

Debruijn example using defaults

    $ python3 LongAmpliconPhasing.py -m debruijn -p outdir/example input[.bam|.fastq] mySampleName 

## Sequence outputs

### BAM
BAM outputs (with BAM input) can be single or one per cluster (`--splitBam`). Use the `-d` option to drop any read which is not assigned a cluster number (reads outside of the region, if any, secondary/supplementary alignments, partial coverage reads, etc). To turn off all bam export, use the `-x` flag.

### FASTQ
Use the `-F` option to exort one fastq file per cluster.

### clusters.txt
Listing of read assignments by cluster in faux-fasta format.  

    >cluster0_numreads111
    <readname>
    <readname>
    ...

### clusterSplits.txt
Simple graph output defining algorithm decision for phasing variants. 

## Variant Reports

When using a reference with either `align` or `debruijn` method and an aligned bam as input, use the `-v` option to return a set of reports and draft consensus sequences for clustered reads.  When using the debruijn method and/or compressing homopolymers, the input aligned bam will be used for generating variant calls per cluster after phasing.

### alleleClusterSummary.csv
This table summarizes the read clusters identified by the program.  Cluster `0` will _always_ represent reads that match the reference sequence at all candidate positions (see `minSignal`).  Phased reads with 1 or more variants with respect to the reference begin with cluster number `1` and are sorted in descending number of reads.   

Variant calls in each column of this table represent a simple plurality given the reads in the cluster.  

    $ column -ts, output/example.alleleClusterSummary.csv
    contig   HTT_region_hg19  HTT_region_hg19         HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  nReads  frequency
    pos      915              2723                    2774             2779             2782             2785             2792             3495
    cluster
    1        .                .-6CAGCAG               .                .                .                .                .                .                55      0.4135
    2        G                .-18CAGCAGCAGCAGCAGCAG  A                C                C                C                G                G                42      0.3157
    3        .                .-3CAG                  .                .                .                .                .                .                21      0.1578
    4        G                .                       *                C                .                .                G                G                15      0.1127

### sampleVariantSummary.csv
This table is a count of all unique combinations of variants in the dataset, given the candidate positions as described above for `minSignal`.  Reads matching the reference at all positions are labeled cluster `0`, and all other combinations are sorted by number of reads.  There is a hard cut-off of `3` reads such that all unique combinations whith fewer than 3 reads are labeled as a single noise group `-1`.  Actual combinations assigned `-1` can be viewd in the log file.  

    $ column -ts, output/example.sampleVariantSummary.csv
    contig   HTT_region_hg19  HTT_region_hg19         HTT_region_hg19     HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  HTT_region_hg19  nReads  frequency
    pos      915              2723                    2758                2761             2763             2773             2774             2779             2782             2785             2792             3495
    cluster
    0        .                .                       .                   .                .                .                .                .                .                .                .                .                5       0.0375
    1        .                .-6CAGCAG               .                   .                .                .                .                .                .                .                .                .                47      0.3533
    2        G                .-18CAGCAGCAGCAGCAGCAG  .                   .                .                .                A                C                C                C                G                G                40      0.3007
    3        .                .-3CAG                  .                   .                .                .                .                .                .                .                .                .                8       0.0601
    4        .                .-6CAGCAG               .                   .                .                .                .                .                .                .                .                .+1G             5       0.0375
    5        .                .-3CAG                  .                   .                .                .                .                .                .                .                .                .+1G             4       0.0300
    6        G                .                       .-14GCAGCAGCAGCAGC  *                *                .-1G             *                C                C                C                G                G                4       0.0300
    7        G                .                       .                   C                .-11AGCAGCAGCAG  *                *                C                *                *                G                G                3       0.0225
    -1       G                .                       .                   .                .                .                .                .                .                .                .+9CCGCCGCCG     G                17      0.1278

### variantFraction.csv
(Draft) Total variant fraction for each position, regardless of cluster.

### entropy.csv
This table shows the entropy score for each position for each cluster, where entropy is calculated as shannon diversity of calls at a position within a cluster.  High values in this table indicate position/clusters which may be incompletely separated, or include noisy calling regions.  This table can be visualized as a heatplot with the option `-e`.

### laphase.log
Log of algorithm execution.  See for details of splitting.

# Cluster Amplicons

## Quickstart

Examples:
 * [No-Amp Repeat Expansion walkthrough](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/no_amp/README.md)
 * [HTT Repeat Region from WGS HiFi reads](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/no_amp/README.md#clustering-wgs-data-by-region)
 * [HLA Alleles walkthrough](https://github.com/PacificBiosciences/pbampliconclustering/tree/master/examples/hla/README.md)

The clustering tool has two sub-tools.  The first, `describe`, is used for describing the available clustering algorithms and the mapping between command-line options and tool options.  

The second tool, `cluster`, is the primary clustering tool for grouping and labeling CCS reads.

    $ py3 ClusterAmplicons.py -h
    usage: ClusterAmplicons.py [-h] {cluster,describe} ...
    
    Clustering by kmer counts
    
    optional arguments:
      -h, --help          show this help message and exit
    
    subcommands:
      {cluster,describe}
        cluster           cluster reads
        describe          describe models

### Describe Model Inputs

Describe defaults and CL => KW argument map.  Us the -M option for a specific tool, or no arguments to see rules for all clustering algorithms.  Details of what each algorithm accepts can be found on the [scikit-learn web site](https://scikit-learn.org/stable/modules/clustering.html#clustering). 

    $ py3 ClusterAmplicons.py describe -M dbscan
    -----------------DBSCAN-----------------
    ArgMap:
                     eps => eps
                minReads => min_samples
                   njobs => n_jobs
    Defaults:
                     eps => 0.01
             min_samples => 3
                  metric => euclidean
                  n_jobs => 2

The full set of options for any clustering algorithm can be accessed using a `.json` configuration file passed to the option `-P` (see below).

### Cluster Reads
Options and examples discussed below.

    $ py3 ClusterAmplicons.py cluster -h
    usage: ClusterAmplicons.py cluster [-h] [-j,--njobs NJOBS] [-k KMER]
                                       [-z MINIMIZER] [-H] [-T TRIM]
                                       [-M {dbscan,optics,aggcluster,affprop,meanshift,kmeans}]
                                       [-a {pca,featagg}] [-c COMPONENTS] [-e EPS]
                                       [-m MINREADS] [-n {l1,l2,none}]
                                       [-i IGNOREENDS] [-P PARAMS] [-r REGION]
                                       [--extractReference REFERENCE] [-q MINQV]
                                       [-l MINLENGTH] [-L MAXLENGTH]
                                       [-w WHITELIST] [-N NREADS] [-f FLANKS] [-A]
                                       [-s SEED] [-p PREFIX] [-S] [-x] [-F] [-d]
                                       [-t] [-g PLOTREADS] [-X]
                                       [inBAM]
    
    positional arguments:
      inBAM                 input BAM of CCS alignments. Default stdin
    
    optional arguments:
      -h, --help            show this help message and exit
      -j,--njobs NJOBS      j parallel jobs (only for some models). Default 1
    
    kmers:
      -k KMER, --kmer KMER  kmer size for clustering. Default 11
      -z MINIMIZER, --minimizer MINIMIZER
                            group kmers by minimizer of length z. Default 0 (no
                            minimizer)
      -H, --noHPcollapse    do not compress homopolymers. Default collapse HP
      -T TRIM, --trim TRIM  Trim kmers with frequency < trim. Default 0.10
    
    cluster:
      -M {dbscan,optics,aggcluster,affprop,meanshift,kmeans}, --model {dbscan,optics,aggcluster,affprop,meanshift,kmeans}
                            clustering model. See https://scikit-
                            learn.org/stable/modules/clustering.html. Default
                            dbscan
      -a {pca,featagg}, --agg {pca,featagg}
                            Feature reduction method. Default pca
      -c COMPONENTS, --components COMPONENTS
                            Use first c components of PCA/FeatAgg for clustering.
                            Set to 0 for no reduction. Default 2
      -e EPS, --eps EPS     eps cluster tolerance. Default None
      -m MINREADS, --minReads MINREADS
                            Minimum reads to be a cluster. Default 5
      -n {l1,l2,none}, --normalize {l1,l2,none}
                            normalization of kmer counts. Default l2
      -i IGNOREENDS, --ignoreEnds IGNOREENDS
                            ignore i bases at ends of amplicons for clustering.
                            Default 0
      -P PARAMS, --params PARAMS
                            json file of parameters for specific model. Order of
                            precedence: json > CL-opts > defaults. Default None
    
    filter:
      -r REGION, --region REGION
                            Target region for selection of reads, format
                            '[chr]:[start]-[stop]'. Example '4:3076604-3076660'.
                            Default all reads (no region)
      --extractReference REFERENCE
                            Extract subsequence at region coordinates for
                            clustering using fasta reference (must have .fai).
                            Maps 100nt on either side of region to each read and
                            extracts sequence inbetween for kmer counting. Default
                            None (use full read)
      -q MINQV, --minQV MINQV
                            Minimum quality [0-1] to use for clustering. Default
                            0.99
      -l MINLENGTH, --minLength MINLENGTH
                            Minimum length read to use for clustering. Default 500
      -L MAXLENGTH, --maxLength MAXLENGTH
                            Maximum length read to use for clustering. Default
                            25000
      -w WHITELIST, --whitelist WHITELIST
                            whitelist of read names to cluster. Default None
      -N NREADS, --nReads NREADS
                            Randomly downsample to nReads after filtering. Default
                            0 (all avail reads)
      -f FLANKS, --flanks FLANKS
                            fasta of flanking/primer sequence. Reads not mapping
                            to both will be filtered. Default None
      -A, --noArtifactFilter
                            Turn off palindromic-artifact filtering. Default use
                            artifact filter
      -s SEED, --seed SEED  Random seed for downsampling. Default 17
    
    output:
      -p PREFIX, --prefix PREFIX
                            Output prefix. Default ./clustered
      -S, --splitBam        split clusters into separate bams (noise and no-
                            cluster dropped). Default one bam
      -x, --noBam           Do not export HP-tagged bam of clustered reads
      -F, --fastq           Export one fastq per cluster
      -d, --drop            Drop reads with no cluster in output bam. Default keep
                            all reads.
      -t, --testPlot        Plot reads vs dist to nearest m-neighbors without
                            clustering
      -g PLOTREADS, --plotReads PLOTREADS
                            Write pairplot of first g reduced axes for each read.
                            Default None (no plot)
      -X, --exportKmerTable
                            Export kmer count table after trimming. Default False

## Region Selection
Clustering can occur for all reads, a subset of reads, or over a defined reference window spanned by a subset of reads.  By default, all sequence in the input bam will be characterized by kmer counts and clustered.  

If a *region* is provided without an _extractReference_, then all reads intersecting the region (returned by pysam _fetch_ method) will be characterized and clustered.

If a *region* and *extractReference* are both provided, then only the sequence between the reference coordinates is clustered from reads completely spanning the region.  Sequence between region coordinates is extracted by mapping 100bp of flanking sequence from the reference to each mapped read returned by pysam _fetch_.   

### Filtering
Reads are filtered by minimum read quality `-q` [0-1], default `0.99`.  For extracted sequence, the QV filter is applied to the extracted sequence only. 

Primer sequences can be supplied to filter artifacts.  Reads will only be included in clustering analysis if both primers occur in the read.

Potential sequencing artifacts with missing adapters are automatically removed.  To turn off this filter, use the `-A` flag.

## Clustering
Clustering is based on kmer count vectors for each read in the input dataset, following region selection and filtering.  

### Kmers
By default homopolymer stretches (n>=2) are compressed prior to kmer counting.  This step reduces noise caused by one of the primary sources of error in PB sequencing.  This option can be turned off with the `-H` option.  

Kmers can be grouped by a _minimizer_ of size `-z`.  This is a naive implementation that labels all kmers by the first lexicographically sorted substring of length _z_.  

Kmers of frequency less than `T` or greater than `1 - T` in the dataset will be removed prior to clustering.

### Feature Reduction
[PCA](https://scikit-learn.org/stable/modules/decomposition.html#principal-component-analysis-pca) or [feature agglomeration](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.FeatureAgglomeration.html#sklearn.cluster.FeatureAgglomeration) can be used to reduce the number of clustering features.  The option `-a,--agg` sets the method, and `-c` determines the number of used components (PCA) or output features (featagg).  Setting the number of components to 0 turns off feature reduction.

### Normalize
Kmer counts are [normalized](https://scikit-learn.org/stable/modules/preprocessing.html#normalization) _within samples_ unless `-n` is set to `none`.

### Ignore Ends
To avoid clustering reads based on degenerate primers, this option can be set to ignore sequence `-i` bases from the ends of each read.

### Minimum Cluster size
Clusters must have at least `-m` reads.  Clusters with less than `-m` reads will be reclassified as _noise_.  

### Custom Parameters
A simple json file can be provided to set all options for any clustering algorithm.  The json config file trumps all other input parameters (ie defaults and CL options).  See [example json file](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/optics_config.json) for the OPTICS algorithm.

## Outputs
### clusters.txt
The primary output is a text file listing reads in each output cluster.  Reads have their original name, unless the `--extracReference` option is used to extract a subsequence from each read, in which case the extraction coordinates will be appended to the read names.

    >cluster0_numreads42
    m54309_190912_232752/34538433/ccs
    m54309_190912_232752/64291805/ccs
    m54309_190912_232752/70058377/ccs
    ...
    >cluster1_numreads31
    m54309_190912_232752/8847473/ccs
    m54309_190912_232752/40436366/ccs
    m54309_190912_232752/41288675/ccs
    ...
    >Noise_numreads2
    m54309_190912_232752/45744303/ccs
    m54309_190912_232752/47055558/ccs

Reads filtered prior to clustering are *not* listed.

### HP-tagged BAM
Cluster numbers are inserted into each row of the output BAM using the `HP` tag.  If the `-d` option is passed, only clustered reads will be included in the output.  Otherwise, filtered reads are labeled `999` and reads that enter the clustering process but are classified as _noise_ are labeled `-1`.  All output reads also have an RGB color defined by cluster in the `YC` tag for visualization in IGV.  The option `-S` generates a single BAM output per cluster, and `-x` will prevent any bam output from being written.

![Cluster alleles](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/igv_3Alleles.png)

### Fastq per cluster
Use the `-F` option to export a fastq file per cluster.  This can be used as input for [consensus](https://github.com/armintoepfer/c3s).

### Nearest Neighbor plot
For some clustering algorithms (e.g. DBSCAN), it can be useful to view a plot of sorted nearest neightbor distances to set the _eps_ value.  The option `-t` generates such a plot for a given parameter set and read input.

![EPS Estimator](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/no_amp/allTargets50.eps_estimator.png)

### Cluster Plot
The option `-g` generates a plot of each read position given the first 2 reduced components from the input matrix.

![DRB split](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/clusterDRB.clusters.png)



THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
