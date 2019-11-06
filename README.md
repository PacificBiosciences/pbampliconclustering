# Tools for clustering of PacBio CCS reads

This repo contains python3 tools to cluster PB CCS reads using kmer counts and clustering algorithms provided by the [Python Scikit-learn machine learning toolset](https://scikit-learn.org/stable/index.html).  The primary use case is for amplicon data, where reads cover a specific region in a reference dataset. For non-targeted data, options are provided to cluster any mapped sequence data covering a defined region in a reference sequence (i.e. WGS data).

## Dependencies

Python 3 is used to take full advantage of the scikit-learn library.  The following packages are required:
 * [scikit-learn](https://scikit-learn.org/stable/install.html)
 * [numpy](https://numpy.org/)
 * [pandas](https://pandas.pydata.org/)
 * [pysam](https://pysam.readthedocs.io/en/latest/api.html)
 * [mappy](https://pypi.org/project/mappy/)
 * [matplotlib](https://matplotlib.org/)
 * [seaborn](https://seaborn.pydata.org/)

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

The full set of options for any clustering algorithm can be accessed using `.json` a configuration file passed to the option `-P` (see below).

### Cluster Reads
Options and examples discussed below.

    $ py3 ClusterAmplicons.py cluster -h
    usage: ClusterAmplicons.py cluster [-h] [-j,--njobs NJOBS] [-k,--kmer KMER]
                                       [-z,--minimizer MINIMIZER]
                                       [-H,--noHPcollapse] [-T TRIM]
                                       [-M,--model {dbscan,optics,aggcluster,affprop,meanshift,kmeans}]
                                       [-a,--agg {pca,featagg}]
                                       [-c,--components COMPONENTS] [-e,--eps EPS]
                                       [-m,--minReads MINREADS]
                                       [-n,--normalize {l1,l2,none}]
                                       [-i,--ignoreEnds IGNOREENDS]
                                       [-P,--params PARAMS] [-r,--region REGION]
                                       [--extractReference REFERENCE]
                                       [-q,--minQV MINQV]
                                       [-w,--whitelist WHITELIST]
                                       [-f,--flanks FLANKS] [-p,--prefix PREFIX]
                                       [-S,--splitBam] [-x,--noBam] [-d,--drop]
                                       [-t,--testPlot] [-g,--plotReads]
                                       [inBAM]
    
    positional arguments:
      inBAM                 input BAM of CCS alignments. Default stdin
    
    optional arguments:
      -h, --help            show this help message and exit
      -j,--njobs NJOBS      j parallel jobs (only for some models). Default 1
    
    kmers:
      -k,--kmer KMER        kmer size for clustering. Default 11
      -z,--minimizer MINIMIZER
                            group kmers by minimizer of length z. Default 0 (no
                            minimizer)
      -H,--noHPcollapse     do not compress homopolymers. Default collapse HP
      -T TRIM, --trim TRIM  Trim kmers with frequency < trim. Default 0.10
    
    cluster:
      -M,--model {dbscan,optics,aggcluster,affprop,meanshift,kmeans}
                            clustering model. See https://scikit-
                            learn.org/stable/modules/clustering.html. Default
                            dbscan
      -a,--agg {pca,featagg}
                            Feature reduction method. Default pca
      -c,--components COMPONENTS
                            Use first c components of PCA/FeatAgg for clustering.
                            Set to 0 for no reduction. Default 2
      -e,--eps EPS          eps cluster tolerance. Default None
      -m,--minReads MINREADS
                            Minimum reads to be a cluster. Default 5
      -n,--normalize {l1,l2,none}
                            normalization of kmer counts. Default l2
      -i,--ignoreEnds IGNOREENDS
                            ignore i bases at ends of amplicons for clustering.
                            Default 0
      -P,--params PARAMS    json file of parameters for specific model. Order of
                            precedence: json > CL-opts > defaults. Default None
    
    filter:
      -r,--region REGION    Target region for selection of reads, format
                            '[chr]:[start]-[stop]'. Example '4:3076604-3076660'.
                            Default all reads (no region)
      --extractReference REFERENCE
                            Extract subsequence at region coordinates for
                            clustering using fasta reference (must have .fai).
                            Maps 100nt on either side of region to each read and
                            extracts sequence inbetween for kmer counting. Default
                            None (use full read)
      -q,--minQV MINQV      Minimum quality [0-1] to use for clustering. Default
                            0.99
      -w,--whitelist WHITELIST
                            whitelist of read names to cluster. Default None
      -f,--flanks FLANKS    fasta of flanking/primer sequence. Reads not mapping
                            to both will be filtered. Default None
    
    output:
      -p,--prefix PREFIX    Output prefix. Default ./clustered
      -S,--splitBam         split clusters into separate bams (noise and no-
                            cluster dropped). Default one bam
      -x,--noBam            Do not export HP-tagged bam of clustered reads
      -d,--drop             Drop reads with no cluster in output bam. Default keep
                            all reads.
      -t,--testPlot         Plot reads vs dist to nearest m-neighbors without
                            clustering
      -g,--plotReads        Plot first 2 axes of PCA for each read. Default no
                            plot generated
## Region Selection
Clustering can occur for all reads, a subset of reads, or over a defined reference window spanned by a subset of reads.  By default, all sequence in the input bam will be characterized by kmer counts and clustered.  

If a *region* is provided without an _extractReference_, then all reads intersecting the region (returned by pysam _fetch_ method) will be characterized and clustered.

If a *region* and *extractReference* are both provided, then only the sequence between the reference coordinates is clustered from reads completely spanning the region.  Sequence between region coordinates is extracted by mapping 100bp of flanking sequence from the reference to each mapped read returned by pysam _fetch_.   

### Filtering
Reads are filtered by minimum read quality `-q` [0-1], default `0.99`.  For extracted sequence, the QV filter is applied to the extracted sequence only. 

Primer sequences can be supplied to filter artifacts.  Reads will only be included in clustering analysis if both primers occur in the read.

## Clustering
Clustering is based on kmer count vectors for each read in the input dataset, following region selection and filtering.  

### Kmers
By default homopolymer stretches (n>=2) are compressed prior to kmer counting.  This step reduces noise caused by one of the primary sources of error in PB sequencing.  This option can be turned off with the `-H` option.  

Kmers can be grouped by a _minimizer_ of size `-z`.  This is a naive implementation that labels all kmers by the first lexicographically sorted substring of length _z_.  

Kmers of frequency less than `-T` in the dataset will be removed prior to clustering.

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

### Nearest Neighbor plot
For some clustering algorithms (e.g. DBSCAN), it can be useful to view a plot of sorted nearest neightbor distances to set the _eps_ value.  The option `-t` generates such a plot for a given parameter set and read input.

![DRB split](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/allTargets50.eps_estimator.png)

### Cluster Plot
The option `-g` generates a plot of each read position given the first 2 reduced components from the input matrix.

![DRB split](https://github.com/PacificBiosciences/pbampliconclustering/blob/master/examples/hla/clusterDRB.clusters.png)



THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
