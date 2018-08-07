# germline variant calling

### This pipeline is builting on open softwares! Please ckeck the requirment!


## 1. Requirement
   1.1 Envir
   1.1.1 python
       python == 3.6
       ./configure --prefix=/opt/python3.6
       make
       sudo make install
       #python3.6 -m pip install pkg --user
       packages == Cython, Pandas, Numpy, Scipy, pysam, bx-python, setuptools, future, regex, matplotlib, pyyaml
      
   1.1.2 Java
       java == 1.8

   1.1.3 R
       R == 3.4.4
       packages = ggplot2 , optagrse , reshape2 , splines
       
   1.1.4 Git
       Git == 2.7.4
       git-lfs == 2.4.2
       # download https://github.com/git-lfs/git-lfs/releases/download/v2.4.2/git-lfs-linux-amd64-2.4.2.tar.gz
         git lfs install

   1.1.5 CMAKE
       cmake > 2.8
       
   1.1.6 GCC/G++
       GCC/G++ >= 4.8+
   
   1.1.7 Boost
       Boost >= 1.55
   
   1.1.8 zlib1g-dev
       sudo apt-get zlib1g-dev

   1.2 Bioinformation softwares 
   1.2.1 BWA (alignment via Burrows-Wheeler transformation)
       bwa == 0.7.12-r1039
       
   1.2.2 samtools
       samtools == 1.7
        autoheader     # If using configure, generate the header template...
        autoconf       # ...and configure script (or use autoreconf to do both)
        ./configure    # Optional, needed for choosing optional functionality
        make
        make install

   1.2.3 GATK
       GATK == 4.0.5
       
   1.2.4 UMI-tools
       git clone https://github.com/CGATOxford/UMI-tools.git
       cd UMI-tools-master
       python3.6 setup.py install --user
        
   1.2.5 hay.py--benchmarking
       git clone https://github.com/sequencing/hap.py
       mkdir hap.py-build
       cd hap.py-build
       cmake ../hap.py
       make


## 2. Datasets
   2.1 reference genome
       Download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
       ucsc.hg19.fasta.gz
      
   2.2 known-sites
       Download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/
       1000G_phase1.snps.high_confidence.hg19.sites.vcf
       Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
       dbsnp_138.hg19.vcf

   2.3 truth vcf
       Download from ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/
       ConfidentRegions.bed.gz
       NA12878/NA12878.vcf.gz
       NA12877/NA12877.vcf.gz
   
   2.4 target bem 
       target_breast.bed
   2.5 primer file
       DHS-001Z_primers_target.csv

##3. USAGE 
   3.1 install 
     # clone from the github
     git clone https://github.com/WXyanengbio/pip_germlineVC.git
   3.2 get the installed path of the bioinformation software and modify the run1.py or input the parameters by command line
     line63 in run1.py :  default = 'fastqc'
     line68 in run1.py :  default = 'bwa'
     line73 in run1.py :  default = 'samtools'
     line78 in run1.py :  default = 'gatk'
     line83 in run1.py :  default = '/home/dell/.local/bin/umi_tools'
     line88 in run1.py :  default = '/home/dell/Works/Softwares/bioapps/benchmarking-tools/tools/hap.py_build/bin/hap.py'
   3.3 put the files into the paths of the datasets---Datasets and modify the run1.py or input the parameters by command line
     line52 in run1.py : default = '/home/dell/Works/Projects/Datasets'
     3.3.1 reference genome
         ucsc.hg19.fasta  --- Datasets/genome/ucsc.hg19.fasta
         ucsc.hg19.fasta.fai  --- Datasets/genome/ucsc.hg19.fasta.fai 
         ucsc.hg19.dict  --- Datasets/genome/ucsc.hg19.dict
         ucsc.hg19.amb  --- Datasets/genome/ucsc.hg19.amb
         ucsc.hg19.ann  --- Datasets/genome/ucsc.hg19.ann
         ucsc.hg19.bwt  --- Datasets/genome/ucsc.hg19.bwt
         ucsc.hg19.pac --- Datasets/genome/ucsc.hg19.pac
         ucsc.hg19.sa --- Datasets/genome/ucsc.hg19.sa
     3.3.2 known-sites
         1000G_phase1.snps.high_confidence.hg19.sites.vcf  --- Datasets/known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf 
         Mills_and_1000G_gold_standard.indels.hg19.sites.vcf  --- Datasets/known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf 
         dbsnp_138.hg19.vcf  --- Datasets/known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf 
     3.3.3 truth vcf 
         ConfidentRegions.bed Datasets/truthDB/confidentregion/ConfidentRegions.bed
         NA12878.vcf.gz ---  Datasets/truthDB/samll_variant/NA12878.vcf.gz
         NA12877.vcf.gz ---  Datasets/truthDB/samll_variant/NA12877.vcf.gz
     3.3.4 target files
         target_breast.bed  --- Datasets/genome/target_breast/target_breast.bed 
         target_breast.refSeq.fa ---Datasets/genome/target_breast/target_breast.refSeq.fa
     3.3.5 primer file
         put the primer file into the path of the raw reads
   3.4 running
     python3.6 pip_germlineVC/run1.py --source  dir of raw reads --sample_name sampleID
   3.4.1 the optional arguments
        -h, --help                         show this help message and exit
       --source                            Path to input reads in FASTA format
       --sample_name                       the sample name of raw reads
       --tailname                          the tailname of sample raw read
       --datasets_dir                      the tailname of sample raw read
      --common_seq1                        the common seq1 of QIAseq Targeted DNA Panel
      --common_seq2                        the common seq2 of QIAseq Targeted DNA Panel
      --output                             Path of output file
      --fastqc_dir                         the install path of fastQC
      --bwa_dir                            the install path of bwa
      --samtools_dir                       the install path of samtools
      --gatk_dir                           the install path of GATK4
      --umitools_dir                       the install path of umitools
      --benchmark_dir                      the install path of benchmark
      --ref_index_name                     the path of ref index--if there isn't a ref index, it will make a index in the path of ref fasta by bwa
      --ref_fa_file                        the path of ref fasta
      --total_ref_fa_file                  the path of ref total ref fa
      --total_ref_fa_dict                  the path of ref total ref fa dict
      --truth_vcf                          the path of truth VCF of variant for benchmarking
      --confident_region_bed               the path of confident region bed file of the truth VCF of variant for benchmarking
      --min_read_len                       the cutoff of the min read length
      --num_threads                        the number of threads to align
      --min_mapq                           the parameter of filter alignment_sam
      --max_soft_clip                      the parameter of filter alignment_sam
      --max_dist                           the parameter of filter alignment_sam
      --primers_file                       Load all primer sequences in the panel
      --edit_dist                          the parameter of edit distance between barcodes
      --memory_size                        the cutoff of Java memory
      --known_sites                        the list of --known-sites , sep by: ,
      --exome_target_bed                   the bed file of exome intervals
      --erc                                switch to running HaplotypeCaller in GVCF mode
      --read_filter                        add a read filter that deals with some problems
      --snp_filter                         add parameters for filtering SNPs
      --indel_filter                       add parameters for filtering Indels
      --db_cosmic                          add cosmic databases of variants
      --db_clinvar                         add clinvar databases of variants
      --db_g1000                           add g1000 databases of variants
      --anno_geneID                        add annotation gene ID of variants
      -v,                                  show program's version number and exit
      --test                               the subprocess of the script
