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
       packages == Cython, Pandas, Numpy, Scipy, pysam, bx-python, setuptools, future, regex, matplotlib,
      
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
   