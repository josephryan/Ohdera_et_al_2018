A. Sequenced Reads Trimming and Quality Filtering

A1. FASTQ files generated from IlluminaHiseq2000 sequencer were trimmed to remove adapters and filtered by base quality. 
Dependencies: Trimmomatic v0.36 (Bolger et al. 2014)

java -jar trimmomatic-0.33.jar PE -phred33 reads_forward.fq.gz reads_reverse.fq.gz trim_reads.For.paired.fq.gz trim_reads.Rev.paired.fq.gz trim_reads.For.unp.fq.gz trim_reads.Rev.unp.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDO:4:15 MINLEN:50

B. Error Correction with ErrorCorrectReads.pl from Allpaths-LG

Dependencies: Allpaths-LG version 52488 (Gnerre et al. 2011) 

perl /usr/local/allpathslg-52488/src/ErrorCorrectReads.pl PAIRED_READS_A_IN=trim_reads.For.paired.fq.gz PAIRED_READS_B_IN=trim_reads.Rev.paired.fq.gz UNPAIRED_READS_IN=trim_reads.For.Rev.unp.fq.gz PAIRED_SEP=100 THREADS=30 PHRED_ENCODING=33 READS_OUT=Reads.ErrorCorrect

C. De novo Genome Assembly and gap closing

Dependencies: Platanus version 1.2.1 (Kajitani et al. 2014)

plat.pl —out=Assembly.Plat.89 —k=89 —m=500 —left=trim_reads.For.paired.ec.fq.gz —right=trim_reads.Rev.paired.ec.fq.gz —unp=trim_reads.For.Rev.unp.ec.fq.gz
