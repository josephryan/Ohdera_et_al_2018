A. Sequenced Reads Trimming and Quality Filtering

A1. FASTQ files generated from IlluminaHiseq2000 sequencer were trimmed to remove adapters and filtered by base quality. 
Dependencies: Trimmomatic v0.32 (Bolger et al. 2014)

java -jar /usr/local/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 32 -phred33 Lcru.R1.fq Lcru.R2.fq Lcru.R1.trimmo.fq Lcru.unp1.fq Lcru.R2.trimmo.fq Lcru.unp2.fq ILLUMINACLIP:/usr/local/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:12:1:true MINLEN:36

A2. Error Correction with ErrorCorrectReads.pl from Allpaths-LG

Dependencies: Allpaths-LG version 52488 (Gnerre et al. 2011) 


perl /usr/local/allpathslg-44837/src/ErrorCorrectReads.pl PAIRED_READS_A_IN=../01-TRIMMOMATIC/Lcru.R1.trimmo.fq.gz PAIRED_READS_B_IN=../01-TRIMMOMATIC/Lcru.R2.trimmo.fq.gz UNPAIRED_READS_IN=../01-TRIMMOMATIC/Lcru.unp1.fq.gz PAIRED_SEP=100 THREADS=30 PHRED_ENCODING=33 READS_OUT=Lcru

B. Mitochondria Reads Removal

Dependencies: FastqSifter v0.03 (formerly remove_mt v0.02 - https://github.com/josephryan/FastqSifter)

remove_mt --out=CcruxMT --mt=CcruxMT.fa --left=Lcru.A.fq --right=Lcru.B.fq --unp=Lcru.unp.fq

C. De novo Assembly using Platanus

Dependencies: Platanus v1.2.1 (Kajitani et al. 2014)
	      plat.pl

- 32 kmer
plat.pl --threads=20 --out=Ccrux_MTfiltered.32 --k=32 --m=500 --left=../04-REMOVE_MT/CcruxMT.mt_filtered.A.fq --right=../04-REMOVE_MT/CcruxMT.mt_filtered.B.fq --unp=../04-REMOVE_MT/CcruxMT.mt_filtered.unp.fq 

-45 kmer
plat.pl --threads=20 --out=Ccrux_MTfiltered.45 --k=45 --m=500 --left=../04-REMOVE_MT/CcruxMT.mt_filtered.A.fq --right=../04-REMOVE_MT/CcruxMT.mt_filtered.B.fq --unp=../04-REMOVE_MT/CcruxMT.mt_filtered.unp.fq

D. Construct artificial metaphor libraries using MateMaker

Dependencies: MateMaker (https://github.com/josephryan/matemaker)

- Suboptimal Assembly 1 (kmer = 45)  
matemaker --matesperkb=50 --readlen=50 --insertsize=1000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.1000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=2000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.2000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=3000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.3000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=4000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.4000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=5000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.5000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=7500 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.7500.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=10000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.10000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=15000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.15000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=20000 --assembly=./03-PLAT45.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat45k.20000.mpkb

- Suboptimal Assembly 2 (kmer = 32)  
matemaker --matesperkb=50 --readlen=50 --insertsize=1000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.1000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=2000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.2000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=3000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.3000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=4000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.4000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=5000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.5000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=7500 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.7500.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=10000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.10000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=15000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.15000.mpkb
matemaker --matesperkb=50 --readlen=50 --insertsize=20000 --assembly=./03-PLAT32.NO_MT_CLEAN/out_gapClosed.fa --out=Lc.plat32k.20000.mpkb

E. Scaffolding Platanus assembly (kmer = 45) with artificial matepair libraries.

Dependencies: SSPACE VERSION: SSPACE_Standard_v3.0_linux (Boetzer et al. 2011)

SSPACE -l libraries.txt -s Ccrux.01.fa -b run

# LIBRARIES.TXT used by SSPACE (/data1/jfryan/03-Lucernariopsis/25-MFLS/01-MULTIPLE_LIBS/run3)
Lib1 bowtie Lc.plat32k.1000.mpkb.A.fq Lc.plat32k.1000.mpkb.B.fq 1000 0.25 RF
Lib2 bowtie Lc.plat32k.2000.mpkb.A.fq Lc.plat32k.2000.mpkb.B.fq 2000 0.25 RF
Lib3 bowtie Lc.plat32k.3000.mpkb.A.fq Lc.plat32k.3000.mpkb.B.fq 3000 0.25 RF
Lib4 bowtie Lc.plat32k.4000.mpkb.A.fq Lc.plat32k.4000.mpkb.B.fq 4000 0.25 RF
Lib5 bowtie Lc.plat32k.5000.mpkb.A.fq Lc.plat32k.5000.mpkb.B.fq 5000 0.25 RF
Lib6 bowtie Lc.plat32k.7500.mpkb.A.fq Lc.plat32k.7500.mpkb.B.fq 7500 0.25 RF
Lib7 bowtie Lc.plat32k.10000.mpkb.A.fq Lc.plat32k.10000.mpkb.B.fq 10000 0.25 RF
Lib8 bowtie Lc.plat32k.15000.mpkb.A.fq Lc.plat32k.15000.mpkb.B.fq 15000 0.25 RF
Lib9 bowtie Lc.plat32k.20000.mpkb.A.fq Lc.plat32k.20000.mpkb.B.fq 20000 0.25 RF
Lib10 bowtie Lc.plat48k.1000.mpkb.A.fq Lc.plat48k.1000.mpkb.B.fq 1000 0.25 RF
Lib11 bowtie Lc.plat48k.2000.mpkb.A.fq Lc.plat48k.2000.mpkb.B.fq 2000 0.25 RF
Lib12 bowtie Lc.plat48k.3000.mpkb.A.fq Lc.plat48k.3000.mpkb.B.fq 3000 0.25 RF
Lib13 bowtie Lc.plat48k.4000.mpkb.A.fq Lc.plat48k.4000.mpkb.B.fq 4000 0.25 RF
Lib14 bowtie Lc.plat48k.5000.mpkb.A.fq Lc.plat48k.5000.mpkb.B.fq 5000 0.25 RF
Lib15 bowtie Lc.plat48k.7500.mpkb.A.fq Lc.plat48k.7500.mpkb.B.fq 7500 0.25 RF
Lib16 bowtie Lc.plat48k.10000.mpkb.A.fq Lc.plat48k.10000.mpkb.B.fq 10000 0.25 RF
Lib17 bowtie Lc.plat48k.15000.mpkb.A.fq Lc.plat48k.15000.mpkb.B.fq 15000 0.25 RF
Lib18 bowtie Lc.plat48k.20000.mpkb.A.fq Lc.plat48k.20000.mpkb.B.fq 20000 0.25 RF




