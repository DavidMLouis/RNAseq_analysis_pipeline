# RNAseq_analysis_pipeline
Modified and updated version analysis from FASTQ to Differential Expression


This is an updated version of the 'RNA-seq: Basic Bioinformatics Analysis' paper by Fei Ji and Ruslan Sadreyev.

Download and install the required tools:
STAR: https://github.com/alexdobin/STAR
Picard: https://broadinstitute.github.io/picard/
HTseq: https://htseq.readthedocs.io/en/release_0.9.1/install.html R: https://www.r-project.org

Download Reference set of DNA and GTF files:
https://useast.ensembl.org/info/data/ftp/

1. Create reference database
STAR --genomeDir /Users/kingdavid/Desktop/stanford/stem_transplant_GB/RNA_seq/human_GRCh38 \
--readFilesIn /Users/kingdavid/Desktop/stanford/stem_transplant_GB/RNA_seq/data/SRR18938586_2.fastq \
--sjdbGTFfile /Users/kingdavid/Desktop/stanford/stem_transplant_GB/RNA_seq/ref/Homo_sapiens.GRCh38.106.gtf\
--genomeLoad NoSharedMemory --outReadsUnmapped Fastx

2. ./RNAseq_align.pl sample_list.txt configure.txt
3. ./RNAseq_qc.pl sample_list.txt configure.txt
4. ./RNAseq_count.pl sample_list.txt configure.txt

Then follow the sequences from:

5. RNAseq.r
