1. download sra
prefetch --option-file SRA.list
2. sra to fastq
fastq-dump $1.sra --gzip --split-3 -O /work/home/shuziqiang/xiaotai/test_miRNA/data/fastq
3. adapter trimming
java -cp /work/home/shuziqiang/xiaotai/test_miRNA/0.soft/TBtools_JRE1.6.jar biocjava.sRNA.Tools.sRNAseqAdaperRemover --minLen 17 --inFxFile /work/home/shuziqiang/xiaotai/test_miRNA/data/fastq/$1.fastq.gz --outFaFile /work/home/shuziqiang/xiaotai/test_miRNA/2.analysis/01.filter/$1.fa
