#!/usr/bin/perl -w
use strict;
my $mi=$ARGV[0];#miRNA seq in fa format
my $tr=$ARGV[1];#tRNA seq in fa format
system("/NAS/zl/software/blast-2.2.16/bin/blastall -p blastn -a 32 -i $mi -d $tr -F F -o ./map2tRNA.blast");
system("perl ./blast2table.pl map2tRNA.blast > map2tRNA.tab");
system("perl ./filter-tab.pl map2tRNA.tab > map2tRNA.filter");
system("perl ./aln-v1.pl tRNA.fa ./map2tRNA.filter miRNA.fa > map2tRNA.aln");



