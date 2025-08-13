#!/usr/bin/perl -w
use strict;
#get tRNA id and seq;
my $name;
my %tRNAid_seq;
my %hash;
open (IN,$ARGV[0]);
while(<IN>){
chomp;
if (/^>/){
my $ss=$_;my @dd=split/\s+/,$ss;
$name=$dd[0];
$name=~s/>//;
$ss=~s/>//;
$hash{$name}=$ss;
}
else{
$tRNAid_seq{$name}=$_;
}
}
close IN;


#get tRNA_id=> tag pos
my %tRNA_tag;
open (INN,$ARGV[1]);
while(<INN>){
chomp;
my @d=split(/\t/);
if ($d[3]==$d[5] && $d[4]==$d[5] && $d[4]/$d[5]>0.9 && $d[9]<$d[10]){
my $tag={'tagId',$d[0],'beg',$d[9],'end',$d[10],'strand',"+"};
push @{$tRNA_tag{$d[1]}},$tag;
}
}
close INN;

#get tag information;
my $tag_file=$ARGV[2];
my %tagId_seq;
my %tagId_number;
&read_tag_file($tag_file,\%tagId_seq,\%tagId_number);


#print aln
foreach my $id (sort keys %tRNA_tag){
#if (exists $short_tRNAid{$id}){
print "\>$hash{$id}\n";
print "$tRNAid_seq{$id}\n";
my $len=length $tRNAid_seq{$id};
my @tag=@{$tRNA_tag{$id}};
foreach my $tag (sort {$a->{beg}<=>$b->{beg}} @tag){
my $tagstrand=$tag->{strand};
my $tagId=$tag->{tagId};
my $tagseq=$tagId_seq{$tagId};
my $taglen=length $tagseq;
my $number=$tagId_number{$tagId};
my $str="." x $len;
substr($str,$tag->{beg}-1,$taglen)=$tagseq;
print "$str $tagId $taglen $number $tagstrand $tag->{beg}\n";
#}
}
}

sub read_tag_file{
my $infile=shift;
my $tagId_tag=shift;
my $tagId_number=shift;
open NN,$infile || die $!;
while(<NN>){
chomp;
if (/^>(\S+)\s+(\d+)/){
my $id=$1;
my $number=$2;
my $tag=<NN>;
chomp $tag;
$tagId_tag->{$id}=$tag;
$tagId_number->{$id}=$number;
}
}
close NN;
}
