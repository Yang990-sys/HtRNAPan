#!/usr/bin/perl -w
use strict;
open(IN,$ARGV[0]);

while(<IN>){
chomp;
my @a=split/\t/;

if ($a[3]==$a[5] && $a[4]==$a[5]){
print "$_\n";
}
}
close IN;


