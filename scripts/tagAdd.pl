#!/usr/bin/perl -w
use warnings;
#use strict;
use Getopt::Long;

my ($bam,$file,$out,$tag_check,$tag_add,$debug,$root);

GetOptions(
   "bam:s" => \$bam,
   "file:s" => \$file,
   "out:s" => \$out,
   "tag_check:s" => \$tag_check,
   "tag_add:s" => \$tag_add,
   "debug" => \$debug,
    "root:s" => \$root,
);

$tag_check ||= "CB:Z:";
$tag_add   ||= "DB:Z:";

if(!defined $bam || !defined $file || !defined $out){
   print "\tperl $0 -bam xxx.bam -file xxx.txt -out xxx.add.bam -tag_check CB:Z: -tag_add DB:Z:\n\n";
   exit(1);
}

my %barcode_ref;
open FILE, "$file" || die "$!";
while(<FILE>){
   chomp;
   my($a,$b) = split(/\s+/,$_);
   $barcode_ref{$a} = $b;
}close FILE;

my $NUM=0;
open OUT,"| samtools view -bS >$out" || die "$!";
open BAM,"samtools view -h $bam |" || die "$!";
while(<BAM>){
   chomp;
   if($_ =~ /^@/){print OUT "$_\n";next}
    $NUM++;
   my @a = split(/\s+/,$_);
   foreach my $a(@a){
      if($a =~ $tag_check){
         my ($tag) = $a =~ /$tag_check(.*)/;
         if(defined $barcode_ref{$tag}){print OUT "$_\t$tag_add$barcode_ref{$tag}\n";}
      }
   }
last if($debug && $NUM == 1000);
}close BAM;
close OUT; 
