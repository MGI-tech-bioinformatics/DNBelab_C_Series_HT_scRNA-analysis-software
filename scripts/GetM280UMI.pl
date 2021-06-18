#!/usr/bin/perl
use warnings;
#use strict;
use Getopt::Long;

my ($fq, $summary, $outfile, $beadBarcode);

GetOptions(
    "Fq:s" => \$fq,
    "Out:s" => \$outfile,
    "BeadBarcode:s" => \$beadBarcode,
);

if(!defined $fq || !defined $outfile || !defined $beadBarcode){
    print "\tperl $0 -BeadBarcode xxxx.txt -Fq xxx.fq -Out xxx.txt\n\n";
    exit(1);
}

my %bar;
open BAR, "$beadBarcode" || die "$!";
while(<BAR>){
    chomp;
    $bar{$_} = "";
}close BAR;

my(%hash, %count, %read);
my $num = 0; 
my $total = 0;
open IN, "$fq" || die "$!";
while(<IN>){
    chomp;
    $total+=1;
    my $r1_1 = $_;
    my $r1_2 = <IN>;
    my $r1_3 = <IN>;
    my $r1_4 = <IN>;
    
    my ($barcode) = $r1_1 =~ /CB:Z:(.*)/;
    my $umi = substr($r1_2, 0, 10);
    my $m280 = substr($r1_2, 10, 20);

    next if(!defined $bar{$barcode});
    $read{$barcode}{$m280} += 1;

    if(!defined $hash{$barcode}{$m280}{$umi}){
    	$hash{$barcode}{$m280}{$umi} = '';
        $count{$barcode}{$m280}+=1;
        $num+=1;
    }

}close IN;

open OUT, ">$outfile" || die "$!";
foreach my $barcode(keys %count){
    foreach my $m280(keys %{$count{$barcode}}){
		print OUT "$count{$barcode}{$m280}\t$m280\t$barcode\n";
    }
}
close OUT;
