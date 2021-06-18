#!/usr/bin/perl -w
use warnings;
use Getopt::Long;
GetOptions(
	"c:s" => \$cell_stat,
	"m:s" => \$matrix,
	"o:s" => \$out,
);






open IN,"$cell_stat" || die "$!";
my %beads;
while(<IN>){
	chomp;
	my @a=split;
	$beads{$a[0]}{Raw}=$a[1];
	$beads{$a[0]}{UB}=$a[2];
}
close IN;

open IN,"$matrix" || die "$!";
my %cell;
while(<IN>){
	chomp;
	my @a=split;
	if(!exists $cell{$a[1]}){
		$cell{$a[1]}{Raw}=$beads{$a[0]}{Raw};
		$cell{$a[1]}{UB}=$beads{$a[0]}{UB};
	}	
	else{
		$cell{$a[1]}{Raw}+=$beads{$a[0]}{Raw};
		$cell{$a[1]}{UB}+=$beads{$a[0]}{UB};
	}
}
close IN;

open OUT,"> $out/merge_cell.stat";

print OUT "Cell\tRaw\tUB\n";
foreach my $key(keys %cell){
	print OUT "$key\t",$cell{$key}{Raw},"\t",$cell{$key}{UB},"\n";
}
close OUT;
