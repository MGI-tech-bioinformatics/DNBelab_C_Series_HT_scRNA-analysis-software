#!/usr/bin/perl

=head1 Name

fishInWinter.pl  -- get out all the fishes which match to the given baits 

=head1 Description

Here I use "fishInWinter" to name this program, which do like fishing with baits.
The name is for memory of a college classmate met on Nankai bbs, and 
"fishInWinter" is just her bbs id. I wonder all the time what it means, 
a fish in the winter, or fishing in winter. Although all has passed, but
I just can't forget this bbs id "fishInWinter". 

This program works just like fishing, it helps you to get things that you 
wanted from a target file.  The things you wanted is stored in the bait file,
and the target file(fish file) is just like a fish pool. You can get out all 
the fishes which match the given baits. "Match" here has two meaning, same
or similar.

For the file formats, currently it supports fasta, table, gff, and fq. You can
specify the bait file format and the fish file format manually, or let the 
program detect it automatically. Note that table format is a default format,
you need to specify the key column for fishing, default is the first column.
For fasta, fq, and gff, it will use sequence id for fishing.

For the table format, the default seperator of each column is \s+, you can
change it inside the program if needed.

In the pattern mode, it will not require the same between bait and fish, but
work as:  fish =~ /bait/; 

In the gene mode, BM0001-TA,  BM0001-PA, and BM0001 will be take as the same
thing, because they belong to the same gene. "-TA" means transcprit A, "-PA"
means protein A. This is for convenience in our annotation system, you don't
need to care about it unless has the same demand as us.

The --except option will work as what it means, get all the fishes that not
match to baits. Both pattern mode and gene mode can be co-used with the --except 
option.

=head1 History

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 5.0,  Date: 2007-12-25

=head1 Usage
  
  % fishInWinter.pl <bait_file> <fish_file>
  --bformat <str>     set bait file format, fasta|table|gff|fq|fqgz
  --fformat <str>     set fish file format, fasta|table|gff|fq|fqgz
  --bcolumn <num>     set bait file column, default=1
  --fcolumn <num>     set fish file column, default=1
  --except            get things not in the bait file
  --patternmode       change to pattern mode, do not need exact same
  --genemode          change to gene mode, get things belonged to same gene
  --verbose           output running information to screen  
  --help              output help information to screen  

=head1 Exmple

  perl ./fishInWinter.pl -bf table -bc 10 -ff fasta  test-data/wanted.psl test-data/chr2.cds > test-data/chr2.cds.wanted
  perl ./fishInWinter.pl -bf table -bc 10 -ff gff test-data/wanted.psl test-data/chr2.gff > test-data/chr2.gff.wanted

  perl ./fishInWinter.pl -bf table -ff fasta -gene  test-data/needed.list test-data/chr2.cds > test-data/chr2.cds.needed
  perl ./fishInWinter.pl -bf table -ff fasta -gene -except test-data/needed.list test-data/chr2.cds > test-data/chr2.cds.not.needed


=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my $BaitSeperator = '\s+'; ## default seperator of bait file
my $FishSeperator = '\s+'; ## default seperator of fish file
my ($Baitformat,$FishFormat,$BaitColumn,$FishColumn,$Except,$PatternMode,$GeneMode);
my ($Verbose,$Help);
GetOptions(
	"bformat:s"=>\$Baitformat,
	"fformat:s"=>\$FishFormat,
	"bcolumn:s"=>\$BaitColumn,
	"fcolumn:s"=>\$FishColumn,
	"except"=>\$Except,
	"patternmode"=>\$PatternMode,
	"genemode"=>\$GeneMode,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$BaitColumn ||= 1;
$FishColumn ||= 1;
die `pod2text $0` if (@ARGV < 2 || $Help);

my $bait_file=shift;
my $fish_file=shift;

my %Bait;

if ($Baitformat eq "fasta" || (!$Baitformat && $bait_file=~/\.fa$/)) {
	read_fasta($bait_file,\%Bait,$BaitColumn,$BaitSeperator);
}elsif($Baitformat eq "gff" || (!$Baitformat && $bait_file=~/\.gff$/)) {
	read_gff($bait_file,\%Bait);
}elsif($Baitformat eq "fq" || (!$Baitformat && $bait_file=~/\.fq$/)) {
	read_fq($bait_file,\%Bait,$BaitColumn,$BaitSeperator);
}elsif($Baitformat eq "fqgz" || (!$Baitformat && $bait_file=~/\.fq\.gz$/)) {
	read_fqgz($bait_file,\%Bait,$BaitColumn,$BaitSeperator);
}else{
	read_table($bait_file,\%Bait,$BaitColumn,$BaitSeperator);
}

##print Dumper \%Bait;

print STDERR "read bait done\n" if ($Verbose);

if ($FishFormat eq "fasta" || (!$FishFormat && $fish_file=~/\.fa$/)) {
	output_fasta($fish_file,\%Bait,$FishColumn,$FishSeperator);
}elsif($FishFormat eq "gff" || (!$FishFormat && $fish_file=~/\.gff$/)) {
	output_gff($fish_file,\%Bait);
}elsif($FishFormat eq "fq" || (!$FishFormat && $fish_file=~/\.fq$/)) {
	output_fq($fish_file,\%Bait,$FishColumn,$FishSeperator);
}elsif($FishFormat eq "fqgz" || (!$FishFormat && $fish_file=~/\.fq\.gz$/)) {
	output_fqgz($fish_file,\%Bait,$FishColumn,$FishSeperator);
}else{
	output_table($fish_file,\%Bait,$FishColumn,$FishSeperator);
}

print STDERR "Fish out done\n" if ($Verbose);


####################################################
################### Sub Routines ###################
####################################################


sub read_fq{
	my ($file,$bait_hp,$bait_colum,$bait_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		if (/^\@(\S+)/) {
			my $id = $1;
			<IN>;<IN>;<IN>;
			my @temp=split(/$bait_seperator/,$id);
			$$bait_hp{$temp[$bait_colum-1]}=1;	## the fq file do not has gene mode, in fact
		}
	}
	close(IN);
}

sub read_fqgz{
	my ($file,$bait_hp,$bait_colum,$bait_seperator) = @_;
#	open(IN,$file)||die("fail to open $file\n");
	open IN,"gzip -dc $file|"||die("fail to open $file\n");
	while (<IN>) {
		if (/^\@(\S+)/) {
			my $id = $1;
			<IN>;<IN>;<IN>;
			my @temp=split(/$bait_seperator/,$id);
			$$bait_hp{$temp[$bait_colum-1]}=1;	## the fq file do not has gene mode, in fact
		}
	}
	close(IN);
}
sub read_gff{
	my ($file,$bait_hp) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/\t/,$_);
		my $id = $2 if($temp[8] =~ /(ID|Parent)=([^;]+);/);
		if (!$GeneMode) {
			$$bait_hp{$id}=1;
		}else{
			$id = $1 if($id =~ /^(\w+)-\w+$/);
			$$bait_hp{$id} = 1;
		}
		
	}
	close(IN);
	
}


sub read_table{
	my ($file,$bait_hp,$bait_colum,$bait_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/$bait_seperator/,$_);
		my $id = $temp[$bait_colum-1];

		if (!$GeneMode) {
			$$bait_hp{$id}=1;
		}else{
			$id = $1 if($id =~ /^(\w+)-\w+$/);
			$$bait_hp{$id} = 1;
			
		}
		
	}
	close(IN);
}


sub read_fasta{
	my ($file,$bait_hp,$bait_colum,$bait_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		chomp $title;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/$bait_seperator/,$title);
		my $id = $temp[$bait_colum-1];

		if (!$GeneMode) {
			$$bait_hp{$id}=1;
		}else{
			$id = $1 if($id =~ /^(\w+)-\w+$/);
			$$bait_hp{$id} = 1;
		}
	}
	close(IN);
}

## do not have gene mode in fq format
sub output_fq{
	my ($file,$fish_hp,$fish_colum,$fish_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		if (/^\@(\S+)/) {
			my $id = $1;
			my $content = $_;
			$content .= <IN>;
			$content .= <IN>;
			$content .= <IN>;

			my @temp=split(/$fish_seperator/,$id);
			
			if (!$Except) {
				print $content if (!$PatternMode && exists $$fish_hp{$temp[$fish_colum-1]});
				print $content if ($PatternMode && word_pattern_hash($temp[$fish_colum-1],$fish_hp));

			}else{
				print $content if (!$PatternMode && !exists $$fish_hp{$temp[$fish_colum-1]});
				print $content if ($PatternMode && !word_pattern_hash($temp[$fish_colum-1],$fish_hp));
			}

		}
		
	}
	close(IN);
}

sub output_fqgz{
	my ($file,$fish_hp,$fish_colum,$fish_seperator) = @_;
#	open(IN,$file)||die("fail to open $file\n");
	open IN,"gzip -dc $file|"||die("fail to open $file\n");
	while (<IN>) {
		if (/^\@(\S+)/) {
			my $id = $1;
			my $content = $_;
			$content .= <IN>;
			$content .= <IN>;
			$content .= <IN>;

			my @temp=split(/$fish_seperator/,$id);
			
			if (!$Except) {
				print $content if (!$PatternMode && exists $$fish_hp{$temp[$fish_colum-1]});
				print $content if ($PatternMode && word_pattern_hash($temp[$fish_colum-1],$fish_hp));

			}else{
				print $content if (!$PatternMode && !exists $$fish_hp{$temp[$fish_colum-1]});
				print $content if ($PatternMode && !word_pattern_hash($temp[$fish_colum-1],$fish_hp));
			}

		}
		
	}
	close(IN);
}

sub output_gff{
	my ($file,$fish_hp) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/\t/,$_);
		my $id = $2 if($temp[8] =~ /(ID|Parent)=([^;]+);/);
		$id = $1 if($GeneMode && $id =~ /^(\w+)-\w+$/);

		if (!$Except) {
			print $_."\n" if (!$PatternMode && exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && word_pattern_hash($id,$fish_hp));

		}else{
			print $_."\n" if (!$PatternMode && !exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && !word_pattern_hash($id,$fish_hp));
		}
		
	}
	close(IN);
}


sub output_table{
	my ($file,$fish_hp,$fish_colum,$fish_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//g;
		my @temp=split(/$fish_seperator/,$_);
		my $id = $temp[$fish_colum-1];
		$id = $1 if($GeneMode && $id =~ /^(\w+)-\w+$/);

		if (!$Except) {
			print $_."\n" if (!$PatternMode && exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && word_pattern_hash($id,$fish_hp));

		}else{
			print $_."\n" if (!$PatternMode && !exists $$fish_hp{$id});
			print $_."\n" if ($PatternMode && !word_pattern_hash($id,$fish_hp));
		}
		
	}
	close(IN);
}


sub output_fasta{
	my ($file,$fish_hp,$fish_colum,$fish_seperator) = @_;
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/$fish_seperator/,$title);
		my $id = $temp[$fish_colum-1];
		$id = $1 if($GeneMode && $id =~ /^(\w+)-\w+$/);
		
		if (!$Except) {
			print ">".$title.$seq if (!$PatternMode && exists $$fish_hp{$id});
			print ">".$title.$seq if ($PatternMode && word_pattern_hash($id,$fish_hp));

		}else{
			print ">".$title.$seq if (!$PatternMode && !exists $$fish_hp{$id});
			print ">".$title.$seq if ($PatternMode && !word_pattern_hash($id,$fish_hp));
		}
	}
	close(IN);
}


sub word_pattern_hash{
	my $word = shift;
	my $hash_p = shift;
	my $value = 0;
	foreach my $hash_key (keys %$hash_p) {
		if ($word =~ /$hash_key/){
			$value = 1;
			last;
		}
	}
	return $value;
}
