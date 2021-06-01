#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use Parallel::ForkManager;
my $usage=<<USAGE;

*************************************************************************************************************************************************
Program: Depth evaluating and filtering

Usage: perl DepthEvaluAndFilter.pl 

Options:
	-i	<file> The fitting curve parameter and read number for each sample; 
	-k	<int> The formula canstant[default:1];
	-o	<string> The output directory[default:current];
	-h	<string> output help information;

Example:perl $0 -i para_readNum.txt -k 1 

**************************************************************************************************************************************************
 
USAGE

my($para_num,$outDir,$k,$h);

GetOptions(
		'i=s'=>\$para_num,
		'o=s'=>\$outDir,
		'k=i'=>\$k,
		'h'=>\$h

);

die $usage if(!$outDir && !$para_num);
$outDir=getcwd if (!$outDir);
$k=1 if(!$k);

my $depth;my %hash;my %hash1;
`mkdir -p $outDir/Result`;

open IN,"$para_num" or die $!;
open D,">$outDir/Result/depth.txt"or die $!;
while(<IN>){
	chomp;
	next if $_=~/^Sample/;
	my($sample,$a,$b,$readNum)=(split /\t/,$_)[0,1,2,4];
	my $depth=int($readNum**(1-$b)/($a*$k));
	print D "$sample\t$depth\n";
	$hash{$sample}=$depth;
}
close IN;
open S,">$outDir/Result/sampleGene.txt"or die $!;
open M,">$outDir/Result/matrix.txt"or die $!;
open P,">$outDir/Result/result.txt"or die $!;


for my $file (glob("$outDir/Annotation/*_count.txt")){
	my $sample=basename($file,'_count.txt');
	if (exists $hash{$sample})
		{	
		open IN1,"$file"or die $!;
		my $header=<IN1>;
		while(<IN1>){
			chomp;
			my ($gene,$count)= (split /\t/,$_)[4,6];
			if ($count>$hash{$sample})
				{print M "$sample\t$_\n";print S "$sample\t$gene\n";
			open IN2,"$outDir/Annotation/$sample.txt"or die $!;
			while(<IN2>){
				chomp;
				my $line=$_;
				my($sample,$chr,$r_pos,$r_strand,$start,$end,$gene1,$strand,$num)=(split /\s+/,$line)[0,1,2,3,4,5,6,7,8];
				if($gene1 eq $gene){print P "$sample\t$chr\t$r_pos\t$r_strand\t$start\t$end\t$gene1\t$strand\t$num\n";}
				else{next;}	
				}
			close IN2;
				}
			}			
				
		}

		close IN1;
	
}

close S;
close M;
close P;
