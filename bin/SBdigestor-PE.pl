#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd;
my $usage=<<USAGE;

*******************************************************************************************************************************
Program: SB Digestor
Version: Version1(06/01/2020)
Contact: yb97647\@umac.mo

Usage: perl SBdigestor.pl [options] 

Options:
	-i		<file> config file with sample name,enzyme type, fq1 and fq2;
	-o		<string> Output directory;
	-t		<int> thread number[default:10];
	-h		<string> output help information;	


Example: perl $0 -i fq.list -o outDir
******************************************************************************************************************************

USAGE

my($config,$outDir,$MAX_processes,$seqMethod,$h);

GetOptions(
		'i=s'=>\$config,
		'o=s'=>\$outDir,
		't=i'=>\$MAX_processes,
		'h'=>\$h

);

die $usage if(!$outDir && !$config);
$outDir=getcwd if (!$outDir);
$seqMethod="paired-end" if(!$seqMethod);
my $bin="/home/aipingzhang/project/Kai/newtool/FinalTest/bin_test";
my $database="/home/aipingzhang/project/Kai/newtool/FinalTest/database";

$MAX_processes=10 if(!$MAX_processes);
my $pm = Parallel::ForkManager->new($MAX_processes);
#my $pm1= Parallel::ForkManager->new($MAX_processes);
`mkdir -p $outDir/trim`;
`mkdir -p $outDir/align`;
`mkdir -p $outDir/align/all`;
`mkdir -p $outDir/align/gradient`;
`mkdir -p $outDir/Result`;

open IN, "$config" or die $!;
while(<IN>){
	my $pid = $pm->start and next;
	chomp;
	my($sampleName,$enzyme,$fq1,$fq2)=(split /\t/,$_)[0,1,2,3];
	#cut adapter and filter.
	`$bin/cutadapt -a GTCCCTTAAGCGGAGCCCTA  -g GTATGTAAACTTCCGACTTCAACTG -A CAGTTGAAGTCGGAAGTTTACATAC -G TAGGGCTCCGCTTAAGGGAC -m 20  -n 2 -pair-filter=both -o $outDir/trim/$sampleName\_R1.fastq.gz -p $outDir/trim/$sampleName\_R2.fastq.gz $fq1 $fq2`;
	#Calculating total reads number and randonly extract reads by 50 gradient for each sample.
	`perl $bin/GradientExtra.pl -sampleName $sampleName -fq1 $outDir/trim/$sampleName\_R1.fastq.gz -fq2 $outDir/trim/$sampleName\_R2.fastq.gz -o $outDir -bin $bin`;
	#alignment of all reads.
	`$bin/bowtie2 -x $database/Mus_musculus.GRCm38.dna.primary_assembly -1 $outDir/trim/$sampleName\_R1.fastq.gz -2 $outDir/trim/$sampleName\_R2.fastq.gz -S $outDir/align/all/$sampleName.sam`;
	`$bin/samtools sort -m 4G -o $outDir/align/all/$sampleName.bam $outDir/align/all/$sampleName.sam`;
	`$bin/samtools index $outDir/align/all/$sampleName.bam`;
	`$bin/samtools flagstat $outDir/align/all/$sampleName.bam > $outDir/align/all/$sampleName.flagstat`;
	`rm $outDir/align/all/$sampleName.sam`;
	#annotation and reads count for each gene; do binomial test to determine significant insertions(real SB insertions).
	`perl $bin/GetSigInsertions.pl --bamfile $outDir/align/all/$sampleName.bam --enzyme $enzyme --genesbed $database/finalRegion.txt --OutPath $outDir`;
	open POS,"$outDir/Annotation/$sampleName\_Pos.txt"or die $!;
	open OUT,">$outDir/Annotation/$sampleName.txt"or die $!;
	my %res=();
	while(<POS>){
		chomp;
		$res{$_}++;
		}
	close POS;
	foreach my $ele(keys(%res))
		{print OUT "$ele\t$res{$ele}\n";}
	close OUT;
	#alignment of each gradient for saturation evaluation and curve fitting.
	`mkdir -p $outDir/align/gradient/bowtie2/$sampleName`;
	open CON, "$outDir/align/gradient/config/$sampleName.align" or die $!;
	while (<CON>){
	#	my $pid1 = $pm1->start and next;	
		chomp;
		my ($subsample,$subfq1,$subfq2)=(split /\t/,$_)[1,2,3];
		`$bin/bowtie2 -x $database/Mus_musculus.GRCm38.dna.primary_assembly -1 $subfq1 -2 $subfq2 -S $outDir/align/gradient/bowtie2/$sampleName/$subsample.sam`; 
		`$bin/samtools sort -m 4G -o $outDir/align/gradient/bowtie2/$sampleName/$subsample.bam $outDir/align/gradient/bowtie2/$sampleName/$subsample.sam`;
		`$bin/samtools index $outDir/align/gradient/bowtie2/$sampleName/$subsample.bam`;
		`$bin/samtools flagstat $outDir/align/gradient/bowtie2/$sampleName/$subsample.bam >$outDir/align/gradient/bowtie2/$sampleName/$subsample.flastat`;
		`rm $outDir/align/gradient/bowtie2/$sampleName/$subsample.sam`;
		#$pm1->finish;
		}
		#$pm1->wait_all_children;
	close CON;
	#annotation for each gradient;
	`perl $bin/GetGradientGeneNum.pl   --bamDir $outDir/align/gradient/bowtie2/$sampleName --genesbed $database/finalRegion.txt --sampleName $sampleName --OutDir $outDir`;
	$pm->finish;	
}
$pm->wait_all_children;
close IN;

`cat $outDir/Depth/ReadCount/*txt >$outDir/Depth/ReadCount/readNum.txt`; 
#fitting curve by annotated genes and gradient reads number of each sample.
`perl $bin/FitCurve.pl $outDir`;

`perl $bin/CombineParaReadNum.pl  $outDir/Depth/formula/para.txt $outDir/Depth/ReadCount/readNum.txt $outDir/Depth/formula/para_readNum.txt`;
#calculating depth for each sample.
`perl $bin/DepthEvaluAndFilter.pl -i $outDir/Depth/formula/para_readNum.txt -o $outDir`;
#call common genes
`perl $bin/callCIS.pl $outDir/Result/matrix.txt`;
