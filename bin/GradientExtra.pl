#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use Cwd;
my $usage=<<USAGE;

************************************************************************************************************************************

Function: Calculating total reads number and randonly extract reads by 50 gradient for each sample.  

Usage: perl readsCountAndExtract.pl  [options]
Options:
	-sampleName	<string> Sample Name;
	-fq1		<string> The clean fq1;
	-fq2            <string> The clean fq2;
	-o		<string> Output directory;
	-bin		<string> Program directory;
	-h		<string> Output help information;

Example: perl $0 -sampleName -fq1 fq1 -fq2 fq2 -o outDir -bin ./

*************************************************************************************************************************************

USAGE

my($sampleName,$fq1,$fq2,$outDir,$bin,$h);

GetOptions(
	'sampleName=s'=>\$sampleName,
	'fq1=s'=>\$fq1,
	'fq2=s'=>\$fq2,
	'o=s'=>\$outDir,
	'bin=s'=>\$bin,
	'h=s'=>\$h
);

die $usage if(!$outDir && !$bin);
$outDir=getcwd if(!$outDir);

`mkdir -p $outDir/ExtraReads`;
`mkdir -p $outDir/Depth`;
`mkdir -p $outDir/Depth/ReadCount`;
`mkdir -p $outDir/align`;
`mkdir -p $outDir/align/gradient`;
`mkdir -p $outDir/align/gradient/config`;
`mkdir -p $outDir/align/gradient/bowtie2`;
open OUT1, ">$outDir/Depth/ReadCount/$sampleName.txt"or die $!;
`mkdir -p $outDir/ExtraReads/$sampleName`;
	my $read1_counter=0;
	my $read2_counter=0;
	if($fq1=~/gz$/){
		open FQ1, "gunzip -dc $fq1|" or die $!;
	}else{
		open FQ1, "$fq1" or die $!;} 

	if($fq2=~/gz$/){
		open FQ2, "gunzip -dc $fq2|" or die $!;
	}else{
		open FQ2, "$fq2" or die $!;}

	while(<FQ1>){
		$read1_counter++;
		}
	close FQ1;
	while(<FQ2>){
		$read2_counter++;
		}
	close FQ2;
	my $line=$read1_counter+$read2_counter;
	my $readnum=$line/4;
	my $readpair=$readnum/2;
	print OUT1 "$sampleName\t$readpair\n";
	
	if($readpair>50){
		my $average=int($readpair/50);
		my @array=();
		my $i;
		for ( $i=$average;$i<$readpair;$i=$i+$average){
			my $gradient=$i;
			print "$gradient\n";
			`$bin/seqtk-master/seqtk sample -s 100 $fq1 $gradient > $outDir/ExtraReads/$sampleName/$sampleName\_$gradient\_R1.fq`;
			`gzip $outDir/ExtraReads/$sampleName/$sampleName\_$gradient\_R1.fq`;
			`$bin/seqtk-master/seqtk sample -s 100 $fq2 $gradient > $outDir/ExtraReads/$sampleName/$sampleName\_$gradient\_R2.fq`;
			`gzip $outDir/ExtraReads/$sampleName/$sampleName\_$gradient\_R2.fq`;
			}
		}
		
`mkdir -p $outDir/ExtraReads/$sampleName/$sampleName`;
`find $outDir/ExtraReads/$sampleName -name "*R1.fq.gz"|sort >$outDir/ExtraReads/$sampleName/$sampleName/fq1.list`;
`find $outDir/ExtraReads/$sampleName -name "*R2.fq.gz"|sort >$outDir/ExtraReads/$sampleName/$sampleName/fq2.list`;
`paste $outDir/ExtraReads/$sampleName/$sampleName/fq1.list $outDir/ExtraReads/$sampleName/$sampleName/fq2.list>$outDir/ExtraReads/$sampleName/$sampleName/fq.list`;

open IN,"$outDir/ExtraReads/$sampleName/$sampleName/fq.list"or die $!;
open OUT2,">$outDir/align/gradient/config/$sampleName.align"or die $!;
while(<IN>){
	chomp;
	my $col1=(split /\t/,$_)[0];
	my $subSample=basename($col1,'_R1.fq.gz');
	print OUT2 "$sampleName\t$subSample\t$_\n";
}

close IN;
close OUT1;
close OUT2;
