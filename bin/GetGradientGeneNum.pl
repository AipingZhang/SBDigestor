#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
#use Parallel::ForkManager;
my $start = time;
my $usage=<<USAGE;

*******************************************************************************************************************************

Function: Significant insertional gene number calculation for each gradient for later saturation analysis.

Usage:	perl $0 [options]
	--bamDir <string>  bam file director;
	--genesbed <string> gene information with expected SB insertion probability;
	--MAPQ <int> quality thresold of the aligned reads[default:10] 
	--sampleName sample name;
	--h output help information;

Example: perl $0 --bamfile PE_reads.sort.bam --bedfile final.bed --OutPath ./		

***************************************************************************************************************************
USAGE
my($bamdir,$bed,$h,$dir,$cutoff,$sample,%hash);
GetOptions(
		'bamDir=s'=>\$bamdir,
		'genesbed=s'=>\$bed,
		'sampleName=s'=>\$sample,
		'MAPQ=i'=>\$cutoff,
		'h'=>\$h
	);  
die $usage if (!$bamdir or $h);
$cutoff=10 if(!$cutoff);
my $current=getcwd;

#my $sample=(split /\//,$bamdir)[-1];
`mkdir -p $current/Depth`;
`mkdir -p $current/Depth/GeneNum`;

#Loading genes information.
open OUT,">$current/Depth/GeneNum/$sample.txt"or die $!;
print OUT "Sample\tReadNum\tGeneNum\n";
open IN,"$bed" or die $!;
while(<IN>){
	chomp;
	my($chr,$start,$end,$strand,$gene,$fre)=(split /\s+/,$_)[0,1,2,3,4,5];
	$chr=~s/chr//;
	next if $chr!~/^(\d+)|X|MT|Y/ ;	
	$chr=$chr eq "X"?20:$chr eq "Y"?21:$chr eq "MT"?22:$chr;
	$hash{$chr}{"$start\t$end\t$strand\t$gene\t$fre"}=1;
}

close IN;

#Loading alignment results and determine mapping orientation.
for my $file (glob("$bamdir/*bam")){
	my %hash1=();
	my %hash2=();
	my $dirname=dirname($file);
	my $sampleName=(split /\//,$dirname)[-1];
	my $subsampleName=basename($file,'.bam');
	my $readNum=(split /_/,$subsampleName)[-1];
	open BAM,"samtools view $file |" or die $!;
	while (<BAM>){
		chomp;
		my($ID,$chr,$start,$mapq,$insertSize,$seq)=(split /\t/,$_)[0,2,3,4,8,9];
		next if $mapq<$cutoff;
		$chr=~s/chr//;
		next if $chr!~/^(\d+)|X|MT|Y/ ;
		$chr=$chr eq "X"?20:$chr eq "Y"?21:$chr eq "MT"?22:$chr;
		$hash1{$chr}{"$start\t$ID"}=1;
}

#Reads count for each gene.
for my $chr (sort {$a<=>$b}keys %hash)
	{for my $ele(keys %{$hash{$chr}})
		{my $count=0;
		my($start,$end,$strand,$gene,$fre)=(split /\t/,$ele)[0,1,2,3,4];
		for my $chr1(sort {$a<=>$b} keys %hash1)
			{next if $chr != $chr1 ;
			for my $con(keys %{$hash1{$chr1}})
				{
				my $pos=(split /\t/,$con)[0];
				if(($pos>$start)&&($pos<$end)){$count++;}
				}
			}
		
		if($count>0){$hash2{$gene}=1;}
	}

}

#Once a gene was identified in the previous all reads annotation results, we consider it is a true signal as long as the read count more than 0.
my $geneNum=0;
open RD,"$current/Annotation/$sampleName.txt"or die $!;
while(<RD>){
	chomp;
	my($gene)=(split /\t/,$_)[6];
	if(exists $hash2{$gene}){$geneNum++;}
}

print OUT "$sampleName\t$readNum\t$geneNum\n";

}
`sort -nk 2 $current/Depth/GeneNum/$sample.txt>$current/Depth/GeneNum/$sample\_1.txt`;
`mv $current/Depth/GeneNum/$sample\_1.txt $current/Depth/GeneNum/$sample.txt`;
my $duration = time - $start;
print "$duration\n";
close RD;
close OUT
