#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Statistics::R;
my $usage=<<USAGE;

*******************************************************************************************************************************

Function: Gene Annotation And Read Count

Usage:	perl $0 [options]
	--bamfile <string> bam file;
	--genesbed <string> gene information with expected SB insertion probability;
	--seqMethod <string> single-end or paired-end sequencing[paired-end];
	--enzyme type <string> N or B;
	--MAPQ <int> quality thresold of the aligned reads[default:10];
	--pvalueThresold<int> p value of binomial distribution test for the reads number of each gene[defalt:0.05]; 
	--OutPath <string> Output direction[default:./];
	--h output help information;

Example: perl $0 --bamfile PE_reads.sort.bam --enzyme N --bedfile genes.bed --OutPath ./		

***************************************************************************************************************************
USAGE

my($bam,$bed,$h,$dir,$cutoff,$outpath,%hash,%hash1,$p_cutoff,$seqMethod,$enzyme);
GetOptions(
		'bamfile=s'=>\$bam,
		'genesbed=s'=>\$bed,
		'seqMethod=s'=>\$seqMethod,
		'enzyme=s'=>\$enzyme,
		'MAPQ=i'=>\$cutoff,
		'pvalueThresold=i'=>\$p_cutoff,
		'OutPath=s'=>\$outpath,
		'h'=>\$h
	);  
die $usage if (!$bam or $h);
$cutoff=10 if(!$cutoff);
$p_cutoff=0.05 if(!$p_cutoff);
$seqMethod="paired-end" if(!$seqMethod);
my $current=getcwd;
$outpath=$current if(!$outpath);
my $sampleName=(split /\//,$bam)[-1];
$sampleName=(split /\./,$sampleName)[0];
#$enzyme=(split /_/,$sampleName)[-1];
`mkdir -p  $current/Annotation`;

#Loading genes information.
open IN,"$bed" or die $!;
while(<IN>){
	chomp;
	my($chr,$start,$end,$strand,$gene,$fre)=(split /\s+/,$_)[0,1,2,3,4,5];
	$chr=~s/chr//;
	next if $chr!~/(\d+)|X|MT|Y/ ;	
	$chr=$chr eq "X"?20:$chr eq "Y"?21:$chr eq "MT"?22:$chr;
	$hash{$chr}{"$start\t$end\t$strand\t$gene\t$fre"}=1;
}

close IN;

#Loading alignment results and determine mapping orientation.
my $mapped;
open BAM,"samtools view  -F 4 $bam |" or die $!;
while (<BAM>){
	chomp;
	my $strand;
	my($ID,$flag,$chr,$start,$mapq,$insertSize,$seq)=(split /\t/,$_)[0,1,2,3,4,8,9];
	next if $mapq<$cutoff;
	$chr=~s/chr//;
	next if $chr!~/^(\d+)|X|MT|Y/ ;
	$chr=$chr eq "X"?20:$chr eq "Y"?21:$chr eq "MT"?22:$chr;
	if($seqMethod eq "paired-end"){
		if($enzyme eq "N"){
 			if($flag=~ /147|99/){$strand="+";}
			elsif($flag=~/83|163/){$strand="-";}	
			else{next;}
		}
		else{	
			if($flag=~ /147|99/){$strand="-";}
			elsif($flag=~/83|163/){$strand="+";}
			else{next;}
		}
	}
	else{
		if($enzyme eq "N"){
			if($flag=~ /0/){$strand="+";}
			elsif($flag=~/16/){$strand="-";}
			else{next;}
				}
		else{
			if($flag=~ /0/){$strand="-";}
			elsif($flag=~ /16/){$strand="+";}
			else{next;}

			}
		}
		
	next if $chr!~/^(\d+)/ ;
	$mapped++;
	$hash1{$chr}{"$start\t$ID\t$strand"}=1;
}
close BAM;

#Do binomial test to determine significant insertions(real SB insertions).
open OUT,">$current/Annotation/$sampleName\_count.txt"or die $!;
open O,">$current/Annotation/$sampleName\_Pos.txt"or die $!;
print O "chr\treadpos\treadstrand\tstart\tend\tgene\tgenestrand\n";
print OUT "chr\tstart\tend\tstrand\tgene\tExpected_p\tcount\tP_value\n";

for my $chr (sort {$a<=>$b}keys %hash)
	{for my $ele(keys %{$hash{$chr}})
		{my $count=0; my %hash2=();
		my($start,$end,$strand,$gene,$fre)=(split /\t/,$ele)[0,1,2,3,4];
		for my $chr1(sort {$a<=>$b} keys %hash1)
			{next if $chr != $chr1 ;
			for my $con(keys %{$hash1{$chr1}})
				{
				my ($pos,$strand)=(split /\t/,$con)[0,2];
				if(($pos>$start)&&($pos<$end)){$count++;$hash2{$gene}{"$con"}=1;}
				}
			}
		
		if($count>0){	my $R = Statistics::R->new( bin => '/home/aipingzhang/anaconda2/bin/R' );
				
				$R->set('x', $count);
				$R->set('n', $mapped);
				$R->set('p', $fre);
				$R->run(qq`library(stats)`);
				my$bionm=$R->run(q`binom.test(x=x,n=n,p=p)`);
				my $pvalue;
				my ($value1,$value2)=(split /\n/,$bionm)[4,5];
				if($value2=~/^\d+/){ $pvalue=$value2;}
				else{
				my $p=(split /,/,$value1)[-1];
				$pvalue=(split /=|>|</,$p)[-1];}
				$pvalue=sprintf ("%.8f",$pvalue);
				if($pvalue < $p_cutoff){
					for my $g (keys %hash2)
						{for my $pos_strand (keys %{$hash2{$g}})
							{my ($read_pos,$read_strand)=(split /\t/,$pos_strand)[0,2];
							print O "$sampleName\t$chr\t$read_pos\t$read_strand\t$start\t$end\t$gene\t$strand\n";}
						}
					print OUT "$chr\t$ele\t$count\t$pvalue\n";}
	
				}
	}

}
close OUT;
