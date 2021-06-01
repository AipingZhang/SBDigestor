#!/usr/bin/perl -w
use strict;
my $infile1=shift;
die "perl $0 <infile1>\n" unless defined $infile1;
my %hash;
open IN1,"$infile1"or die $!;
#my $header=<IN1>;

while(<IN1>){
        chomp;
        my($sample,$gene)= (split /\t/,$_)[0,5];
        print "$gene\n";
	$hash{$gene}{$sample}=1;
}

close IN1;

open OUT, ">CIS.txt" or die $!;
for my $gene(keys %hash)
	{my $count=0;
	print OUT "$gene\t";
	for my $sample (keys %{$hash{$gene}})
		{$count++;print OUT "$sample,";}
		
	print OUT "\t$count\n";
	}
close OUT;
