#!/usr/bin/perl -w

use strict;
my $infile1=shift;
my $infile2=shift;
my $out=shift;
die "perl $0 <parameter> <ReadNum> <outfile>\n" unless defined $infile1;
my %hash;
open OUT,">$out"or die $!;
print OUT "Sample\ta\tb\tr-square\treadNum\n";
open IN1,"$infile1"or die $!;
while(<IN1>){
	chomp;
	next if $_ =~/^sample/;
	my($sample,$a,$b,$r)= (split /\s+/,$_)[0,2,3,4];
	$hash{$sample}="$sample\t$a\t$b\t$r";}
close IN1;
open IN2,"$infile2"or die $!;
while (<IN2>){
        chomp;
	my($sample,$readNum)= (split /\t/,$_)[0,1];
        if (exists $hash{$sample}){print OUT "$hash{$sample}\t$readNum\n";}
}
close IN2;	
close OUT;
