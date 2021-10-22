#!/usr/bin/env perl
use Getopt::Long;
use strict;
use warnings;

my $outdir = "snpEff_by_gene";

my $intergenic = 0;
GetOptions(
	'o|outdir:s'   => \$outdir,
	'intergenic!'  => \$intergenic,
	);

my $header = <>;
unless ( -d $outdir ) {
	mkdir($outdir);
}
my %data;
while(<>) {
	my @row = split(/\t/,$_);
	my $gene = $row[5];
	next if ( ! $intergenic && $row[3] eq 'intergenic' );
	push @{$data{$gene}}, $_;
}

for my $gene ( keys %data ) {
	open(my $outfh => ">$outdir/$gene.snpEff.tab") || die $!;
	print $outfh $header;
	for my $r ( @{ $data{$gene} } ) {
		print $outfh $r;
	}
}
