#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $in = Bio::SeqIO->new(-format => "fasta", -file => shift || die $!);
while(my $seq = $in->next_seq ) {
	my $str = $seq->seq;
	my $gap = $str =~ tr/\.-/.-/;
	print join("\t",$seq->display_id, $gap,sprintf("%.2f",$gap/$seq->length)),"\n";
}
