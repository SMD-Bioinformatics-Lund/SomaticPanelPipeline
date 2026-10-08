#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $id = $ARGV[1];
my $ploidy = $ARGV[2];
open(CNVKIT, $input) or die $!;
if ($ploidy) {
    print "SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy\n";
}
else {
    print "SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\n";
}
while (<CNVKIT>) {
    chomp;
    next if ($_ =~ /^chromosome/);
    my @line = split/\t/,$_;
    ## ignore Y chrom with 0 copies, probably female ##
    if ($line[0] eq 'Y' && $line[8] == 0 ) {
        next;
    }
    print $id."\t";
    print "chr".$line[0]."\t";
    print $line[1]."\t";
    print $line[2]."\t";
    print $line[8]."\t";
    if ($line[10] eq '' or $line[9] eq '') {
        my $line = $line[8]-1;
        print $line."\t";
        print "1\t";
    }
    else {
        print $line[10]."\t";
        print $line[9]."\t";
    }
    if ($ploidy) {
        print $ploidy;
    }
    print "\n";
}