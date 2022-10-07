#!/usr/bin/perl -w
use strict;


my $cr = $ARGV[0];
my $baf = $ARGV[1];
my $purity = $ARGV[2];
open (CR, $cr) or die $!;

while (<CR>) {
    next if ($_ =~ /^@/ or $_ =~ /^C/ );
    chomp;
    my @tmp = split("\t",$_);

    my $foldchange = (2**$tmp[4])/$purity;

    print $tmp[0]."\t".$tmp[4]."\t".$foldchange."\t".$tmp[5]."\n";
}

# /mnt/beegfs/nextflow/21KF00523-1.gmssolidtumorv3-0-hrd/67/8c3916fa6d5667db109cc1abe6d0c4