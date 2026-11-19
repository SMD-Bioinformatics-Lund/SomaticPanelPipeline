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

