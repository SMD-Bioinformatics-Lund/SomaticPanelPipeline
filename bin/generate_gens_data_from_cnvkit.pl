#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;

my ($cnvkit_cnr_fn, $vcf_fn, $sample_id) = ($ARGV[0], $ARGV[1], $ARGV[2]);

die unless -e $cnvkit_cnr_fn;
die unless -e $vcf_fn;
die unless $sample_id;

# Get normalized log2 depth data from cnvkit output 
my @log2_data;
open(my $cnvkit, $cnvkit_cnr_fn);
my $header = <$cnvkit>;
while(<$cnvkit>) {
    chomp;
    my( $chr, $start, $end, $gene, $depth, $log2, $weight ) = split /\t/;
    next if $gene eq "Antitarget";
    my $midpoint = $start + int(($end-$start)/2);
    push @log2_data, $chr."\t".($midpoint-1)."\t".$midpoint."\t".$log2;
}
close $cnvkit;

# Read VCF and calculate BAF data
my @baf_data;
my $vcf = vcf2->new('file'=>$vcf_fn);
while ( my $var = $vcf->next_var() ) {
    for my $gt (@{$var->{GT}}) {
	next unless $gt->{_sample_id} eq $sample_id;
	my $DP = ($gt->{DP} or 0);
	last if $DP < 100 or !$gt->{AD};

	my @AD = split /,/, $gt->{AD};
	last unless @AD == 2;

	my $vaf = $AD[1]/$gt->{DP};
	push @baf_data, $var->{CHROM}."\t".($var->{POS}-1)."\t".$var->{POS}."\t".$vaf;
    }
    
}

# Output depth data
open(my $log2_out, ">".$sample_id.".cov.bed");
for my $resolution (qw(o a b c d)) {
    print $log2_out  $resolution."_".$_."\n" for @log2_data;
}
close $log2_out;

# Output BAF data
open(my $baf_out, ">".$sample_id.".baf.bed");
for my $resolution (qw(o a b c d)) {
    print $baf_out  $resolution."_".$_."\n" for @baf_data;
}
close $baf_out;

# Compress and index files
system("bgzip -f $sample_id.cov.bed");
system("tabix $sample_id.cov.bed.gz");
system("bgzip -f $sample_id.baf.bed");
system("tabix $sample_id.baf.bed.gz");
	    

