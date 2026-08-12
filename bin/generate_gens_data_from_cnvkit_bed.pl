#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname(__FILE__);
use vcf2;

my ($cnvkit_cnr_fn, $vcf_fn, $sample_id, $baf_out, $cov_out);

GetOptions(
    'cnr=s'       => \$cnvkit_cnr_fn,
    'vcf=s'       => \$vcf_fn,
    'sample-id=s' => \$sample_id,
    'baf-out=s'   => \$baf_out,
    'cov-out=s'   => \$cov_out,
) or die "Invalid arguments\n";

die "Missing --cnr\n"       unless $cnvkit_cnr_fn;
die "Missing --vcf\n"       unless $vcf_fn;
die "Missing --sample-id\n" unless $sample_id;
die "Missing --baf-out\n"   unless $baf_out;
die "Missing --cov-out\n"   unless $cov_out;
die "CNR does not exist: $cnvkit_cnr_fn\n" unless -e $cnvkit_cnr_fn;
die "VCF does not exist: $vcf_fn\n"        unless -e $vcf_fn;

# Get normalized log2 depth data from cnvkit output.
my @log2_data;
open(my $cnvkit, $cnvkit_cnr_fn) or die "Cannot open $cnvkit_cnr_fn: $!\n";
my $header = <$cnvkit>;
while (<$cnvkit>) {
    chomp;
    my ($chr, $start, $end, $gene, $depth, $log2, $weight) = split /\t/;
    next if $gene eq "Antitarget";
    my $midpoint = $start + int(($end - $start) / 2);
    push @log2_data, $chr . "\t" . ($midpoint - 1) . "\t" . $midpoint . "\t" . $log2;
}
close $cnvkit;

# Read VCF and calculate BAF data using the legacy vcf2.pm parser behavior.
my @baf_data;
my $vcf = vcf2->new('file' => $vcf_fn);
while (my $var = $vcf->next_var()) {
    for my $gt (@{$var->{GT}}) {
        next unless $gt->{_sample_id} eq $sample_id;
        my $DP = ($gt->{DP} or 0);
        last if $DP < 100 or !$gt->{AO};
        my $AD = $gt->{AO};
        my $vaf = $AD / $gt->{DP};
        push @baf_data, $var->{CHROM} . "\t" . ($var->{POS} - 1) . "\t" . $var->{POS} . "\t" . $vaf;
    }
}

open(my $log2_out, ">", $cov_out) or die "Cannot write $cov_out: $!\n";
for my $resolution (qw(o a b c d)) {
    print $log2_out $resolution . "_" . $_ . "\n" for @log2_data;
}
close $log2_out;

open(my $baf_fh, ">", $baf_out) or die "Cannot write $baf_out: $!\n";
for my $resolution (qw(o a b c d)) {
    print $baf_fh $resolution . "_" . $_ . "\n" for @baf_data;
}
close $baf_fh;
