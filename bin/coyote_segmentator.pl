#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


#####                                                               #####
#                                                                       #
#  Take svdb merged sv-vcf annotated (VEP, maybe db of artefacts).      #
#  Create segments for coyote, perhaps filter really large ones?        #
#                                                                       #
####                                                                #####

# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'normal-id=s');

read_vcf($opt{'vcf'});

sub read_vcf {
    my $fn = shift;

    my $vcf = vcf2->new('file'=>$fn );
    while ( my $var = $vcf->next_var() ) {
        my $chrom = $var->{CHROM};
        my $start = $var->{POS};
        my $type = $var->{INFO}->{SVTYPE};
        next if ($type eq 'BND');
        my $end;
        if ( $var->{INFO}->{END} ) {
            $end = $var->{INFO}->{END};
        }
        else {
            $end = $start + abs($var->{INFO}->{SVLEN});
        }
        ## only in gatk/cnvkit
        my $probes;
        if ($var->{INFO}->{PROBES}) {
            $probes = $var->{INFO}->{PROBES};
        }
        else {
            $probes = 0;
        }
        my $fold;
        if ($var->{INFO}->{PROBES}) {
            $probes = $var->{INFO}->{PROBES};
        }
        else {
            $probes = 0;
        }
        print $chrom."\t".$start."\t".$end."\t".$type."\n";
    }
}