#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my %opt = ();
GetOptions( \%opt, 'vcf=s', 'af=s', "id=s" );

my $vcf = vcf2->new('file'=>$opt{vcf} );

my $bnd = $opt{'id'}."_manta_bnd_filtered.vcf";
my $deldup = $opt{'id'}."_manta_filtered.vcf";
open (BND,'>>', $bnd);
open (DELDUP,'>>', $deldup);

my @header = split/\n/,$vcf->{header_str};
my $count = 1;
foreach my $header (@header) {
    
    if ($header =~ /^##INFO/ && $count == 1) {
        $count++;
    }
    if ($header =~ /##contig/) {
        if ($header =~ /ID=[0-9XYM]{1,2},/ ) {
            print BND $header."\n";
            print DELDUP $header."\n";
        }
    }
    else {
        print BND $header."\n";
        print DELDUP $header."\n";
    }

}
while ( my $a = $vcf->next_var() ) {

    my @str = split/\t/,$a->{vcf_str};
    my $start = $a -> {POS};
    my $paired;
    if (scalar(@{$a->{GT}}) > 1 ) {
        $paired = "true";
    }
    ## print BND separately.
    if ($a->{INFO}->{SVTYPE} eq 'BND') {
        my $pass = AF($a,$opt{'id'},$opt{'af'});
        if ($pass) {        
            print BND $a->{vcf_str}."\n";
        }
    }
    else {
        ## Filter manta based upon PR & SR fields in format
        my $pass = AF($a,$opt{'id'},$opt{'af'});
        if ($pass) {
            print DELDUP $a->{vcf_str}."\n";
        }
    }

}


sub AF {
    my ($var,$id,$AF) = @_;

    my $pass = 1;
    foreach my $ind (@{ $var->{GT}}) {
        if ($ind->{_sample_id} eq $id) {
            my @pr = split(',',$ind->{PR} );
            my @sr = split(',',$ind->{SR} );
            if ($pr[0] == 0 || $sr[0] == 0) {
                $pass = 0;
                next;
            }
            #print $var->{CHROM}."\t".$var->{POS}."\n";
            my $af_pr = $pr[1]/$pr[0];
            my $af_sr = $sr[1]/$sr[0];
            #print "$af_sr  $af_pr\n";
            if ($af_sr < $AF and $af_pr < $AF) {
                $pass = 0;
            }
        }
    }
    return $pass;
}