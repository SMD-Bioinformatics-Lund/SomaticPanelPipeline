#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;


# Open VCF
my $vcf = vcf2->new('file'=>$ARGV[0] );


while ( my $var = $vcf->next_var() ) {

    my $test;
    foreach my $key (keys %{ $var->{INFO}->{CSQ}->[0]}) {
        if ($key =~ /\w+_OID/) {
            
            if ($var->{INFO}->{CSQ}->[0]->{$key} ne "") {
                print $key."\t";
                print $var->{INFO}->{CSQ}->[0]->{$key}."\n";
            }
            
        }
    }


}