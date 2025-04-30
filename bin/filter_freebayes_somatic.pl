#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );

#   Usage : filter_freebayes_somatic.pl tumor.vcf normal.vcf
#   
#   Filters the variants and taggs variants based on differnt filter criteria
# 
#   Dependencies: This script requires vcf2 module. please make sure the module is in the your bin folder. 

my $vcf = vcf2->new('file'=>$ARGV[0] );

my $T = $ARGV[1];
my $N = $ARGV[2];

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
print '##INFO=<ID=SSC,Number=1,Type=Float,Description="Somatic score">'."\n";
print '##FILTER=<ID=FAIL_GT,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=LOH,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=FAIL_QUAL,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=FAIL_LOD,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=WARN_NOVAR,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=FAIL_NVAF,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=WARN_MQ,Description="Record fails the filters">'."\n";
print '##FILTER=<ID=.,Description="Record fails the filters">'."\n";
system("zgrep ^#CHROM $ARGV[0]");

my $SSC_THRES = 30; # Somatic score threshold
my $LOD_THRES = 10; # Likelihood ratio threshold
my $MIN_VAF_RATIO = 3;

while ( my $v = $vcf->next_var() ) {

    my @filters;
    my @status;

    my( %likelihood, %gl_idx, %genotype, %altobs, %depth );
    my $status = "PASS";
    for my $gt (@{$v->{GT}}) {
        my $type = "T";
        $type = "N" if $gt->{_sample_id} eq $N;

	# Fail if GT is 0/0 for tumor
	$status = "FAIL_GT" if $type eq "T" and $gt->{GT} eq "0/0";
	
	my @GL = split /,/, $gt->{GL};
	my @GT = split /\//, $gt->{GT};
	my @AO = split /,/, ($gt->{AO} or "0");

	my $DP = $gt->{DP};
	my $RO = $gt->{RO};

	# This gets the position of the genotype's likelihood value in the GL array
	my $GL_IDX = ($GT[1]*($GT[1]+1)/2) + $GT[0];

	$depth{$type}      = $DP;
	$altobs{$type}     = \@AO;
	$genotype{$type}   = \@GT;
	$likelihood{$type} = \@GL;
	$gl_idx{$type} = $GL_IDX;
    }


    my $LOD_NORM  = $likelihood{N}->[$gl_idx{N}] - $likelihood{N}->[$gl_idx{T}];
    my $LOD_TUMOR = $likelihood{T}->[$gl_idx{T}] - $likelihood{T}->[$gl_idx{N}];
    my $DQUAL = $LOD_TUMOR+$LOD_NORM;

    # Loss of heterozygosity
    if( $genotype{T}->[0] eq $genotype{T}->[1] && ( $genotype{T}->[1] eq $genotype{N}->[0] or $genotype{T}->[1] eq $genotype{N}->[1] ) ) {
	$status = "LOH";
    }

    # Fail if low somatic score
    $status = "FAIL_QUAL" if $DQUAL<$SSC_THRES;

    # Fail if low likelyhood ratio hor either sample
    $status = "FAIL_LOD" if $LOD_NORM<$LOD_THRES || $LOD_TUMOR<$LOD_THRES;

    my $TALT;
    unless( $genotype{T}->[1] == $genotype{N}->[0] or $genotype{T}->[1] == $genotype{N}->[1] ) {
	    $TALT = $genotype{T}->[1];
    }
    else {
	    $TALT = $genotype{T}->[0];
	    $status = "WARN_NOVAR" if $TALT eq "0";
    }

    # Fail if difference between tumor's and normal's VAF is < 3x.
    if( $depth{N} > 0 and $depth{T} > 0 ) {
	    my $NVAF = $altobs{N}->[$TALT-1] / $depth{N};
	    my $TVAF = $altobs{T}->[$TALT-1] / $depth{T};
	    if( $NVAF > 0 and ($TVAF/$NVAF < $MIN_VAF_RATIO) ) {
	        $status = "FAIL_NVAF";
	    }
    }
    push @status,$status;
    if ($v->{INFO}->{MQM} <= 10 ) {
        push @status, "WARN_MQ";
    }
    $v->{FILTER} .= ";".join(";", @status);
    #$v->{FILTER} = $status;
    add_info( $v, "SSC", $DQUAL );    

    vcfstr($v);
}


#   Function: add_info
#
#   Description: Add information from the vcffile

sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}

#   Function:   vcfstr
#
#   Description: This function prints all the information from the vcf string   
sub vcfstr {
    my $v = shift;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	push @all_info, $info_key."=".$v->{INFO}->{$info_key};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}})."\t";

    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
	my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print join(":", @all_gt)."\t";
    }
    print "\n";
}
