#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use JSON;
use Data::Dumper;

my $MAX_VAF_NORMAL = 0.05;
my $MIN_VAF_NORMAL = 0.35;
my $MIN_VAF_TUMOR  = 0.45;
my $MIN_DP         = 100;

# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'normal-id=s' );
check_options( \%opt );

# Default to myeloid genes, for backwards compatibility
my %GENES;

my $vcf = vcf2->new('file'=>$opt{vcf} );

my $nid = ( $opt{'normal-id'} or 0 );
my $tid = ( $opt{'tumor-id'} or 0 );

#print_header($opt{vcf});

my %consequence_score = (
    "frameshift_variant"                => 11,
    "transcript_ablation"               => 10,
    "initiator_codon_variant"           => 9,
    "stop_gained"                       => 8,
    "start_lost"                        => 8,
    "stop_lost"                         => 8,
    "splice_acceptor_variant"           => 8,
    "splice_donor_variant"              => 8,
    "inframe_deletion"                  => 7,
    "transcript_amplification"          => 5,
    "splice_region_variant"             => 5,
    "missense_variant"                  => 5,
    "protein_altering_variant"          => 5,
    "inframe_insertion"                 => 5,
    "incomplete_terminal_codon_variant" => 5,
    "non_coding_transcript_exon_variant"=> 3,
    "synonymous_variant"                => 2,
    "mature_mirna_variant"              => 1,
    "non_coding_transcript_variant"     => 1,
    "regulatory_region_variant"         => 1,
    "upstream_gene_variant"             => 1,
    "regulatory_region_amplification"   => 1,
    "tfbs_amplification"                => 1,
    "5_prime_UTR_variant"               => 1,
    "intron_variant"                    => 1,
    "3_prime_UTR_variant"               => 1,
    "feature_truncation"                => 1,
    "TF_binding_site_variant"           => 1,
    "stop_retained_variant"             => 1,
    "feature_elongation"                => 1,
    "regulatory_region_ablation"        => 1,
    "tfbs_ablation"                     => 1,
    "coding_sequence_variant"           => 1,
    "downstream_gene_variant"           => 1,
    "NMD_transcript_variant"            => 1,
    "intergenic_variant"                => 0
);

my %clinvar_score = (
        "Benign" => -2,
        "Pathogenic" => 1,
        "Likely_pathogenic" => 1
);
my $inclusion_score    = 2;
my $gnomad_cutoff      = 0.01;
my $consequence_cutoff = 5;

while ( my $var = $vcf->next_var() ) {

    my $germline      = 0;
    my $germline_risk = 0;
    my $not_germline  = 0;
    
    # check whether it seems like a germline variant depending on VAF of variant
    for my $gt ( @{$var->{GT}}) {
        # If normal samples. Set GERMLINE filter if VAF > $MIN_VAF_NORMAL and in relevant gene
        if( $nid and $gt->{_sample_id} eq $nid ) {
            $germline = 1 if( $gt->{VAF} > $MIN_VAF_NORMAL and $gt->{DP} > $MIN_DP );
            $not_germline = 1 if $gt->{VAF} < $MAX_VAF_NORMAL and $gt->{DP} > $MIN_DP;
        }
        # Should we really do this? it just creates some false confidence that we can say something is germline in unpaired
        # If tumor samples. Set GERMLINE_RISK filter if VAF > $MIN_VAF_TUMOR
        # if( $tid and $gt->{_sample_id} eq $tid ) {
        #     $germline_risk = 1 if( $gt->{VAF} > $MIN_VAF_TUMOR and $gt->{DP} > $MIN_DP );
        # }
    }

    # try to limit GERMLINE variants
    my $score_results;
    if ($germline) {
        ($germline,$score_results) = mini_rank($var);
    }
    

    # Add GERMLINE filter and remove FAIL_NVAF for confirmed relevant germlines
    if ( $germline ) {
        my @new_filters = ("GERMLINE");
        foreach( split(';', $var->{FILTER}) ) {
            push @new_filters, $_ unless $_ eq "FAIL_NVAF"; 
        }
        $var->{FILTER} = join(";", @new_filters);
        #vcfstr($var);
        #print $var->{CHROM}.":".$var->{POS}
    }
    

    # Add GERMLINE_RISK filter to suspected germline variants in tumor only
    elsif ( $germline_risk and !$nid ) {
        my @filters = split ';', $var->{FILTER};
        push @filters, "GERMLINE_RISK";
        $var->{FILTER} = join(";", @filters);
    }

    #vcfstr($var);
}


sub print_header {
    my $file = shift;

    system("zgrep ^## $file");
    print "##FILTER=<ID=GERMLINE,Description=\"Germline variant, detected in normal sample\">\n";
    print "##FILTER=<ID=GERMLINE_RISK,Description=\"Potential germline variant, from tumor sample\">\n";
    system("zgrep ^#CHROM $file");
        
}


sub read_json {
    my $fn = shift;

    open( JSONFILE, $fn );
    my @json = <JSONFILE>;
    my $decoded = decode_json( join("", @json ) );
    close JSONFILE;

    return $decoded;
}

sub mini_rank {
    my $var = shift;
    my $score = 0;
    my %score_results;
    ## check clinvar ##
    # if only one annotation handle differently, this is to support sole Benign annotations
    # otherwise find highest score acoring to rank->clinvar from JSON. 
    # This let's the user define some value and let other be 0. This is ugly?
    my $clinvar = 0; my $clinvar_effect;
    if ($var->{INFO}->{CLNSIG}) {
        my @clinsig = split("/",$var->{INFO}->{CLNSIG});
        if (scalar(@clinsig) == 1) {
            if ($clinvar_score{$clinsig[0]}) {
                $clinvar = $clinvar_score{$clinsig[0]};
                $clinvar_effect = $clinsig[0];
            }
        }
        else {
            foreach my $match (@clinsig) {
                if ($clinvar_score{$match}) {
                    if ($clinvar_score{$match} > $clinvar) {
                        $clinvar = $clinvar_score{$match};
                        $clinvar_effect = $match;
                    }                    
                }
            }
        }
    }
    # add the clinvar score to score, Benign can disqualify the other two categories completely (is this ok?)
    if ($clinvar) {
        $score = $clinvar + $score;
        $score_results{"clinvar"}{"score"} = $score;
        $score_results{"clinvar"}{"effect"} = $clinvar_effect;
    }
    ## check consequence ## ## check gnomad ##
    # save the most severe consequence defined in JSON
    # check if gnomad passes cutoff from JSON
    my $max_score = 0; my $gnomad = 0; my $max_csq; my $gnomad_af = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        ($max_score,$max_csq) = conseqeunce($tx->{Consequence});
        if ( $tx->{gnomADg_AF} ) {
            my @afs = split(",",$tx->{gnomADg_AF});
            if ($afs[0]) {
                if ($afs[0] < $gnomad_cutoff) {
                    $gnomad = 1;
                    $gnomad_af = $afs[0];
                }
            }
        }
    }
    # add to score if most severe conseqeunce is within JSON cutoff

    if ($max_score >= $consequence_cutoff ) {
        $score++;
        $score_results{"csq"}{"score"} = $max_score;
        $score_results{"csq"}{"effect"} = $max_csq;
    }
    # add to score if gnomad is below cutoff
    if ($gnomad > 0) {
        $score++;
        $score_results{"gnomad"}{"score"} = $gnomad;
        $score_results{"gnomad"}{"af"} = $gnomad_af;
    }
    # if score higher than or equal to inclusion score return true and retain GERMLINE status
    if ($score >= $inclusion_score ) {
        print Dumper(%score_results);
        return 1,\%score_results;
    }
    # else send back 0, and thus GERMLINE is no longer true
    else {
        return 0,\%score_results;
    }
}

sub conseqeunce {
    my $tx = shift;
    my $max_score = 0;
    my $max_effect;
    foreach my $csq (@{$tx} ) {
        if ($consequence_score{$csq}) {
            if ( $consequence_score{$csq} > $max_score ) {
                $max_score = $consequence_score{$csq};
                $max_effect = $csq;
            }
        }
    }
    return $max_score,$max_effect;
}

sub check_options {
    my %opt = %{ $_[0] };

    help_text() unless $opt{vcf};

    die "File does not exist $opt{vcf}..." if ! -s $opt{vcf};

}


sub help_text {
    my $error = shift;
    
    print "\n\$ mark_germlines.pl --vcf INPUT_VCF --tumor-id [--normal-id}\n\n";
    print "   --vcf        Input vcf\n";
    print "   --normal-id  Normal sample ID\n";
    print "   --tumor-id   Tumor sample ID\n";
    print "   --assay      Select an assay. Determines which genes are included\n";
    print "   --assay_json Select an assay according to input json, filter uninteresing GERMLINE for select genes";
    print "\n";
    exit(0);
}


sub vcfstr {
    my( $v ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
    my $key2 = $info_key;
    $key2 = "_CSQ_str" if $info_key eq "CSQ";
    push @all_info, $info_key."=".$v->{INFO}->{$key2};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}});

    
    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
    my @all_gt;
    for my $key ( @{$v->{FORMAT}} ) {
        push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
    }
    print "\t".join(":", @all_gt);
    }
    print "\n";
}
