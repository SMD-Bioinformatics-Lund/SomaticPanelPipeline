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
my $MIN_DP = 100;
my %genes_per_assay = (
    'myeloid'=>{'CEBPA'=>1},
    'lymphoid'=>{'CEBPA'=>1},
    'PARP_inhib'=>{'ALL_GENES'=>1},
    'solid'=>{'CEBPA'=>1, 'BRCA1'=>1, 'BRCA2'=>1}
    );


# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'tumor-id=s', 'normal-id=s', 'assay=s', 'assay_json=s' );
check_options( \%opt );

# Default to myeloid genes, for backwards compatibility
my %GENES = %{$genes_per_assay{'myeloid'}};
my %assay_json;
if( $opt{'assay_json'}) {
    my $assay_json = read_json($opt{'assay_json'});
    %assay_json = %$assay_json;
    #print Dumper(%assay_json);
    %GENES = ();
    foreach my $gene ( @{ $assay_json{"genes"} }) {
        $GENES{$gene} = 1;
    }
    # print Dumper(%GENES);

    # exit;
}

elsif( $opt{'assay'} ) {
    my $assay = $opt{'assay'};
    die "No genes set for assay: $assay\n" if !$genes_per_assay{$assay};
    %GENES = %{$genes_per_assay{$assay}};
}

my $vcf = vcf2->new('file'=>$opt{vcf} );

my $nid = ( $opt{'normal-id'} or 0 );
my $tid = ( $opt{'tumor-id'} or 0 );

print_header($opt{vcf});


while ( my $var = $vcf->next_var() ) {

    # Check if variant is in one of the selected genes
    my $in_relevant_gene = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if( $GENES{ $tx->{SYMBOL} } or $GENES{'ALL_GENES'} ) {
            $in_relevant_gene = 1;
            last;
        }
    }

    my $germline = 0;
    my $germline_risk = 0;
    my $not_germline = 0;
    
    # check whether it seems like a germline variant depending on VAF of variant
    for my $gt ( @{$var->{GT}}) {
        # If normal samples. Set GERMLINE filter if VAF > $MIN_VAF_NORMAL and in relevant gene
        if( $nid and $gt->{_sample_id} eq $nid ) {
            $germline = 1 if( $gt->{VAF} > $MIN_VAF_NORMAL and $gt->{DP} > $MIN_DP and $in_relevant_gene);
            $not_germline = 1 if $gt->{VAF} < $MAX_VAF_NORMAL and $gt->{DP} > $MIN_DP;
        }
        
        # If tumor samples. Set GERMLINE_RISK filter if VAF > $MIN_VAF_TUMOR
        if( $tid and $gt->{_sample_id} eq $tid ) {
            $germline_risk = 1 if( $gt->{VAF} > $MIN_VAF_TUMOR and $gt->{DP} > $MIN_DP );
        }
    }

    # if using exclusion criterias from assay_json reset germline if variant fail check
    if( $opt{'assay_json'} && $germline) {
        $germline = mini_rank($var,\%assay_json);
    }

    # Add GERMLINE filter and remove FAIL_NVAF for confirmed relevant germlines
    if ( $germline ) {
        my @new_filters = ("GERMLINE",$germline);
        foreach( split(';', $var->{FILTER}) ) {
            push @new_filters, $_ unless $_ eq "FAIL_NVAF"; 
        }
        $var->{FILTER} = join(";", @new_filters);
    }
    

    # Add GERMLINE_RISK filter to suspected germline variants in tumor only
    elsif ( $germline_risk and !$nid ) {
        my @filters = split ';', $var->{FILTER};
        push @filters, "GERMLINE_RISK";
        $var->{FILTER} = join(";", @filters);
    }

    vcfstr($var);
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
    my $rank = shift;
    my %rank = %$rank;
    my $score = 0;
    my $gnomad = 0;


    ## check clinvar ##
    my $clinvar = 0;
    if ($var->{INFO}->{CLNSIG}) {
        my @clinsig = split("/",$var->{INFO}->{CLNSIG});
        if (scalar(@clinsig) == 1) {
            if ($rank{"clinvar"}{$clinsig[0]}) {
                $clinvar = $rank{"clinvar"}{$clinsig[0]};
            }
        }
        else {
            foreach my $match (@clinsig) {
                if ($rank{"clinvar"}{$match}) {
                    if ($rank{"clinvar"}{$match} > $clinvar) {
                        $clinvar = $rank{"clinvar"}{$match};
                    }                    
                }
            }
        }


    }
    ## check consequence ## ## check gnomad ##
    my $max_score = 0;
    for my $tx ( @{ $var->{INFO}->{CSQ} } ) {
        if( $tx->{Consequence} ) {
            foreach my $csq (@{$tx->{Consequence}} ) {
                if ($rank{"consequence_score"}{$csq}) {
                    if ( $rank{"consequence_score"}{$csq} > $max_score ) {
                        $max_score = $rank{"consequence_score"}{$csq};
                    }
                }
            }
        }
        if ( $tx->{gnomADg}) {
            my @afs =  split(",",$tx->{gnomADg_AF});
            if ($afs[0] < $assay_json{"gnomad_cutoff"}) {
                $gnomad = 1;
            }
        }
    }
    if ($max_score >= $assay_json{"consequence_cutoff"}) {
        $score++;
    }
    if ($gnomad >= 1) {
        $score++;
    }
    if ($clinvar) {
        $score = $clinvar + $score;
    }
    if ($score > 1) {
        return 1;
    }
    else {
        return 0;
    }
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
