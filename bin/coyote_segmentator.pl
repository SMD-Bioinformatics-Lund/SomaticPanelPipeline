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
GetOptions( \%opt, 'vcf=s', 'id=s', 'panel=s', 'normal');
my $genes = read_panel($opt{'panel'});
my %genes_panel = %$genes;
read_vcf($opt{'vcf'},\%genes_panel);


sub read_vcf {
    my $fn = shift;
    my $genes_panel = shift;
    my $vcf = vcf2->new('file'=>$fn );
    my $callers = callers($vcf->{header_str});
    my @callers = @$callers;
    my $filtered_bed = $opt{'id'}.".cn-segments.panel.bed";
    my $unfiltered_bed = $opt{'id'}.".cn-segments.bed";
    unlink( $filtered_bed );
    unlink( $unfiltered_bed );
    system(" touch $filtered_bed ");
    system(" touch $unfiltered_bed ");
    open (FILTERED, '>>', $filtered_bed);
    open (UNFILTERED, '>>', $unfiltered_bed);
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
        if ($var->{INFO}->{FOLD_CHANGE_LOG}) {
            $fold = $var->{INFO}->{FOLD_CHANGE_LOG};
        }
        else {
            if ($type eq "DEL") {
                $fold = -0.5;
            }
            elsif ($type eq "DUP" || $type eq "INS") {
                $fold = 1.5
            }
        }
        ## callers ##
        if ( $var->{INFO}->{set} eq "Intersection") {
            $callers = join('&',@callers);
        }
        else {
            $callers = $var->{INFO}->{set};
        }
        ## AMP OR DEL
        if ($type eq "DEL") {
            $type = "-";
        }
        elsif ($type eq "DUP" || $type eq "INS") {
            $type = "+";
        }
        ## genes ##
        my $genes = genes($var->{INFO}->{CSQ});
        my ($match,$intresting) = get_panel_genes($genes,$genes_panel);
        
        ## ADD PR AND SR for manta-variants
        my ($PR,$SR) = (0,0);
        if ($callers =~ /manta/ ) {
            ($PR,$SR) = manta_variant($var,$opt{'id'});
        }
        ## print panel-filtered bed
        if ($match > 0) {
            print FILTERED $chrom."\t".$start."\t".$end."\t".$probes."\t".$fold."\t".$type."\t".$genes."\t".$intresting."\t".$callers."\t"."$PR:$SR";
            if ($opt{'normal'}) {
                print FILTERED "\tNORMAL";
            }
            print FILTERED "\n";
        }
        ## print unfiltered bed
        print UNFILTERED $chrom."\t".$start."\t".$end."\t".$probes."\t".$fold."\t".$type."\t".$genes."\t".$callers."\t"."$PR:$SR";
        if ($opt{'normal'}) {
            print UNFILTERED "\tNORMAL";
        }
        print UNFILTERED "\n";

    }
}

sub genes {
    my $csq = shift;
    my %genes;
    foreach my $transcript ( @{ $csq }) {
        $genes{$transcript->{SYMBOL}} = 1;
    }
    my @genes;
    foreach my $gene (keys %genes) {
        unless ($gene eq "") {
            push @genes,$gene;
        }
        
    }
    my $genes = join(',',@genes);
    return $genes;
}

sub read_panel {
    my $fn = shift;
    open(PANEL, $fn) or die $!;
    my %genes_panel;
    while(<PANEL>) {
        chomp;
        my @line = split("\t");
        my $gene = $line[0];
        my $type = $line[1];
        my $question = $line[2];
        $genes_panel{$gene}{TYPE} = $type;
        $genes_panel{$gene}{QUESTION} = $question;
    }
    return \%genes_panel;
}

sub get_panel_genes {
    my ($genes,$genes_panel) = @_;
    my @genes = split(',',$genes);
    my @intresting;
    foreach my $gene ( @genes ) {
        if ( $genes_panel{$gene} ) {
            push @intresting,$gene.":".$genes_panel{$gene}{TYPE}.":".$genes_panel{$gene}{QUESTION};
        }
    }
    my $any_interest = scalar(@intresting);
    my $intresting = join(',',@intresting);   
    return $any_interest,$intresting;
}

sub callers {
    my $header = shift;
    my @callers;
    if ($header =~ /gatk/ ) {
        push @callers,"gatk";
    }
    if ($header =~ /manta/ ) {
        push @callers,"manta";
    }
    if ($header =~ /cnvkit/ ) {
        push @callers,"cnvkit";
    }
    return \@callers;
}

sub manta_variant {
    my ($format, $id) = @_;

    my $PR; my $SR;
    foreach my $ind ( @{ $format->{GT} }) {
        if ($ind->{_sample_id} eq $id ) {
            $PR = $ind->{PR};
            $PR =~ s/,/\//;
            $SR = $ind->{SR};
            $SR =~ s/,/\//;
        }
    }
    return $PR,$SR;
}