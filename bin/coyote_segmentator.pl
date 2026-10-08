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
my $rand = int rand 1000000000;
my %opt = ();
GetOptions( \%opt, 'vcf=s', 'id=s', 'panel=s', 'normal', 'genes=s');

my %genes_panel;
my @panels = split(',',$opt{'panel'});
## read panels in order separated by comma. The base panel + specific. 
## Overlapping genes will have the information from the last in the list
foreach my $panel (@panels) { 
    read_panel($panel);
}
my $bedtools = read_vcf($opt{'vcf'});

annotate_genes($bedtools, $opt{"genes"},\%genes_panel);
sub read_vcf {
    my $fn = shift;
    my $vcf = vcf2->new('file'=>$fn );
    my $callers = callers($vcf->{header_str});
    my @callers = @$callers;
    my $bedtools = "tmp.bed";
    open (BED, ">>", $bedtools);
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
        ## for gatk-normals
        elsif ( $var->{INFO}->{gatkCN}) {
            $fold = $var->{INFO}->{gatkCN}/2;            
        }
        else {
            if ($type eq "DEL") {
                $fold = "DEL";
            }
            elsif ($type eq "DUP" || $type eq "INS") {
                $fold = "AMP"
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

        ## ADD PR AND SR for manta-variants
        my ($PR,$SR) = (0,0);
        if ($callers =~ /manta/ ) {
            ($PR,$SR) = manta_variant($var,$opt{'id'});
        }
       

        print BED $chrom."\t".$start."\t".$end."\t".$probes."\t".$fold."\t".$type."\t".$callers."\t"."$PR:$SR\n";

        

    }
    return $bedtools;
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

    while(<PANEL>) {
        chomp;
        my @line = split("\t");
        my $gene = $line[0];
        my $type = $line[1];
        my $question = $line[2];
        $genes_panel{$gene}{TYPE} = $type;
        $genes_panel{$gene}{QUESTION} = $question;
    }

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

sub annotate_genes {
    my ($bed, $genes, $genes_panel)  = @_;
    my $rnd = int rand 1000000000;
    my $tmp_infile = "input.$rnd.bed";
    my @overlap = `bedtools intersect -a $bed -b $genes -loj`;
    unlink ("tmp.bed");
    my $prev_reg;
    my @genes;
    my $prev_bed_str;

    my $filtered_bed = $opt{'id'}.".cn-segments.panel.bed";
    my $unfiltered_bed = $opt{'id'}.".cn-segments.bed";
    unlink( $filtered_bed );
    unlink( $unfiltered_bed );
    system(" touch $filtered_bed ");
    system(" touch $unfiltered_bed ");
    open (FILTERED, '>>', $filtered_bed);
    open (UNFILTERED, '>>', $unfiltered_bed);
            
    foreach my $line (@overlap) {
        chomp $line;
        my @f = split /\t/, $line;

        my $reg = "$f[0]\t$f[1]\t$f[2]\t";
        my $bed_str = "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]";
        if( $prev_reg and $reg ne $prev_reg ) {
            my $gene_conc = join(',',@genes);
            my ($match,$intresting) = get_panel_genes($gene_conc,$genes_panel);
            my @bed_str = split("\t",$prev_bed_str);
            if ($match > 0) { 
                print FILTERED join("\t",@bed_str[0..5]);
                print FILTERED "\t".$gene_conc."\t".$intresting."\t";
                print FILTERED join("\t",@bed_str[6..7]);
                if ($opt{'normal'}) {
                    print FILTERED "\tNORMAL";
                }
                print FILTERED "\n";
            }
            print UNFILTERED join("\t",@bed_str[0..5]);
            print UNFILTERED "\t".$gene_conc."\t";
            print UNFILTERED join("\t",@bed_str[6..7]);
            print UNFILTERED "\n";
            @genes = ();
        }
        push @genes, $f[-1];
        
        $prev_reg = $reg;
        $prev_bed_str = $bed_str;
    }
    ### PRINT LAST VARIANT, ugly code try to reduce when you have time
    my $gene_conc = join(',',@genes);
    my ($match,$intresting) = get_panel_genes($gene_conc,$genes_panel);
    my @bed_str = split("\t",$prev_bed_str);
    if ($match > 0) {
        
        print FILTERED join("\t",@bed_str[0..5]);
        print FILTERED "\t".$gene_conc."\t".$intresting."\t";
        print FILTERED join("\t",@bed_str[6..7]);
        if ($opt{'normal'}) {
            print FILTERED "\tNORMAL";
        }
        print FILTERED "\n";
    }
    print UNFILTERED join("\t",@bed_str[0..5]);
    print UNFILTERED "\t".$gene_conc."\t";
    print UNFILTERED join("\t",@bed_str[6..7]);
    print UNFILTERED "\n";
}