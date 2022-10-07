#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );

my $vcf = vcf2->new('file'=>$ARGV[0] );

#TODO ADD INFO HEADER STRINGS PROPERLY!
# also this only works for unpaired so far. FIX for paired or do separate script
# system("zgrep ^## $ARGV[0]");
# system("zgrep ^#CHROM $ARGV[0]");

my @header = split/\n/,$vcf->{header_str};
my $count = 1;
foreach my $header (@header) {    
	if ($header =~ /^##INFO/ && $count == 1) {
		print $header."\n";
        print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">\n";
		$count ++;
	}
    elsif ( $header =~ /ID=VAF/) {

    }
    else {
        print $header."\n";
    }

}

while ( my $v = $vcf->next_var() ) {

    #    my $status = "PASS";
    my @status;
    my (%vaf, %dp, %strand_bias, %msi, %msilen);
    my $VAF;
	my $VD;
    for my $gt (@{$v->{GT}}) {
        $VAF = $gt->{VAF};
		$VD = $gt->{VD};
        #delete @{$v->{GT}}[0]->{VAF};
    }
	my $gnomad = 0;
	if (${$v->{INFO}->{CSQ}}[0]->{"gnomADg_AF"}) {
		$gnomad = ${$v->{INFO}->{CSQ}}[0]->{"gnomADg_AF"};
	}
	
	#print Dumper($v);
	unless ($gnomad >= 0.05 && $VD >= 50) {
		next;
	}

    add_info($v, "AF",$VAF);
    my $prstr = vcfstr($v, []);
    print $prstr;
}





sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}


sub vcfstr {
	my( $v, $sample_order ) = @_;
	
	my @all_info;
	my $tot_str = $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

	# Generate and print INFO field
	for my $info_key (@{$v->{INFO_order}}) {
		if($info_key eq "CSQ") {
           
            #print Dumper($v);
			push @all_info, $info_key."=".$v->{INFO}->{_CSQ_str};
		}
		else {
			push @all_info, $info_key."=".$v->{INFO}->{$info_key};
		}
	}
	$tot_str = $tot_str.join(";", @all_info)."\t";

	# Print FORMAT field
	$tot_str = $tot_str.join(":", @{$v->{FORMAT}})."\t";


	my %order;
	my $i=0;
	if( @$sample_order > 0 ) {
		$order{$_} = $i++ foreach @{$sample_order};
	}
	else {
		$order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
	}
	# Print GT fields for all samples
	for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
		my @all_gt;
		for my $key ( @{$v->{FORMAT}} ) {
			push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
		}
		$tot_str = $tot_str.join(":", @all_gt)."\t";
	}
	$tot_str = $tot_str."\n";
	return $tot_str;
}