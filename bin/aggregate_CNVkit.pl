#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my $tum = $ARGV[0];
my $tumid = $ARGV[1];
my $norm = $ARGV[2];
my $normid = $ARGV[3];

my $divergent_size = 100;
my @info_fields_to_comma_merge = ("FOLD_CHANGE","FOLD_CHANGE_LOG","PROBES");


my $tumvar = read_vcf($tum,$tumid,"yes");
my %tumvar = %$tumvar;
#print Dumper(%tumvar);

my $normvar = read_vcf($norm,$normid,"no");
my %normvar = %$normvar;
#print Dumper(%normvar);
my $comb = merge();
my %comb = %$comb;
#print Dumper($comb{"17_43065910_43130251_DEL"});
foreach my $var ( keys %comb ) {
	print_variant($var);
}
	

sub read_vcf {
	my( $fn, $id, $print ) = @_;
	my %vars;
	my $vcf;
	if( is_gzipped($fn) ) {
		open($vcf, "zcat $fn |" ) or die "Could not open file $fn";
	}
	else {
		open($vcf, $fn ) or die "Could not open file $fn";
	}
	
	my %agg_info;
	my @agg;

	my @header;
	my $varcount = 0;
	while(<$vcf>) {
		if( /^#/ ) {
			# Modify header entries to allow multiple values for modified fields
			if( /<ID=(.*?),/ ) {
				if( grep /^$1$/, @info_fields_to_comma_merge ) {
					s/Number=1,/Number=\.,/;
				}
			}
			#push @header, $_;
			if ($print eq "yes" ) { 
				if ( $_=~ /(^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT)/) {
					print $1."\t".$tumid."\t".$normid."\n";
				}
				else {
					print $_; 
				}
				

			}
			next;
		}

	
		$varcount ++;
		chomp;
		my @data = split /\t/;
		my $info = $data[7];
		my @info = split /;/, $info;
		my %info;
		foreach ( @info ) {
			my( $key, $val ) = split /=/;
			$info{$key} = ( defined($val) ? $val : "defined" );
		}
		my $varid = $data[0]."_".$data[1]."_".$info{END}."_".$info{SVTYPE};
		my @gt = $data[9];
		$vars{$varid}->{CHROM} = $data[0];
		$vars{$varid}->{pos} = $data[1];
		$vars{$varid}->{id} = $data[2];
		$vars{$varid}->{REF} = $data[3];
		$vars{$varid}->{ALT} = $data[4];
		$vars{$varid}->{QUAL} = $data[5];
		$vars{$varid}->{FILTER} = $data[6];
		$vars{$varid}->{INFO} = \%info;
		$vars{$varid}->{FORMAT} = $data[8];
		$vars{$varid}->{GT}{$id} = $data[9];

	}
	return \%vars;

}





sub merge {

	my %combvars;
	foreach my $Tvar (keys %tumvar) {
		my $check = 0;
		foreach my $Nvar (keys %normvar) {
			if ($Tvar eq $Nvar) { ## same variant
				$check = 1;
				$tumvar{$Tvar}{GT}{$normid} = $normvar{$Nvar}{GT}{$normid};
				my $agginfo = aggregate_info($tumvar{$Tvar}{INFO},$normvar{$Nvar}{INFO});
				$tumvar{$Tvar}{INFO} = $agginfo;
				$combvars{$Tvar} = $tumvar{$Tvar};
				last;
				### aggregate information
			}
			elsif (!defined $combvars{$Nvar}) { ## not the same variant add normal to combined, with empty tumor
				$normvar{$Nvar}{GT}{$tumid} = "0/0:0";
				$combvars{$Nvar} = $normvar{$Nvar};

			}
		}
		unless ($check) {
			$tumvar{$Tvar}{GT}{$normid} = "0/0:0";
			$combvars{$Tvar} = $tumvar{$Tvar};
		}

	}
	return \%combvars;
}


sub aggregate_info {
	my ($tum, $norm) = @_;
	my %tum = %$tum;
	my %norm = %$norm;
	my @info_fields;
	my %agg_info;
	
	$agg_info{END} = $tum{END};
	$agg_info{SVLEN} = $tum{SVLEN};
	$agg_info{SVTYPE} = $tum{SVTYPE};
	for (@info_fields_to_comma_merge) { 
		$agg_info{$_} = $tum{$_}.",".$norm{$_};
	}


	for ("IMPRECISE", "SVTYPE", "END", "SVLEN", "FOLD_CHANGE", "FOLD_CHANGE_LOG", "PROBES") {
		if ($_ eq "IMPRECISE") {
			push @info_fields, "IMPRECISE";
		}
		else {
			push @info_fields, "$_=$agg_info{$_}";
		}
		
	}
	my %info_hash;
	foreach my $item (@info_fields) {
		my @split = split/=/,$item;
		
		$info_hash{$split[0]} = $split[1];
	}
	#print Dumper(%info_hash);
	return \%info_hash;
}

sub print_variant {
	my $pos = shift;
	#print $pos."\n";
	my @varorder = ("CHROM","pos","id","REF","ALT","QUAL","FILTER");
	my @infoorder = ("IMPRECISE", "SVTYPE", "END", "SVLEN", "FOLD_CHANGE", "FOLD_CHANGE_LOG", "PROBES");
	
	foreach my $vv (@varorder) {
		print $comb{$pos}{$vv}."\t";
	}
	foreach my $vi (@infoorder) {
		if ($vi eq 'IMPRECISE') {
			print $vi.";";
		}
		elsif ($vi eq 'PROBES') {
			print $vi."=".$comb{$pos}{INFO}{$vi}
		}
		else {
			print $vi."=".$comb{$pos}{INFO}{$vi}.";";
		}
	}
	print "\t".$comb{$pos}{FORMAT}."\t";
	print $comb{$pos}{GT}{$tumid}."\t".$comb{$pos}{GT}{$normid};
	print "\n";

}

sub is_gzipped {
	my $fn = shift;
	
	my $file_str = `file -L $fn`;
	return 1 if $file_str =~ /gzip compressed/;
	return 0;
}
	 