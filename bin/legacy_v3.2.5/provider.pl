#!/usr/bin/env perl
use strict;
use threads;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use lib dirname (__FILE__);

my $SAMTOOLS_PATH = "samtools";
my $DEFAULT_SNPBED = dirname($0)."/HPA_1000G_final_38.bed";
my $DEFAULT_XYBED  = dirname($0)."/xy_38.bed";

my %opt = &get_options;
my %meta_data = &get_bampaths( $opt{bam}, $opt{nocheck} );
my ($variant_data, $chromosomes, $allele_freq)  = &read_bed( $opt{bed} );

my ($xy_data, $xy_chromosomes);
if ( $opt{ bedxy } ) {
    ($xy_data, $xy_chromosomes) = &read_bed( $opt{bedxy} );
    foreach my $chr (keys %{$xy_chromosomes}) {	$chromosomes->{$chr}++; }
}

unless ( $opt{ nocheck } or &matching_chr_names( (keys %meta_data)[0], $opt{bed} ) ) {
    error( "Chromosome names do not match between main BED and BAMs", 1 );
}

unless ( $opt{ nocheck } or &matching_chr_names( (keys %meta_data)[0], $opt{bedxy} ) ) {
    error( "Chromosome names do not match between xy BED and BAMs", 1 );
}

my %data = &get_base_freqs_from_bams( \%meta_data, $opt{bed}, $opt{bedxy}, $chromosomes );

my %sample = &determine_sex( \%data, $xy_data );
&do_genotyping( \%data, $variant_data, $allele_freq );

&print_genotype_table( \%data, \%sample, $opt{out}, \%meta_data, $variant_data );

&detect_unexpected( \%data, \%meta_data, \%sample, $allele_freq );



########################


sub detect_unexpected{
    my( $data, $meta_data, $sample_data, $allele_freq ) = @_;

    my (%dist, %seen);
    foreach my $samp1 ( keys %$meta_data ) {
	my $name1 = $meta_data->{$samp1}->{name};
	$seen{$samp1} = 1;

	# If sex was annotated, check if predicted sex agrees.
	if ($meta_data->{$samp1}->{sex} and $sample_data->{$name1}->{sex}->{pred}) {
	    my $anno_sex = $meta_data->{$samp1}->{sex};
	    my $pred_sex = $sample_data->{$name1}->{sex}->{pred};
	    print "UNEXPECTED SEX: $name1 (Predicted:$pred_sex, Annotated:$anno_sex)\n" if ($anno_sex ne "-" and $anno_sex ne $pred_sex);
	}
	

	# Check sample simlarity based on genotypes.
	foreach my $samp2 (keys %$meta_data) {
	    next if $seen{$samp2};

	    my $name2 = $meta_data->{$samp2}->{name};
	    unless (defined( $dist{$name2}->{$name1} )) {
		my $dist = distance( $data, $name1, $name2 );

		next if $dist == -1; # Skip samples where distance could not be calculated (all NAs)

		$dist{$name1}->{$name2} = $dist;
		
		# If individual ID exist for both samples, check if annotation & prediction agrees.
		if ( $meta_data->{$samp1}->{individual} and $meta_data->{$samp2}->{individual} ) {
		    my ($ind1, $ind2) = ( $meta_data->{$samp1}->{individual}, $meta_data->{$samp2}->{individual} );

		    if ($ind1 eq $ind2 and $dist > 0.05) {
			printf "UNEXPECTED DIFFERENT: %s - %s (%.2f%%)\n", $name1, $name2, 100*(1-$dist{$name1}->{$name2});
		    }
		    elsif ($ind1 ne $ind2 and $dist < 0.05) {
			printf "UNEXPECTED IDENTICAL: %s - %s (%.2f%%)\n", $name1, $name2, 100*(1-$dist{$name1}->{$name2});
		    }
		}

		# If individual ID is missing for either sample, check if samples appaer to be from same individual
		else {
		    if ($dist < 0.05) {
			printf "Samples %s and %s appear to be same individual/twins (%.2f%%)\n", 
			       $name1, $name2,100*(1-$dist{$name1}->{$name2});
		    }
		}
	    }
	}
    }
}



sub distance{
    my( $data, $id1, $id2 ) = @_;

    my ( $identical, $tot ) = ( 0, 0 );
    foreach my $loc (keys %$data) {
	if( defined($data->{$loc}->{samples}->{$id1}->{basecall}) and defined($data->{$loc}->{samples}->{$id2}->{basecall}) ) {
	    $identical++ if $data->{$loc}->{samples}->{$id1}->{basecall} eq $data->{$loc}->{samples}->{$id2}->{basecall};
	    $tot++;
	}
    }

    if ($tot > 0) {
	return ($tot-$identical) / $tot;
    }
    return -1;
}



sub read_bed{
    my $bed = shift;
    open( BED, $bed );
    my (%var, %chr, %additional);
    while( <BED> ) {
	chomp;
	my @a = split /\t/;
	$chr{$a[0]} ++;
	if ($a[3]) {
	    $var{$a[0].":".$a[2]} = $a[3];
	}
	else {
	    $var{$a[0].":".$a[2]} = $a[0].":".$a[2];
	}

	if ($a[4]) {
	    $additional{$a[0].":".$a[2]} = $a[4];
	}
    }

    return( \%var, \%chr, \%additional );
}



sub print_genotype_table{
    my( $snp_data, $sample_data, $out, $annotation, $variants ) = @_;

    open( GT, ">".$out.".genotypes" );

    if( ! $opt{'long'} ) {
	print GT "loc";
	print GT "\t$_" foreach sort keys %{ ((values %$snp_data)[0])->{samples} };
	print GT "\n";
    }
  
    my ($tot_callable, $tot_sites, $HW_pass, $tot_loc);
    foreach my $loc ( sort keys %$snp_data ) {
	my $snp_id = $opt{position} ? $loc : $variants->{$loc};

	print GT $snp_id if ! $opt{long};

	foreach my $sid ( sort keys %{ $snp_data->{$loc}->{samples} } ) {
	    my $genotype = $snp_data->{$loc}->{samples}->{$sid}->{basecall};

	    if( $opt{long} ) {
		print GT "$sid\t$sid\t".$snp_id."\t". 
		         plink_gt( $snp_data->{$loc}->{samples}->{$sid}->{basecall} ) ."\n";
	    }
	    else {
		print GT "\t". ( defined( $genotype ) ? $genotype : "NA" ) if ! $opt{long};
	    }
	}

	print GT "\n" if ! $opt{long};
	$tot_callable += $snp_data->{$loc}->{Ncallable};
	$tot_sites    += $snp_data->{$loc}->{Ntotal};
	$tot_loc      ++;
    }   
    close GT;

    open (STATS, ">".$out.".stats");
    printf STATS "Average callability: %.2f%%\n", 100*($tot_callable/$tot_sites);
    close STATS;

    # Output file with sex prediction
    open (SAMPLES, ">".$out.".sex");
    print SAMPLES "sample\tpredicted_sex\tannotated_sex\n";
    foreach my $bam ( sort keys %$annotation ) {
	my $sid = $annotation->{$bam}->{name};
	print SAMPLES $sid."\t". ($sample_data->{$sid}->{sex}->{pred} or "-")."\t".($annotation->{$bam}->{sex} or "-")."\n";
    }
    close SAMPLES;
}



sub plink_gt{
    my $bc = shift;
    return "0\t0" if !defined($bc) or $bc eq "unclear" or $bc eq "lowdata" ;
    return "$bc\t$bc" if length($bc) == 1;
    my ($a1, $a2) = split /\//, $bc;
    return "$a1\t$a2";
}

 

sub determine_sex{
    my( $data, $loci ) = @_;

    # Determine average coverage of normal loci, for comparison with sex determination loci
    my %avg_cov;
    foreach my $loc (keys %$data) {
 	next if $loci->{$loc};
 	foreach my $sid (keys %{ $data->{$loc}->{samples} }) {
 	    push( @{ $avg_cov{$sid} }, $data->{$loc}->{samples}->{$sid}->{depth} );
 	}
    }
    $avg_cov{$_} = mean( @{$avg_cov{$_}} ) foreach keys %avg_cov;

    # Check whether sex loci show high/low coverage compared to normal loci 
    my %sample;
    my( %F, %M );
    foreach my $loc ( keys %$loci ) {
	foreach my $sid (sort keys %{ $data->{$loc}->{samples} }) {
	    my $sex_depth = $data->{$loc}->{samples}->{$sid}->{depth};
	    
	    # FIXME: Arbitrary numbers...
	    if ($sex_depth / $avg_cov{$sid} >= 0.2) { 
		$M{$sid}++;
		$sample{$sid}->{sex}->{loc}->{$loc} = "M";
	    }
	    elsif ($sex_depth / $avg_cov{$sid} <= 0.05) {
		$sample{$sid}->{sex}->{loc}->{$loc} = "F";
		$F{$sid}++;
	    }
	    else {
		$sample{$sid}->{sex}->{loc}->{$loc} = "unclear";
	    }
	}

	# If all loci agree -> predicted sex.
 	foreach my $sid (sort keys %{ $data->{$loc}->{samples} }) {
	    if( $M{$sid} and !$F{$sid} ) {
		$sample{$sid}->{sex}->{pred} = "M";
	    }
	    elsif( $F{$sid} and !$M{$sid} ) {
		$sample{$sid}->{sex}->{pred} = "F";
	    }
	    else {
		$sample{$sid}->{sex}->{pred} = "unclear";
	    }
	}

	delete( $data->{$loc} );
    }
    return %sample;
}



# Run mpileup for a single chromosome in a thread.
sub mpileup_thread {
    my ($chr, $out, $bams) = @_;
    my $id = threads->tid();
    system( $SAMTOOLS_PATH." mpileup -l ". $opt{out}.".tmp.bed -r $chr $bams".
	    " > $out.$chr 2> $opt{out}.mpileup.$chr.log" );
    threads->exit();
}



# Run mpilup on bams and count number of A, C, G and Ts in each position.
sub get_base_freqs_from_bams{
    my( $file_data, $snp_fn, $xy_fn, $chromosomes ) = @_;

    # Define final pileup output file name
    my $pile_out = $opt{out}.".pile.tmp";

    # Get BAM files into an array
    my @bams = sort keys %$file_data;

    # Merge the two BED files
    merge_files( [$snp_fn, $xy_fn], $opt{out}.".tmp.bed" );

    if( $opt{overwrite} or !-s $pile_out ) {

	# Run mpileup in threads, split by chromosomes
	if( $opt{thread} ) {
	    my @threads;
	    my $i = 0;
	    my @pileup_files;
	    foreach my $chr ( sort { $chromosomes->{$b} <=> $chromosomes->{$a} } keys %{$chromosomes} ) {
		$threads[$i] = threads->create( \&mpileup_thread, ( $chr, $pile_out, join( " ", @bams ) ) );
		push( @pileup_files, "$pile_out.$chr" );
		$i++;
	    }
	    $_->join() for @threads;

	    merge_files( \@pileup_files, $pile_out );
	    unlink $_ foreach @pileup_files;
	}

	# Run mpileup unthreaded
	else {
	    system( $SAMTOOLS_PATH." mpileup -l ". $opt{out}.".tmp.bed ". join( " ", @bams ).
		    " > $pile_out 2> $opt{out}.mpileup.log" );
	}

    }
    else {
	print STDERR "WARNING: Prior pileup file found, resuming with that file. Use --overwrite to override\n";
    }

    # Collect base frequencies from pileup file
    my %frq;
    open( PILE, $pile_out );
    while (<PILE>) {
	my @p = split /\t/;
	my $loc = $p[0].":".$p[1];
	
	my $i;
	foreach my $fn ( @bams ) {
	    $i += 3;
	    my $sample_name = $file_data->{$fn}->{name};
	    my ($depth, $base, $qual) = ( $p[$i], $p[$i+1], $p[$i+2] );
	    my $bases = count_bases( $base, $qual );
	    $frq{$loc}->{samples}->{ $sample_name }->{bases} = $bases;
	    $frq{$loc}->{samples}->{ $sample_name }->{depth} = sum( values %{ $bases } );
	}
    }
    close PILE;
    return %frq
}



sub do_genotyping{
    my( $data, $loci, $allele_freq ) = @_;

    foreach my $loc ( keys %$loci ) {
	my ($num_samples, $num_callable) = (0, 0);

	# SNP loci: Do base calling
	my %gt_cnt;
	foreach my $sid (sort keys %{ $data->{$loc}->{samples} }) {
    
	    # Get base frequencies for the locus into %f 
	    my %f = %{ $data->{$loc}->{samples}->{$sid}->{bases} };

	    # Get read depth
	    my $tot = $data->{$loc}->{samples}->{$sid}->{depth};

	    # Find most common and second most common base
	    my ($max, $second) = ( large(\%f, 1), large(\%f, 2) );

	    # Require a depth of at least 6 reads in the locus
	    my $gt;
	    if ($tot >= 6) {
		my $ap = ( $f{$max} - $f{$second} ) / $tot;
	
		if ($ap < 0.6) { # Heterozygote
		    $gt = join "/", sort( $max,$second );
		}
		elsif ($ap > 0.9) { # Homozygote
		    $gt = $max;
		}
		else {
		    $gt = 'unclear';
		}
	    }
	    else {
		$gt = 'lowdata';
	    }
	    if ($gt ne "unclear" and $gt ne "lowdata") {
		$data->{$loc}->{samples}->{$sid}->{basecall} = $gt;
		$gt_cnt{$gt}++;
		$num_callable++;
	    }
	    else {
		$data->{$loc}->{samples}->{$sid}->{uncalled} = $gt;
	    }
	    $num_samples++;
	}
	$data->{$loc}->{Ntotal}    = $num_samples;
	$data->{$loc}->{Ncallable} = $num_callable;


	# Skipping tri-allelic loci
	if( keys %gt_cnt > 3 ) {
	    print STDERR "Skipping $loc: More than more than 2 alleles detected!\n";
	    $data->{$loc}->{skip} = 1;
	    next;
	}

	# Skip loci that are not callable in any sample
	if( $num_callable < 1 ) {
	    print STDERR "Skipping $loc: No samples callable!\n";
	    $data->{$loc}->{skip} = 1;
	    next;
	}

	# If alt allele frequency was specified and --obsaf option was not used,
	# simply add it to the locus hash.
	if( $allele_freq->{ $loc } and !$opt{ obsaf } ) {
	    my( $alt, $alt_af ) = split /\//, $allele_freq->{ $loc };
	    $data->{$loc}->{alt} = $alt;
	    $data->{$loc}->{alt_af} = $alt_af;
	}
	# Otherwise calculate an observed alt allele frequency
	else {
	    my $first = 1;
	    my ($rare_hom, $het) = (0, 0);
	    foreach my $gtype (sort {$gt_cnt{$a} <=> $gt_cnt{$b}} keys %gt_cnt) {
		if ($gtype =~ /\//) {
		    $het = $gtype;
		}
		else {
		    $rare_hom = $gtype if $first;
		    $first = 0;
		}
	    }

	    # If all samples are heterzygotic, assign one of the alleles as rare.
	    $rare_hom = (split /\//, $het)[0] if !$rare_hom;

	    my ($AB, $BB) = ( ($gt_cnt{$het} or 0), ($gt_cnt{$rare_hom} or 0) );
	    $data->{$loc}->{alt} = $rare_hom;
	    $data->{$loc}->{alt_af} = ( $AB + 2 * $BB ) / ( 2 * $num_callable );
	}
    }
}



# Merge a number of files, similar to 'cat file1 file2 > merged'
sub merge_files{
    my ($files, $out) = @_;
    
    open( OUT, ">$out" );
    foreach (@$files) {
	open (F, $_);
	my @a = <F>;
	close F;
	print OUT join( "", @a );
    }
    close OUT;
}



# FIXME: Use qual string and parse base string properly...
sub count_bases{
    my ($b, $q) = @_;
    if ($b) {
	my $a = ($b =~ s/A/A/gi);
	my $t = ($b =~ s/T/T/gi);
	my $g = ($b =~ s/G/G/gi);
	my $c = ($b =~ s/C/C/gi);
	return {'A'=>$a, 'C'=>$c, 'G'=>$g, 'T'=>$t};
    }
    else {
	return {'A'=>0, 'C'=>0, 'G'=>0, 'T'=>0};
    }
}



# Parse and check command line options
sub get_options{
    my %opt = ( 'bed' => $DEFAULT_SNPBED, 'bedxy' => $DEFAULT_XYBED );
    GetOptions( \%opt, 'bam=s', 'bed=s', 'bedxy=s', 'nocheck', 'overwrite', 'out=s', 
		'thread', 'long', 'position', 'obsaf' );

    error( "Parameter --bam required", 1, 1 ) unless $opt{bam};
    error( "Bed file $opt{bed} does not exist.", 1 ) if $opt{bed} and !-s $opt{bed};

    unless ($opt{bedxy} eq "none") {
	error( "Sex determination bed file $opt{bedxy} does not exist.", 1 ) if $opt{bedxy} and  !-s $opt{bedxy};
    }
    else {
	$opt{bedxy} = "";
    }

    error( "Parameter --out required.", 1, 1 ) unless $opt{out};

    if (!$opt{overwrite} and ( -s $opt{out}.".genotypes" or -s $opt{out}.".stats" ) ) {
	error( "Output with prefix $opt{out} already exists. Use --overwrite to overwrite.", 1 );
    }
    return %opt;
}



# Display help text
sub display_usage{
    print "provider.pl --bam [METADATA FILE|'FILEMASK'] --out OUT_FILE_PREFIX\n";

    print " REQUIRED\n".
	  "   --bam        Either the path to a metadata file listing bam files or a\n".
	  "                filemask for bam files\n".
	  "   --out        Prefix of output files\n".
	  " OPTIONAL:\n".
	  "   --bed        Path to BED file of SNPs to use for finger printing\n".
	  "                Default: $DEFAULT_SNPBED\n".
	  "   --overwrite  Overwrite any existing files with same file prefix\n".
          "                Default: OFF\n".
	  "   --nocheck    Don't check if files exist or if chromosomes match.\n".
	  "                Default: OFF\n".
	  "   --long       Output genotypes in long format: sampleID[TAB]snpID[TAB]genotype\n".
	  "                Default: OFF\n".
	  "   --position   Output genomic position as output SNP ID instead of BED-defined\n".
	  "                Default: OFF\n".
	  "   --bedxy      Extra bed file for sex determination. Set to 'none' if not desired.\n".
	  "                Default: $DEFAULT_XYBED\n".
	  "   --thread     Run in threaded mode.\n";
	  "   --obsaf      Use allele frequencies from data for statistics.\n".
          "                Default: OFF if BED files contains allle frequencies, otherwise ON\n\n";

}



# Get BAM file names from file mask or metadata file.
sub get_bampaths{
    my ($mask, $nocheck) = @_;

    my (@files) = sort glob $mask;
    
    # No file, print error
    if (@files == 0) {
	error( "No file(s) found: $mask", 1 );	
    }

    # Single file
    elsif (@files == 1) {
	my $fn = $files[0];
	if (-s $fn or $nocheck) {
	    # If bam/sam file, just return file name
	    return ( $fn=>{ 'name'=>$fn } ) if $fn =~ /\.[bs]am$/;

	    # Otherwise assume metadata table. Read and return data.
	    my %files = read_sample_metadata( $fn, $nocheck );
	    return %files;
	}
	else {
	    error( "No file(s) found: $mask", 1 );
	}
    }

    # Multiple files
    else {
	my %files;
	foreach (@files) {
	    $files{$_}->{ name } = $_;
	}
	return %files;
    }
}



# Parse meta data file
sub read_sample_metadata{
    my($fn, $nocheck ) = @_;
    
    my( %files, %used_names );
    my $cnt;
    open( TABLE, $fn );
    while( <TABLE> ) {
	chomp;

	next if /^#/ or /^\s*$/;  # Skip comments and empty lines

	my @dat = split /\t/;
	if( -s $dat[0] or $nocheck ) {
	    my $name = ($dat[1] or "unknown".++$cnt);
	    if( $used_names{$name} ) {
		error( "ERROR: Sample IDs must be unique ($name found multple times)", 1 );
	    }
	    $used_names{$name} = 1;
	    $files{ $dat[0] }->{ name } = $name;
	    $files{ $dat[0] }->{ individual } = ( $dat[2] or "");
	    $files{ $dat[0] }->{ sex } = $dat[3] if $dat[3];
	}
	else {
	    error( "ERROR: File not found '$dat[0]'", 1 );
	}
    }
    return %files;
}



# Check if all the chromosome names in the BED file are present in the BAM
sub matching_chr_names {
    my ($bam, $bed) = @_;
    my (%bam_chr, %bed_chr);

    # Get all chromosome names from BAM the file header
    my @header = `$SAMTOOLS_PATH view -H $bam`;
    foreach (@header) {
	if( /^\@SQ/ ) {
	    my ($name) = ( $_ =~ /\tSN:(.*?)\t/ );
	    $bam_chr{$name} = 1;
	}
    }

    # Get chromosome names from BED file with defined SNPs
    open( BED, $bed );
    while( <BED> ) {
	my @a = split /\t/;
	$bed_chr{$a[0]} = 1;
    }

    # Fail if any of the chromosomes defined BED in are missing in BAM.
    foreach my $chr (keys %bed_chr) {
	return 0 unless $bam_chr{ $chr };
    }

    # All matching
    return 1;
}



# Print error message and quit program
sub error{
    my ($msg, $error_code, $print_usage) = @_;
    print STDERR "*** ERROR: $_[0] ***\n\n";
    &display_usage if $print_usage;
    exit $error_code;
}



# Sum up values of an array
sub sum{
    my $sum;
    $sum += $_ foreach @_;
    return $sum;
}


# Return the key of the Xth highest value of a hash.
sub large {
    my ($hash, $pos) = @_;

    my $cnt;
    foreach my $i (sort {$hash->{$b}<=>$hash->{$a}} keys %$hash) {
	$cnt++;
	return $i if $cnt == $pos;
    }
    die "Large outside of array";
}



# Calculate the average of values in an array
sub mean {
    my $sum;
    if (@_) {
	$sum += $_ foreach @_;
	return $sum/@_;
    }
    else {
	return 0;
    }
}
