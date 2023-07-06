#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $input = $ARGV[0];
my $input2 = $ARGV[1];
my $id = $ARGV[2];
my $ploidy = $ARGV[3];
if ($ploidy) {
    print "SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy\n";
}
else {
    print "SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\n";
}

my @logr = read_file($input);
my $logr = segments(@logr);
my %logr = %$logr;
my @baf = read_file($input2);
my $baf = segments(@baf);


foreach my $segment (sort { $a <=> $b } keys %logr) {   

    my $cnprint;
    my ($cn,$cn1,$cn2) = find_baf($baf, $logr{$segment});
    if ($cn) {
        $cnprint = "$cn\t$cn1\t$cn2";
    }
    ## if no baf for region, assume cn1 = cn?
    else {
        $cn = 2**$logr{$segment}{VALUE}*2;
        $cn = int($cn + 0.5);
        if ($cn == 0 && $logr{$segment}{CHROM} eq 'Y') {
            next;
        }
        $cnprint = "$cn\t$cn\t0";
    }
    print $id."\t";
    print "chr".$logr{$segment}{CHROM}."\t";
    print $logr{$segment}{START}."\t";
    print $logr{$segment}{END}."\t";
    print $cnprint;
    if ($ploidy) {
        print "\t".$ploidy;
    }
    print "\n";
}





sub find_baf {
    my ($baf,$logr) = @_;
    my %baf = %$baf;
    my %logr = %$logr;
    my $chrom = $logr{CHROM};
    my $start = $logr{START};
    my $end = $logr{END};
    my $logrval = $logr{VALUE};
    #print $logrval."\t";
    my ($cn,$cn1,$cn2);
    foreach my $segment (keys %baf) {

        if ($baf{$segment}{CHROM} eq $chrom ) {
            if ($baf{$segment}{START} >= $start && $baf{$segment}{END} <= $end) {
                
                $cn = 2**$logrval*2;
                $cn2 = $cn * (1 - $baf{$segment}{VALUE});
                $cn1 = $cn * $baf{$segment}{VALUE};
                $cn = int($cn + 0.5);
                $cn1 = int($cn1 + 0.5);
                $cn2 = $cn - $cn1;
                #print $baf{$segment}{VALUE}." $cn $cn1 $cn2";
            }            
        }
        else {
            next;
        }
        return $cn,$cn1,$cn2;
    }
    #print $chrom."\t".$start."\t".$end."\n";
}














sub segments {
    my @array = @_;
    my %segments;
    my $seg = 1;
    for (my $i=0; $i <= scalar(@array); $i++) {
        for (my $j=$i; $j <= scalar(@array); $j++) {
            my @line = split(",",$array[$i]);
            my @line2;
            ## if you reached the end ##
            if ($j == scalar(@array) ) {
                @line2 = ("NA","NA","NA");            
            }
            ## otherwise $j ##
            else {
                @line2 = split(",",$array[$j]);
            }
            ## previous value, for when ends are met ##
            my @last = split(",",$array[$j-1]);
            ## if still on same chrom, continue but shift to next segment, move $i forward to last $j ##
            if ($line[0] eq $line2[0]) {
                ## unless same logR continue segment, continue $i where $j was ##
                unless ( $line[2] == $line2[2]) {
                    $segments{$seg}->{CHROM}=$line[0];
                    $segments{$seg}->{START}=$line[1];
                    $segments{$seg}->{END}=$last[1];
                    $segments{$seg}->{VALUE}=$last[2];
                    $i = $j; 
                    $seg++;
                }
            }
            ## reached new chromosome ##
            else {
                $segments{$seg}->{CHROM}=$line[0];
                $segments{$seg}->{START}=$line[1];
                $segments{$seg}->{END}=$last[1];
                $segments{$seg}->{VALUE}=$last[2];
                $i=$j;
                $seg++;
            }
        }
    }
    #print Dumper(%segments);
    return \%segments;
}




sub read_file {
    my $input = shift;
    open(LOG2R, $input) or die $!;
    my @file;
    while (<LOG2R>) {
        chomp;
        my $seg++;
        my @line = split/\t/,$_;
        my @chrompos = split/_/,$line[0];
        my $chrom = $chrompos[0];
        $chrom =~ s/\"//g;
        my $pos = $chrompos[1];
        $pos =~ s/\"//g;
        $pos." ".$chrom."\n";
        push @file,$chrom.",".$pos.",".$line[1];
    }
    return @file;
}
##         SampleID Chromosome Start_position End_position total_cn A_cn B_cn
## 1 SamplePatient1       chr1          14574       952448        5    0    5
## 2 SamplePatient1       chr1         953394      1259701        3    0    3
## 3 SamplePatient1       chr1        1278085      4551743        2    0    2
## 4 SamplePatient1       chr1        4551885     14124232        2    0    2
## 5 SamplePatient1       chr1       14161231     31062374        3    1    2
## 6 SamplePatient1       chr1       31074785     47428120        4    2    2