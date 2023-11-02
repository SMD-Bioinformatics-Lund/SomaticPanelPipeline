import os
import sys
from pysam import VariantFile
from pprint import pprint
import cmdvcf

"""
data that needs support to be loaded:
SNV/Indels

CNVs

Translocations(dna-fusions)

fusions(rna, maybe hold off for now?)

Sample-information (should be a json)
    id
    group
    build
    gens
    subpanel
    purity
    clarity-id
    clarity-pool

biomarkers (json)
    key:value

cnvprofile
    just path to an image

lowcov
    bed-file

"""


###################
# Insert variants #
###################
bcf_in = "/fs1/results_dev/solid_hg38_masked/vcf/4-ffpe-blod1.agg.pon.vep.markgerm.vcf.gz"
#bcf_in = VariantFile("/fs1/results_dev/solid_hg38_masked/vcf/23PH12395-231010.agg.vcf.gz")
#bcf_in = "/fs1/results_dev/wgs/vcf/23MD09843.snv.rescored.sorted.vcf.gz"
#bcf_in = "test.vcf"
data = cmdvcf.open_vcf_stream(bcf_in)
count = 0
for var in data["vcf_object"].fetch():
    var_dict = cmdvcf.parse_variant(var,data["meta"],data["samples"])
    pprint(var_dict)
    if count == 5:
        exit()
    count+=1

## loads full vcf into memory
#meta,variants = cmdvcf.parse_vcf(bcf_in)

#pprint(vcf)
data_filtered = {}

# Fix/add various things to data structure
# for variant in vcf:

#     # Add sample ID to mongo-import-object
#     #$data->[$_]->{SAMPLE_ID} = $SAMPLE_ID; #bless(\$SAMPLE_ID, "MongoDB::BSON::String");

#     # Add type for pindel (special rule. should be done in pipeline for god sake!)
#     #$data->[$_]->{INFO}->{TYPE} = $data->[$_]->{INFO}->{SVTYPE} if $data->[$_]->{INFO}->{SVTYPE};

#     # fetch filters from filter-field add to import
#     #my @filters = split /;/, $data->[$_]->{FILTER};
#     #$data->[$_]->{FILTER} = \@filters;

#     # continue if any set criteria from config for assay
#     #next if grep /^(FAIL_NVAF|FAIL_LONGDEL)$/, @filters;
#     #next if grep /^FAIL_PON_/, @filters;
    
#     # add variant callers as a field into mongo-import-object
#     #my @found_in = split /\|/, $data->[$_]->{INFO}->{variant_callers};
#     #$data->[$_]->{INFO}->{variant_callers} = \@found_in;
    
    
#     # in pipe, make sure T-N samples always have Tumor as first sample. This WONT be adjusted for legacy shiet
#     first = 1
#     #for my $i ( 0.. scalar( @{ $data->[$_]->{GT} } )-1 ) {
	
# 	#my $gt = $data->[$_]->{GT}->[$i];
    
#     # check that data is compiant with coyote. VD DP AF AND GT required
#     #if( ( !defined($gt->{AF}) and !defined($gt->{VAF})) or !defined($gt->{DP}) or !defined($gt->{VD}) or !defined($gt->{GT}) ) {
    
#     # translate VAF to AF, don't know why this isnt part of pipeline?
# 	#if( defined $gt->{VAF} ) {
# 	#    $gt->{AF} = $gt->{VAF};
# 	#    delete $gt->{VAF};
# 	#}
#     # if first assign case, otherwise assign control, then reset first and continue with control
#     #$gt->{type} = ( $first ? "case" : "control" );

#     first = 0
#     # remove dataparts no longer needed
#     #delete $data->[$_]->{vcf_str};
#     #delete $data->[$_]->{INFO}->{'technology.illumina'};

#     # save variants surviving filters
#     #push @{$data_filtered}, $data->[$_];


    