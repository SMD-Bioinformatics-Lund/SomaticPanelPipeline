from pysam import VariantFile
from typing import Any


def parse_vcf(infile: object):
    """
    Store VCF records as list of dicts. This function might cause memory issues
    Use vcf_object.fetch() to go through one by one

    This function returns header,var-list
    """
    # Initialize an empty list to store dictionaries
    bcf_in = VariantFile(infile)
    header = bcf_in.header
    vcf_records = []
    # Iterate through the VCF records and convert them to dictionaries
    for record in bcf_in.fetch():
        vcf_dict = parse_variant(record, header)
        vcf_records.append(vcf_dict)
    return header, vcf_records


def fix_gt(gt_dict: dict, samples: list):
    """
    Stores FORMAT-field information and adds _samlple_id to the GT-field dict. Return array of individuals
    with GT dicts. It is important that the sample-order maintains intact, tumor-sample usually is index 0
    for coyote imports.
    Stringifies values of pysam dicts. Special rule for GT, 0/1 is format we use to decide
    genotype. pysam uses (0,1).
    """
    gt_list = []
    for sample in samples:
        format_list = list(dict(gt_dict[sample]).keys())
        sample_dict = dict(gt_dict[sample])
        for key in sample_dict:
            if key == "GT":
                sample_dict["GT"] = "/".join(str(x) for x in list(sample_dict["GT"]))
            else:
                sample_dict[key] = unravel_tuples(sample_dict[key])
        sample_dict["_sample_id"] = sample
        gt_list.append(sample_dict)
    return gt_list, format_list


def parse_variant(record, header):
    """
    Changes from Pysam-keys to CMD-keys per variant and returns dict of variant
    """
    varid = record.id
    if varid == None:
        varid = "."
    gt_list, format_list = fix_gt(dict(record.samples), list(header.samples))
    vcf_dict = {
        "CHROM": record.chrom,
        "POS": record.pos,
        "ID": varid,
        "REF": record.ref,
        "ALT": ",".join(list(record.alts)),
        "QUAL": record.qual,
        "FILTER": ";".join(list(record.filter)),
        "INFO": fix_info(dict(record.info), header),
        "FORMAT": format_list,
        "GT": gt_list,
    }
    return vcf_dict


def fix_info(infodict, header):
    """
    make CSQ string a list of key:value pairs per transcript.
    """
    new_info_dict = {}
    for key in infodict:
        if key == "CSQ":
            transcripts = list(infodict[key])
            new_info_dict["CSQ"] = csq(transcripts, header)
        elif key == "ANN":
            transcripts = list(infodict[key])
            new_info_dict["ANN"] = snpeff(transcripts, header)
        else:
            new_info_dict[key] = unravel_tuples(infodict[key])
    return new_info_dict


def csq(transcripts: dict, header: object):
    """
    Store CSQ fields as list of dicts per transcript. Get Key-names from
    vcf_object.header (meta) info['CSQ'].description
    """
    csq = header.info["CSQ"]
    csq = str(csq.description).split(" ")
    csq_keys = csq.pop().split("|")
    csq_list = []
    for trans in transcripts:
        csq_dict = {}
        trans_anno = trans.split("|")
        for key, anno in zip(csq_keys, trans_anno):
            if key == "Consequence":
                conq_list = anno.split("&")
                csq_dict[key] = conq_list
            else:
                csq_dict[key] = anno
        csq_list.append(csq_dict)
    return csq_list


def snpeff(transcripts: dict, header: object):
    """
    Store ANN fields as list of dicts per transcript. Get Key-names from
    vcf_object.header (meta) info['ANN'].description
    """
    snpeff = header.info["ANN"]
    snpeff_keys = snpeff.description.split(" | ")
    snpeff_keys[0] = "Allele"
    snpeff_list = []
    for trans in transcripts:
        snpeff_dict = {}
        trans_anno = trans.split("|")
        for key, anno in zip(snpeff_keys, trans_anno):
            if key == "Annotation":
                conq_list = anno.split("&")
                snpeff_dict[key] = conq_list
            else:
                snpeff_dict[key] = anno
        snpeff_list.append(snpeff_dict)
    return snpeff_list


def unravel_tuples(value: Any):
    """
    Objects in pysam can have different classes, CMD-style VCF dicts should mostly be strings
    If this is correct is not relevant, this module is trying to emulate what vcf2.pm does.
    This function will make values appear as strings.

    Tuples -> join with "," except if not str and make it such
    """
    if isinstance(value, tuple):
        try:
            value = ",".join(list(value))
        except:
            value = ",".join(str(x) for x in list(value))
    return value
