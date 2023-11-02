from pysam import VariantFile

def header_dict(vcf_object):
    """
    interpret meta information from header. Especially for CSQ-fields
    """
    meta = {}
    header = vcf_object.header
    
    csq = header.info['CSQ']
    csq = str(csq.description).split(" ")
    csq = csq.pop().split("|")
    meta['CSQ'] = csq

    header_string = ""
    for x in vcf_object.header.records:
        header_string = header_string+str(x)
    meta["header_str"] = header_string
    meta["pysam_obj"] = vcf_object.header
    return meta

def parse_vcf(infile):
    """
    Store VCF records as list of dicts. This function might cause memory issues
    Use next_var() to go through one by one

    This function needs to return meta-dict,var-list,samples(?)
    """
    # Initialize an empty list to store dictionaries
    bcf_in = VariantFile(infile)
    meta = header_dict(bcf_in)
    samples = list((bcf_in.header.samples))
    vcf_records = []
    # Iterate through the VCF records and convert them to dictionaries
    for record in bcf_in.fetch():
        vcf_dict = parse_variant(record,meta,samples)
        vcf_records.append(vcf_dict)
    return meta,vcf_records

def fix_gt(gt_dict,samples):
    """
    combine format with sample info of format. Return array of individuals with GT dicts
    """
    gt_list = []
    for sample in samples:
        format_list = list(dict(gt_dict[sample]).keys())
        sample_dict = dict(gt_dict[sample])
        sample_dict["_sample_id"] = sample
        sample_dict['GT'] = "/".join(str(x) for x in list(sample_dict['GT']))
        gt_list.append(sample_dict)
    return gt_list,format_list

def parse_variant(record,meta,samples):
    """
    Changes from Pysam-keys to CMD-keys per variant and returns dict of variant
    """
    varid = record.id
    if varid == None:
        varid = "."
    gt_list,format_list = fix_gt(dict(record.samples),samples)
    vcf_dict = {
        "CHROM": record.chrom,
        "POS": record.pos,
        "ID": varid,
        "REF": record.ref,
        "ALT": ','.join(list(record.alts)),
        "QUAL": record.qual,
        "FILTER": ';'.join(list(record.filter)),
        "INFO": fix_info(dict(record.info),meta),
        "FORMAT": format_list,
        "GT": gt_list
    }
    return vcf_dict

def fix_info(infodict,meta):
    """
    make CSQ string a list of key:value pairs per transcript.
    """
    new_info_dict = {}
    for key in infodict:
        if key == 'CSQ':
            transcipts = list(infodict[key])
            #new_info_dict['CSQ'] = csq(transcipts,meta)
        else:
            # just pass all other key:value to new dict
            # if pysam made a real list of comma separated values store as string not array
            if isinstance(infodict[key],bool):
                value = 1
            elif len(list(infodict[key])) > 1:
                value = ",".join(list(infodict[key]))
            else:
                value = "".join(list(infodict[key]))
            new_info_dict[key] = value
    return new_info_dict

def csq(transcipts,meta):
    """
    Store CSQ fields as list of dicts per transcript
    """
    # split csq on ; save as append as list to CSQ key  
    csq_list = []
    csq_dict = {}
    for trans in transcipts:
        trans_anno = trans.split('|')
        for anno in range(0,len(trans_anno)-1):
            if meta['CSQ'][anno] == "Consequence":
                conq_list = trans_anno[anno].split('&')
                csq_dict[meta['CSQ'][anno]] = conq_list
            else:
                csq_dict[meta['CSQ'][anno]] = trans_anno[anno]
        csq_list.append(csq_dict)
    return csq_list

def open_vcf_stream(infile):
    bcf_in = VariantFile(infile)
    data = {}
    meta = header_dict(bcf_in)
    data["meta"] = meta
    data["vcf_object"] = bcf_in
    data["samples"] = list((bcf_in.header.samples))
    return data
