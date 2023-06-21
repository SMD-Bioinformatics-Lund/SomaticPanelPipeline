    #!/usr/bin/python3

import argparse
import pprint
import json

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--msi_s')
parser.add_argument('-p', '--msi_p')
parser.add_argument('-d', '--hrd')
parser.add_argument('-t', '--tmb')
parser.add_argument('-o', '--out')
parser.add_argument('-i', '--id')
args = parser.parse_args()
json_out = open(str(args.out), "w")

def read_msi(file):
    next(file)
    for line in file:
        info = line.split("\t")
        tot = info[0]
        som = info[1]
        perc = info[2]
        msi_dict = {
            'tot': int(tot),
            'som': int(som),
            'perc': float(perc.strip())
        }
    return msi_dict

def read_hrd(file):
    next(file)
    for line in file:
        info = line.split("\t")
        tai = info[2]
        lst = info[3]
        hrd = info[1]
        sum = info[4]
        hrd_dict = {
                'tai': int(tai),
                'hrd': int(hrd),
                'lst': int(lst),
                'sum': int(sum.strip())
        }
        return hrd_dict
        

biomarkers = {}
sample_id = args.id
# MSI-S single sample predicted MSI-score
if (args.msi_s):
    in_file = open(args.msi_s)
    msi_dict = read_msi(in_file)
    biomarkers['MSIS'] = msi_dict 

# MSI-S paired sample predicted MSI-score
if (args.msi_p):
    in_file = open(args.msi_p)
    msi_dict = read_msi(in_file)
    biomarkers['MSIP'] = msi_dict 

# HRD
if (args.hrd):
    in_file = open(args.hrd)
    hrd_dict = read_hrd(in_file)
    biomarkers['HRD'] = hrd_dict 

# TMB 
# To be done?
biomarkers['name'] = sample_id
json_object = json.dumps(biomarkers, indent = 4) 
json_out.write(json_object)