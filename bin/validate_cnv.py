import sys
import logging
import json
from pprint import pprint
"""
read json containing called variants. The interesting genes should be annotated
with class and cnv_type

coordinates("chrom:start-end"): {
callers: string
ratio: float(log2-ratio, optional)
size: int,
PR: string(optional),
SR: string(optional),
genes: [
    { 
    gene: string(name of gene),
    class, string(somatic/research),
    cnv_type: string(unspec/amp/del)
    },
    {
    gene: string(name of gene)
    }
    ]
},
...
 read json containing variants to be recalled
      gene-name: {
        type : amp/del,
        estimated_copies: int,
        size: int(optional)
      },
      ...
"""
def read_json(json_file:str):
    """
    read json return dict
    """
    with open(json_file, 'r') as json_object:
        json_dict = json.load(json_object)
    return json_dict


def recall(calls:dict,definitions:dict):
    """
    iterate over call if existing in definitions add evidence to
    definintions dict. Present how well recall performed
    """
    for variant in calls:
        for gene in calls[variant]["genes"]:
            if "cnv_type" in gene:
                if gene["gene"] in definitions:
                    recall_dict = analyze_recall(calls[variant],definitions[gene["gene"]])
                    if recall_dict:
                        definitions[gene["gene"]]["recalled"] = recall_dict


    return definitions

def analyze_recall(variant,defined_vals):
    """
    summerize recall data. DUP=DUP, CN~=CN
    what callers: add as list
    optional size, is the size within expected sizes
    """
    callers = variant['callers'].split('-')
    
    if "recalled" in defined_vals:
        summary = defined_vals["recalled"]
    else:
        summary = {}
    for caller in callers:
        tmp = {}
        if "ratio" in variant:
            cn = 2*(2**variant["ratio"])
            cn = round(cn)
            if cn >= 2:
                if defined_vals["type"] == "amp":
                    if cn == defined_vals["estimated_copies"][0]:
                        tmp["cn_score"] = 3
                    elif cn > defined_vals["estimated_copies"][1]:
                        tmp["cn_score"] = 2
                    else:
                        tmp["cn_score"] = 1
                else:
                    return False
            else:
                if defined_vals["type"] == "del":
                    tmp["cn_score"] = 3
                else:
                    return False
        summary[caller] = tmp        
    return summary
def summerize_results(results:dict):
    """
    results summerized. % recalled, score per variant + average score
    """
    
    #count = sum(specific_key in inner_dict for inner_dict in results.values())
    to_be_recalled = len(results.keys())
    recalled = 0
    tot_var_score = 0
    for var in results:
        var_score = 0
        varcallers = list(results[var]["recalled"].keys())
        if "recalled" in results[var]:
            recalled +=1
            for caller in results[var]["recalled"]:
                var_score +=1
                var_score = var_score + results[var]["recalled"][caller]["cn_score"]
                tot_var_score += var_score
        print(var,"&".join(varcallers),var_score)
    print(str(100*recalled/to_be_recalled),"%",str(tot_var_score),str(tot_var_score/recalled))




definitions = read_json(sys.argv[1])
calls = read_json(sys.argv[2])
results = recall(calls,definitions)
summerize_results(results)