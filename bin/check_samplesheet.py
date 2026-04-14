#!/usr/bin/env python3

# if no header -> stop the pipeline
# if no assay -> stop the pipeline
# 
# should  have valid sampleid
# should have valid fastq files and link 
# if paired check for tumor and normal 
'''
./check_samplesheet.py -c headerless.csv -o samplecheck.txt
'''

import csv
import argparse
import json
from pathlib import Path


def load_allowed_assays(json_path):
    with Path(json_path).open() as handle:
        dataset = json.load(handle)

    if isinstance(dataset, dict):
        assays = dataset.get("cmd_assays", [])
    elif isinstance(dataset, list):
        assays = dataset
    else:
        raise ValueError("CMD assay JSON must be a list or contain a 'cmd_assays' list")

    return {assay.strip() for assay in assays if assay and assay.strip()}

def process_linescsv_file(test):
    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)
        nrows = len(list(csvFile))
        return (nrows)
    
def process_csv_file(test, nrows, allowed_assays):
    testAssay = []
    testId = []
    testType = []

    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)

        if nrows == 0:
            return None
        else:
            for row in csvFile:
                assay = row["assay"].strip()
                if len(assay) == 0 or assay not in allowed_assays:
                    return None
                else:
                    testAssay.append(assay)

                if len(row["id"]) == 0:
                    return None
                else:
                    testId.append(row["id"])

                if len(row["type"]) == 0:
                    return None
                else:
                    testType.append(row["type"])
        return nrows, testAssay, testId, testType

def writeFile (result,output):
    if result is not None:
        with open(output, mode= 'w') as outFile:
            lineCount, testAssay, testId, testType = result
            outFile.write(str(lineCount))
            outFile.write(str(testAssay))
            outFile.write(str(testId))
            outFile.write(str(testType))
            return outFile

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--csv', dest = 'csv',  default = "test.csv", help = "Sample Sheet csv for the nextflow")
    parser.add_argument('-o', '--out', dest = 'output', default = "result", help = "Input csv structure and content signal")
    parser.add_argument(
        '--cmd-assays-json',
        dest='cmd_assays_json',
        help='Path to a JSON file containing allowed CMD assay names.',
    )

    args = parser.parse_args()
    inputCsv = args.csv
    outputCsv = args.output
    allowed_assays = load_allowed_assays(args.cmd_assays_json)
    count = process_linescsv_file(inputCsv)
    
    result = process_csv_file(inputCsv, count, allowed_assays)

    writeFile(result, outputCsv)

if __name__ == "__main__":
    Main()


