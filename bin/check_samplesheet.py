#!/usr/bin/env python3

# Released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.


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

CMDAssay = [    
                "gmslymphomav3-0", 
                "GMSMyeloidv1-0", 
                "gmssolidtumorv3-0", 
                "PARPinhibv1-0" 
            ]
def process_linescsv_file(test):
    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)
        nrows = len(list(csvFile))
        return (nrows)
    
def process_csv_file(test,nrows):
    testAssay = []
    testId = []
    testType = []

    with open(test, mode='r') as file:
        csvFile = csv.DictReader(file)

        if nrows == 0:
            return None
        else:
            for row in csvFile:
                if len(row["assay"]) == 0 and row["assay"][0] not in CMDAssay:
                    return None
                else:
                    testAssay.append(row["assay"])

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

    args = parser.parse_args()
    inputCsv = args.csv
    outputCsv = args.output
    count = process_linescsv_file(inputCsv)
    
    result = process_csv_file(inputCsv,count)

    writeFile(result, outputCsv)

if __name__ == "__main__":
    Main()





