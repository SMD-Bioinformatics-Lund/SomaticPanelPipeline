#!/usr/bin/env python3

import argparse
import json

# Function to load JSON data from a file
def load_json_file(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Function to save JSON data to a file
def save_json_file(data, file_path):
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

# Main function to combine the Info JSON and the genotype JSON
def combine_json_files(info_json_file, genotype_json_file, partner_run_json_file, output_file):
    # Load the Info and genotype JSON files
    info_data = load_json_file(info_json_file)
    genotype_data = load_json_file(genotype_json_file)
    
    # Combine the data by appending the genotype JSON
    combined_data = {
        **info_data,
        "id_snp_genotypes": genotype_data
    }

    # If partner_run_json_file is provided, include it in the combined data
    if partner_run_json_file:
        partner_run = load_json_file(partner_run_json_file)
        combined_data["partner_info"] = partner_run

    # Save the combined data to the output file
    save_json_file(combined_data, output_file)

# Argument parser setup
def parse_args():
    parser = argparse.ArgumentParser(description='Append genotype JSON into Info JSON and save the result.')

    # Positional arguments are kept for backward compatibility with existing calls.
    parser.add_argument('info_json_file_pos', nargs='?', help='Path to the Info JSON file')
    parser.add_argument('genotype_json_file_pos', nargs='?', help='Path to the genotype JSON file')
    parser.add_argument(
        '--partner-run-json-file',
        '--partner_run_json_file',
        dest='partner_run_json_file',
        help='Path to the partner sample run JSON file (optional)',
        default=None,
    )
    parser.add_argument('output_file_pos', nargs='?', help='Path to the output combined JSON file')
    parser.add_argument('--info-json', dest='info_json_file', help='Path to the Info JSON file')
    parser.add_argument('--genotype-json', dest='genotype_json_file', help='Path to the genotype JSON file')
    parser.add_argument('--out', dest='output_file', help='Path to the output combined JSON file')

    args = parser.parse_args()
    args.info_json_file = args.info_json_file or args.info_json_file_pos
    args.genotype_json_file = args.genotype_json_file or args.genotype_json_file_pos
    args.output_file = args.output_file or args.output_file_pos
    missing = [
        name for name in ("info_json_file", "genotype_json_file", "output_file")
        if getattr(args, name) is None
    ]
    if missing:
        parser.error(f"missing required arguments: {', '.join(missing)}")
    return args

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # Call the function to combine the JSON files
    combine_json_files(args.info_json_file, args.genotype_json_file, args.partner_run_json_file, args.output_file)
