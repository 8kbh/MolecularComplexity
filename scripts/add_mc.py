import json
import argparse
import pandas as pd
from utils import get_predictor
from tqdm.autonotebook import tqdm
import logging
import csv
import os
import re
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog('rdApp.*')

logging.basicConfig(filename='add_mc.log', level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

predictor = get_predictor()

def process_single_smiles(smiles):
    is_valid_smiles = validate_smiles(str(smiles))
    if not is_valid_smiles:
        return None
    try:
        mc_prediction = predictor.predict([smiles,])
        return mc_prediction[0]
    except Exception as e:
        logging.info(smiles)
        logging.error(e)
        return None


def add_to_csv(filepath, rows):
    with open(filepath, "a", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def add_to_json(filepath, items):
    if os.path.exists(filepath):
        with open(filepath, "r", encoding="utf-8") as f:
            init_data = json.load(f)
        init_data += items
    else:
        init_data = items
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(init_data, f)


def ask_append_replace(filename):
    while True:
        user_input = input(f"File {filename} already exists. Replace file (R) or append to file (A)? ").strip().lower()
        if user_input == 'r':
            return True
        elif user_input == 'a':
            return False
        else:
            print("Invalid input. Please enter 'R' or 'A'. Or interrupt to exit")


def validate_smiles(smiles: str) -> bool:
    # Step 1: Regex for allowed SMILES symbols (simplified version)
    pattern = r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+$'
    if len(smiles) == 0 or not re.fullmatch(pattern, smiles):
        logging.warning(f"Invalid characters in SMILES string: '{smiles}'")
        return False

    # Step 2: Try parsing with RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.warning(f"RDKit could not parse the SMILES string: '{smiles}'")
        return False
    else:
        return True


def main():

    parser = argparse.ArgumentParser(description="Add Molecular Complexity value to existing database")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--db_csv", "-t", type=str, help="Path to the CSV database file")
    group.add_argument("--db_json", "-d", type=str, help="Path to the JSON database file")

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument("--smiles_column", "-c", type=int, help="Number of the column in which SMILES are stored (for CSV)")
    group2.add_argument("--smiles_field", "-f", type=str, help="Key of the dictionary where SMILES are stored (for JSON)")

    parser.add_argument("--skip_header", "-s", type=int, help="Number of header rows to be skipped (for CSV)")
    parser.add_argument("--writing_batch", "-w", type=int, default=10, help="Batch size for writing results in file. Set 0 to write at the end")
    # Define output argument
    parser.add_argument("--output", "-o", type=str, help="Name of the output file (optional)")

    args = parser.parse_args()

    # Validate SMILES column/field based on input type
    if args.db_csv and not args.smiles_column:
        parser.error("--smiles_column (-c) is required when using --db_csv (-t)")
    if args.db_json and not args.smiles_field:
        parser.error("--smiles_field (-f) is required when using --db_json (-d)")

    # Generate output filename if not provided
    if not args.output:
        if args.db_csv:
            args.output = args.db_csv.replace('.csv', '_mc.csv')
        elif args.db_json:
            args.output = args.db_json.replace('.json', '_mc.json')

    if not args.skip_header:
        args.skip_header = 0

    # Your processing logic here
    print(f"Input file: {args.db_csv or args.db_json}")
    print(f"SMILES column/field: {args.smiles_column or args.smiles_field}")
    print(f"Output file: {args.output}")
    if os.path.exists(args.output):
        to_replace = ask_append_replace(args.output)
        if to_replace:
            os.remove(args.output)
    print("Start processing")

    # process csv
    if args.db_csv:
        df = pd.read_csv(args.db_csv, skiprows=args.skip_header, header=None)
        if args.writing_batch == 0:
            tqdm.pandas()
            df[df.shape[1]] = df[args.smiles_column - 1].progress_apply(process_single_smiles)
            df.to_csv(args.output, index=False, na_rep='N/A', encoding='utf-8', header=False)
        else:
            temp_data = []
            for index, row in tqdm(df.iterrows(), total=df.shape[0]):
                row_list = row.to_list()
                if index + 1 >= args.writing_batch and index % args.writing_batch == 0:
                    add_to_csv(args.output, temp_data)
                    temp_data = []
                prediction = process_single_smiles(row_list[args.smiles_column - 1])
                temp_data.append(row_list + [prediction, ])
            if len(temp_data) > 0:
                add_to_csv(args.output, temp_data)


    # process json
    if args.db_json:
        with open(args.db_json, "r", encoding="utf-8") as f:
            data = json.load(f)
        temp_data = []
        
        if args.writing_batch == 0:
            args.writing_batch = len(data)

        for index, sample in enumerate(tqdm(data)):
            if index + 1 >= args.writing_batch and index % args.writing_batch == 0:
                add_to_json(args.output, temp_data)
                temp_data = []

            smiles = sample.get(args.smiles_field, None)
            temp_data.append(data[index])

            if not smiles:
                temp_data[-1]["_MolecularComplexity"] = None
            else:
                temp_data[-1]["_MolecularComplexity"] = process_single_smiles(smiles)

        if len(temp_data) > 0:
            add_to_json(args.output, temp_data)

    print("Done")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()
