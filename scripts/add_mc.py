import argparse
import csv
import json
import logging
import os
import re

import pandas as pd
from rdkit import Chem, RDLogger
from tqdm.autonotebook import tqdm
from utils import get_predictor

tqdm.pandas()

RDLogger.DisableLog('rdApp.*')

logging.basicConfig(filename='add_mc.log', level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

predictor = get_predictor()


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


def join_and_save_json(fp, dict, df):
    json_joined = dict.copy()

    for i, row in enumerate(df.itertuples(), start=0):
        json_joined[i]["_MolecularComplexity"] = row.mc

    with open(fp, "w", encoding="utf-8") as f:
        json.dump(json_joined, f, ensure_ascii=False)


def ask_replace(filename):
    user_input = input(f"File {filename} already exists. Replace file? (N/y) ").strip().lower()
    if user_input == 'y':
        return True
    else:
        print("Exiting")
        print("NOTE: you can specify output file name with option --output/-o")
        exit()


def validate_smiles(smiles: str) -> bool:
    if not isinstance(smiles, str):
        return False

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

    parser.add_argument("--skip_header", "-s", type=int, default=0, help="Number of header rows to be skipped (for CSV)")
    parser.add_argument("--processing_batch", "-p", type=int, default=10, help="Batch size for processing. Set 0 to process all simultaneously")
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

    # Your processing logic here
    print(f"Input file: {args.db_csv or args.db_json}")
    print(f"SMILES column/field: {args.smiles_column or args.smiles_field}")
    print(f"Output file: {args.output}")

    if os.path.exists(args.output):
        to_replace = ask_replace(args.output)
        if to_replace:
            os.remove(args.output)
    print("Start processing")

    # process csv
    if args.db_csv:
        df = pd.read_csv(args.db_csv, skiprows=args.skip_header, header=None)
        print("Validating SMILESes...")

        df = df.rename(columns={args.smiles_column - 1: "smiles"})
        df["is_valid_smiles"] = df["smiles"].progress_apply(validate_smiles)

    if args.db_json:
        with open(args.db_json, "r", encoding="utf-8") as f:
            data_json = json.load(f)

        print("Validating SMILESes...")
        valid_smiles = []
        for entry in tqdm(data_json):
            smiles = entry.get(args.smiles_field, "")
            valid_smiles.append((smiles, validate_smiles(smiles)))

        df = pd.DataFrame(valid_smiles, columns=["smiles", "is_valid_smiles"])

    valid_smiles_df = df[df["is_valid_smiles"]].copy()
    valid_smiles_df["mc"] = 0.0

    print()
    print(f"Predicting molecular complexity (batched: {args.processing_batch}its)")
    for i in tqdm(range(0, valid_smiles_df.shape[0], args.processing_batch)):
        prediction = predictor.predict(valid_smiles_df["smiles"][i:i + args.processing_batch].to_list())
        valid_smiles_df.iloc[i:i + args.processing_batch, valid_smiles_df.columns.get_loc("mc")] = prediction
        last_calculated_row = valid_smiles_df.iloc[i:i + args.processing_batch].index[-1]

        df_joined = df.join(valid_smiles_df[["mc"]], how="left").loc[:last_calculated_row]
        # print(df_joined)
        if args.db_csv:
            df_joined.to_csv(args.output, index=False, header=None)

        if args.db_json:
            join_and_save_json(args.output, data_json, df_joined)

    if (last_calculated_row < df.shape[0] - 1):
        if args.db_csv:
            df.join(valid_smiles_df[["mc"]], how="left").to_csv(args.output, index=False, header=None)

        if args.db_json:
            join_and_save_json(args.output, data_json, df.join(valid_smiles_df[["mc"]], how="left"))


if __name__ == "__main__":
    main()
