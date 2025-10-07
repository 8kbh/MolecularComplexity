import json
import argparse
import pandas as pd
from utils import get_predictor
from tqdm.autonotebook import tqdm


predictor = get_predictor()
def process_single_smiles(smiles):
    try:
        mc_prediction = predictor.predict([smiles,])
        return mc_prediction[0]
    except Exception as e:
        print(e)
        return None


def main():
    parser = argparse.ArgumentParser(description="Add Molecular Complexity value to existing database")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--db_csv", "-t", type=str, help="Path to the CSV database file")
    group.add_argument("--db_json", "-d", type=str, help="Path to the JSON database file")

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument("--smiles_column", "-c", type=int, help="Number of the column in which SMILES are stored (for CSV)")
    group2.add_argument("--smiles_field", "-f", type=str, help="Key of the dictionary where SMILES are stored (for JSON)")

    parser.add_argument("--skip_header", "-s", type=int, help="Number of header rows to be skipped (for CSV)")
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
    print("Start processing")

    # process csv
    if args.db_csv:
        tqdm.pandas()
        df = pd.read_csv(args.db_csv, skiprows=args.skip_header, header=None)
        # print(df)
        df[df.shape[1]] = df[args.smiles_column - 1].progress_apply(process_single_smiles)
        df.to_csv(args.output, index=False, na_rep='N/A', encoding='utf-8', header=False)

    # process json
    if args.db_json:
        with open(args.db_json, "r", encoding="utf-8") as f:
            data = json.load(f)
        for i, sample in enumerate(tqdm(data)):
            smiles = sample.get(args.smiles_field, None)
            if not smiles:
                data[i]["_MolecularComplexity"] = None
            else:
                data[i]["_MolecularComplexity"] = process_single_smiles(smiles)
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    print("Done")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()
