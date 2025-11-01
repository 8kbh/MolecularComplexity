from flask import Flask, request, jsonify, Response, send_file
import json
import re
from flask_cors import CORS
import pandas as pd
from rdkit import Chem, RDLogger
from tqdm import tqdm
from utils import get_predictor
import logging
import random


RDLogger.DisableLog('rdApp.*')
logging.basicConfig(filename='add_mc.log', level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
predictor = get_predictor()

app = Flask(__name__)
CORS(app)


def validate_smiles(smiles: str) -> bool:
    if not isinstance(smiles, str):
        return False
    pattern = r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+$'
    if len(smiles) == 0 or not re.fullmatch(pattern, smiles):
        logging.warning(f"Invalid characters in SMILES string: '{smiles}'")
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.warning(f"RDKit could not parse the SMILES string: '{smiles}'")
        return False
    else:
        return True


@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    file_type = request.form.get('file_type')
    smiles_column = request.form.get('smiles_column')
    smiles_field = request.form.get('smiles_field')
    skip_header = int(request.form.get('skip_header', 0))
    processing_batch = int(request.form.get('processing_batch', 10))

    if file_type == 'csv' and not smiles_column:
        return jsonify({'error': 'SMILES column is required for CSV'}), 400
    if file_type == 'json' and not smiles_field:
        return jsonify({'error': 'SMILES field is required for JSON'}), 400

    temp_file_path = f'temp_{"".join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890") for _ in range(5))}.{file_type}'
    file.save(temp_file_path)

    output_file_path = temp_file_path.replace('.csv', '_mc.csv').replace('.json', '_mc.json')

    def generate():
        # Validation Progress
        if file_type == 'csv':
            df = pd.read_csv(temp_file_path, skiprows=skip_header, header=None)
            df = df.rename(columns={int(smiles_column) - 1: "smiles"})
            yield "data: Validating SMILES...\n\n"
            df["is_valid_smiles"] = df["smiles"].apply(validate_smiles)
            valid_smiles_df = df[df["is_valid_smiles"]].copy()
            valid_smiles_df["mc"] = 0.0
            yield "data: Validation complete\n\n"
        elif file_type == 'json':
            with open(temp_file_path, "r", encoding="utf-8") as f:
                data_json = json.load(f)
            yield "data: Validating SMILES...\n\n"
            valid_smiles = []
            for entry in tqdm(data_json):
                smiles = entry.get(smiles_field, "")
                valid_smiles.append((smiles, validate_smiles(smiles)))
            df = pd.DataFrame(valid_smiles, columns=["smiles", "is_valid_smiles"])
            valid_smiles_df = df[df["is_valid_smiles"]].copy()
            valid_smiles_df["mc"] = 0.0
            yield "data: Validation complete\n\n"

        # Processing Progress
        yield "data: Predicting molecular complexity...\n\n"
        for i in range(0, valid_smiles_df.shape[0], processing_batch):
            prediction = predictor.predict(valid_smiles_df["smiles"][i:i + processing_batch].to_list())
            valid_smiles_df.iloc[i:i + processing_batch, valid_smiles_df.columns.get_loc("mc")] = prediction
            progress = min(100, int((i + processing_batch) / valid_smiles_df.shape[0] * 100))
            yield f"data: {progress}\n\n"

        # Save the final result
        if file_type == 'csv':
            df.join(valid_smiles_df[["mc"]], how="left").to_csv(output_file_path, index=False, header=None)
        elif file_type == 'json':
            join_and_save_json(output_file_path, data_json, df.join(valid_smiles_df[["mc"]], how="left"))

        yield "data: Processing complete\n\n"
        yield f"data: {output_file_path}\n\n"

    return Response(generate(), mimetype='text/event-stream')


@app.route('/download/<filename>', methods=['GET'])
def download_file(filename):
    return send_file(filename, as_attachment=True)


def join_and_save_json(fp, dict_data, df):
    json_joined = dict_data.copy()
    for i, row in enumerate(df.itertuples(), start=0):
        if i < len(json_joined):
            json_joined[i]["_MolecularComplexity"] = row.mc
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(json_joined, f, ensure_ascii=False)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
