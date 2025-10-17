# Code for the paper "Digitization of molecular complexity with machine learning"

This repository contains the code for the paper ["Digitization of molecular complexity with machine learning"](https://doi.org/10.1039/D4SC07320G).

## Installation

### Conda

After cloning the repository, install conda environment with all the dependencies.
```bash
conda env create -f environment.yml
conda activate mc
pip install -e .
```

### [uv](https://docs.astral.sh/uv/)
Alternatively you can use uv project manager
```bash
uv sync
uv pip install -e .
```
To launch scripts use `uv run`

Ex.:
```bash
uv run scripts/add_mc.py --db_csv data/example_database.csv --skip_header 1 --smiles_column 5
```

### Necessary data

To use the code you will also need to download the zip file containing all the relevant data from [here](https://drive.google.com/file/d/1JatMaqwZvjAmyE-rznIihQa6WF8U4poy/view?usp=sharing)

Or using
```bash
curl -O https://pubfs.8kbh.ru/files/data.zip
```

And and unzip it in the project's folder
```bash
unzip data.zip
rm data.zip
```

## Usage

### Simple calculate
To calculate the molecular complexity of a list of smiles, run the following command:
```bash
python scripts/calculate.py --txt_with_smiles data/example_smiles.txt
```
Where `data/example_smiles.txt` is a text file with one SMILES string per line. The output will be saved in `example_mc.json` file.

To plot the results from the paper, run the scripts from `scripts` folder that begin with `plot_`:
```bash
python scripts/plot_fda.py
python scripts/plot_synthesis.py
```

### Add Molecular Complexity value to existing CSV/JSON file

```bash
uv run scripts/add_mc.py [OPTIONS]
```

| Option                     | Description                                      | Default Value |
|----------------------------|--------------------------------------------------|---------------|
| `--db_csv -t`              | Path to the CSV database file                    |               |
| `--db_json -d`             | Path to the JSON database file                   |               |
| `--smiles_column -c`       | Column number where SMILES are stored (for CSV)  |               |
| `--smiles_field -f`        | Dictionary key where SMILES are stored (for JSON)|               |
| `--skip_header -s`         | Number of header rows to skip (for CSV)          | 0             |
| `--processing_batch -p`    | Batch size for processing                        | 10            |
| `--output -o`              | Name of the output file (optional)               |               |

