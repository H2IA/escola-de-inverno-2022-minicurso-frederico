from rdkit import Chem
from rdkit.Chem import AllChem
from argparse import ArgumentParser
from tqdm import tqdm
import numpy as np
import pandas as pd

def main():
    argument_parser = ArgumentParser()
    argument_parser.add_argument('--input', default='data/raw/b3db.tsv')
    argument_parser.add_argument('--output', default='data/features/b3db.csv')
    arguments = argument_parser.parse_args()

    print('Computing morgan fingerprints')

    df_raw = pd.read_csv(arguments.input, sep='\t')
    df_features = compute_features(df_raw)
    df_features.to_csv(arguments.output, index=False)

def compute_features(df_raw:pd.DataFrame) -> pd.DataFrame:

    X = []
    y = []

    for r, row in tqdm(df_raw.iterrows(), total=df_raw.shape[0]):
        
        fingerprints = smiles_to_morganfingerprints(row.SMILES)
        X.append(fingerprints)
        y.append(1 if row['BBB+/BBB-'] == 'BBB+' else 0)

    X = np.array(X)
    y = np.array(y)

    df_features = pd.DataFrame(X)
    df_features['label'] = y
    return df_features

def smiles_to_morganfingerprints(smiles:str) -> np.array:
    mol = Chem.MolFromSmiles(smiles)
    fingerprints = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits=512, bitInfo={})
    return np.array(fingerprints)

if __name__ == '__main__':
    main()