from features import smiles_to_morganfingerprints
from bentoml.io import JSON, NumpyNdarray
from pydantic import BaseModel
from rdkit.Chem import Mol
import numpy as np
import bentoml

class Molecule(BaseModel):
    smiles:str = ''

bbb_clf_runner = bentoml.sklearn.get("bbb-model:latest").to_runner()
svc = bentoml.Service("bbb-model", runners=[bbb_clf_runner])

input_spec = JSON(pydantic_model=Molecule)

@svc.api(input=input_spec, output=NumpyNdarray())
def classify(molecule: Molecule) -> np.ndarray:
    features = smiles_to_morganfingerprints(molecule.smiles)
    result = bbb_clf_runner.predict.run(np.array([features]))
    return result