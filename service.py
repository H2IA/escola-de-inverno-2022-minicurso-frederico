from features import smiles_to_morganfingerprints
from bentoml.io import Text, NumpyNdarray
from rdkit.Chem import Mol
import numpy as np
import bentoml

bbb_clf_runner = bentoml.sklearn.get("bbb-model:latest").to_runner()

svc = bentoml.Service("bbb-model", runners=[bbb_clf_runner])

@svc.api(input=Text(), output=NumpyNdarray())
def classify(smiles: str) -> np.ndarray:
    features = smiles_to_morganfingerprints(smiles)
    result = bbb_clf_runner.predict.run(np.array([features]))
    return result