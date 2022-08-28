# H2IA: Botando pra rodar

Neste projeto utilizaremos dados do B3DB, um banco de dados de moléculas com diferentes níveis de permeabilidade para a barreira cérebro-sangue (*blood-brain barrier*), para produzir um modelo preditivo capaz de estimar esta propriedade. Além disso, disponibilizaremos o modelo na forma de um *API REST* no *Google Cloud Platform* usando o serviço *App Run*. 

## Setup

```
$ conda env create --file environment.yml
$ conda activate escola-de-inverno-2022-minicurso-frederico
```

## Download data

```
$ make download
```



## Produzindo as features

```
$ make features
```

Esta regra executará o *script* `features.py`, que processará o arquivo `data/raw/b3db.tsv` e produzirá um arquivo novo contendo as *Morgan Fingerprints* das moléculas. Os resultados serão salvos no arquivo `data/features/b3db.csv`. As *Morgan Fingerprints* são *features* binárias que presentam a presença ou ausência de alguma subestrutura dentro das moléculas. Este processo é realizado pela função `smiles_to_morganfingerprints(...)`, que recebe uma representação da molécula em SMILES e retorna um `numpy.array` com as *features* computadas.

```python
def smiles_to_morganfingerprints(smiles:str) -> np.array:
    mol = Chem.MolFromSmiles(smiles)
    fingerprints = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits=512, bitInfo={})
    return np.array(fingerprints)
```

## Treinando o modelo

```
$ make train
```

Esta regra executará o *script* `train.py`, que treinará um modelo do tipo *ExtraTreeClassifier* do *SciKit-Learn* a partir do conjunto de dados do
B3DB. Após o treinamento e uma validação, o modelo será serializado em um
arquivo `.pickle` localizando na pasta `data/models` e também salvo
através do *BentoML* com o nome `bbb-model`.

```python
with open(arguments.output, 'wb') as writer:
    writer.write(pickle.dumps(model))
print(f"Model file: {arguments.output}")

saved_model = bentoml.sklearn.save_model('bbb-model', model)
print(f"Model tag: {saved_model}")
```

## Empacotando o modelo com BentoML

```
$ make build
```

```
$ make container
```

## GCP Deployment

### Authentication

```
$ gcloud auth login
```

### Create project

```
$ gcloud projects create escola-de-inverno-2022-minicurso-frederico-test
```

### Select active project

```
$ gcloud config set project escola-de-inverno-2022-minicurso-frederico
```

### Set
