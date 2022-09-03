# H2IA: Botando pra rodar

Neste projeto utilizaremos dados do [B3DB](https://github.com/theochem/B3DB), um banco de dados de moléculas com diferentes níveis de permeabilidade para a barreira cérebro-sangue (*blood-brain barrier*), para produzir um modelo preditivo capaz de estimar esta propriedade. Além disso, disponibilizaremos o modelo na forma de um *API REST* no *Google Cloud Platform* usando o serviço *Cloud Run* e as ferramenta *BentoML* e *Terraform*. 

E-mail: [fred.s.kremer@gmail.com](mailto:fred.s.kremer@gmail.com)

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

## Rodando o modelo localmente

## Empacotando o modelo com BentoML e rodando localmente

```
$ make build serve
```

A regra `build` executará o comando `bentoml build`, salvando a versão
atual do modelo treinado. Já a `serve` irá iniciar um servidor local
para disponibilizar o modelo em um API REST.

## Deployment

**Nota:** Para enviar o modelo para o GCP é necessário ter o Terraform instalado. 

```
$ bash scripts/deploy_to_gcp.sh
```

Este *script* executará os seguintes comandos:

```bash
#!/usr/bin/env bash

set -e

# GCP configuration

gcloud auth login
gcloud projects create h2ia-2022-frederico || echo "Project already exists"
gcloud config set project h2ia-2022-frederico

# Terraform config generation

bentoctl operator install google-cloud-run
bentoctl init
bentoctl build -b bbb-model:latest -f deployment_config.yaml

# GCP deployment with terraform

terraform init
terraform apply -var-file=bentoctl.tfvars -auto-approve
```

## Testando a API

Agora vamos testar a API usando a estrutura do fármaco Levodopa, usando no tratamento 
da doença de Parkinson. Este fármaco precisa agir no sistema nervoso central, sendo
assim necessária a sua passagem pela BBB.

```
$ molecule='N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O' # SMILES do levodopa
$ curl \
    --header "Content-Type: application/json" \
    --request POST \
    --data "{\"smiles\": \"$molecule"}" \
    -s
```
Resultado:

```
[1]
```


