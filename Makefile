setup:
	@conda env create || conda env update

dirs:
	@mkdir -p data/raw
	@mkdir -p data/features
	@mkdir -p data/models

download:
	@bash scripts/download_data.sh

features:
	@python features.py

train:
	@python train.py

build:
	@bentoml build

serve:
	@bentoml serve bbb-model:latest --production

container:
	@bentoml containerize bbb-model:latest --platform=linux/amd64

deploy:
	@bash scripts/deploy_to_gcp.sh

all: download features train build container deploy
