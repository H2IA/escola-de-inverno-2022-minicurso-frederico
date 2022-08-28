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