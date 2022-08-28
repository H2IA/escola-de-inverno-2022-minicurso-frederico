#!/usr/bin/env bash

saved_path=$(bentoml get bbb-model:latest -o path)
cd $saved_path
gcloud builds submit --tag gcr.io/bbb-model-gcloud-run/bbb-classifier