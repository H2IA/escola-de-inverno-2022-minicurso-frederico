#!/usr/bin/env bash

set -e

# get molecule from CLI argument

if [ -z "$1" ]; then
    molecule='CCCC'
else
    molecule=$1
fi

# send request and parse JSON response

URL=$(terraform output -json | jq -r .endpoint.value)/classify
prediction=$(curl \
    --header "Content-Type: application/json" \
    --request POST \
    --data "{\"smiles\": \"$molecule\"}" \
    -s \
    $URL | jq '.[0]')

# show results

echo -e "molecule,+/-"

if [ "$prediction" -eq "1" ]; then
    echo -e "$molecule,BBB+"
else
    echo -e "$molecule,BBB-"
fi