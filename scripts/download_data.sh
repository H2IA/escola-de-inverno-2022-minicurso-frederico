#!/usr/bin/env bash

set -e

wget \
    -O data/raw/b3db.tsv \
    'https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_classification.tsv'

    