#!/bin/bash

mkdir -p models/u_0.001_s_0.1
python generate_models.py data/cosmic_signatures_probabilities.txt data/Table_SignatureContribution__SupplTab21.csv 0.001 0.1 models/u_0.001_s_0.1
