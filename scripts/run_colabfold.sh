#!/bin/bash
set -euo pipefail

INPUTFILE="${1:?Provide a FASTA path}"
OUTDIR="${2:?Provide an output directory}"
NUMRECYCLE="${3:?Provide a number of recycles}"
NUMMODELS="${4:?Provide a numer of models}"
RANDOMSEED="${5:?Provide random seed}"

mkdir -p "$OUTDIR"

colabfold_batch \
  --num-recycle ${NUMRECYCLE} \
  --amber \
  --templates \
  --num-models ${NUMMODELS} \
  --random-seed ${RANDOMSEED} \
  ${INPUTFILE} \
  ${OUTDIR}
