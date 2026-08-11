#!/bin/bash
set -euo pipefail

PROCID=$1
NJOBS=$2

JOB=$((PROCID + 1))

CMSSW_DIR=/afs/cern.ch/user/c/charlesf/cmssw-dev/CMSSW_15_0_15/src
WORKDIR=${CMSSW_DIR}/Run3DimuonAnalysisTools/histograms

LIST=${CMSSW_DIR}/Run3DimuonAnalysisTools/Crab/all_mmgTree_files.txt
OUTDIR=/eos/user/c/charlesf/www/dark_photon/root

echo "Host: $(hostname)"
echo "Date: $(date)"
echo "PROCID=${PROCID}"
echo "JOB=${JOB}"
echo "NJOBS=${NJOBS}"
echo "LIST=${LIST}"
echo "OUTDIR=${OUTDIR}"

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd "${CMSSW_DIR}"
eval "$(scram runtime -sh)"

cd "${WORKDIR}"

mkdir -p "${OUTDIR}"

python3 fillDimuHistograms-binned.py \
  -l "${LIST}" \
  -o "${OUTDIR}/histos_" \
  -n "${NJOBS}" \
  -j "${JOB}"

echo "Finished job ${JOB}/${NJOBS} at $(date)"
