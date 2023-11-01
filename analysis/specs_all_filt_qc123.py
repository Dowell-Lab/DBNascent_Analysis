#!/usr/bin/env python
# coding: utf-8

import itertools
import datetime
import pandas as pd
import numpy as np
import scipy.stats as sc
from sys import argv
import random
import math

## Define inputs/outputs
datafile = argv[1]
metafile = argv[2]
output = '/'.join(metafile.split('/')[0:-1]) + '/results/'

## Define functions
def readfile(infile):
  data = []
  with open(infile) as f:
    for line in f:
      inline = tuple(line.strip().split("\t"))
      data.append(inline)
  df = pd.DataFrame(data[1:], columns=data[0])
  return df

def find_current_samples(tissue, metadata):
  # Grab sample names for a specific tissue (and optionally disease state)
  if len(tissue.split("_")) > 1:
  #Disease state is encoded as boolean in metadata as 'tissue_diseasestate'
    samples = list(
      metadata[
        (metadata["tissue"] == tissue.split("_")[0]) &
        (metadata["disease"] == tissue.split("_")[1])
      ]["sample_name"]
    )
  else:
    samples = list(
      metadata[metadata["tissue"] == tissue]["sample_name"]
    )
  return samples

def tpm_stat_calc(category_list, tpm, metadata, outdir, label):
  # Define empty dfs
  tpm_mean = pd.DataFrame(columns=category_list, index=tpm.index)
  tpm_median = pd.DataFrame(columns=category_list, index=tpm.index)
  tpm_std = pd.DataFrame(columns=category_list, index=tpm.index)
  
  # Split out sample TPMs by tissue and add stats to dfs
  for tissue in category_list:
    current = find_current_samples(tissue, metadata)
    countdata = tpm[current]
    tpm_mean[tissue] = countdata.mean(axis=1)
    tpm_median[tissue] = countdata.median(axis=1)
    tpm_std[tissue] = countdata.std(axis=1)

  # Write out data
  tpm_mean.to_csv(outdir + label + '_tpm_mean.txt', sep='\t', index=True, header=True)
  tpm_median.to_csv(outdir + label + '_tpm_median.txt', sep='\t', index=True, header=True)
  tpm_std.to_csv(outdir + label + '_tpm_stdev.txt', sep='\t', index=True, header=True)

def specs_calc(category_list, tpm, metadata):
  # Initialize variables
  specs = pd.DataFrame(columns=category_list, index=tpm.index)
  weight = 1/(len(category_list)-1)
  nrna = len(tpm)

  # Run through loop for each tissue type (d)
  for tissue in category_list:
    # Find comparison tissues and samples for this tissue
    noncurrent = list(np.setdiff1d(category_list,tissue))
    curr_samps = find_current_samples(tissue, metadata)
    numd = len(curr_samps)
    # Initialize SPECS aggregator
    ptot = np.array([0]*nrna)

    # Compare with each tissue type not d (k)
    for comp_tiss in noncurrent:
      # Find sample TPMs for tissue k
      comp_samples = find_current_samples(comp_tiss, metadata)
      tpmk = tpm[comp_samples]
      numk = len(comp_samples)

      # Do pairwise comparison between all samps of d and k
      for samp in curr_samps:
        # Make df with TPM values for curr sample of d
        # Same size as tissue k sample TPM df for easy comparison
        tpmd = pd.concat([tpm[samp]]*numk, axis=1)
        tpmd.columns = tpmk.columns

        # Calculate single-sample SPECS component and add to aggregator
        dgreater = pd.DataFrame(np.where(tpmd > tpmk,1,0))
        I_kidj = dgreater.sum(axis=1)
        I_kidj_norm = I_kidj/(numk*numd)
        ptot = ptot + np.array(I_kidj_norm * weight)
    specs[tissue] = list(ptot)
  return specs

def specs_write(specs, outdir, label):
  specs.to_csv(outdir + label + '_specs_all.txt', sep='\t', index=True, header=True)
  specs_max = pd.concat([specs.max(axis=1),specs.idxmax(axis=1)], axis=1)
  specs_max = specs_max.sort_values(by=[0,1], axis=0, ascending=False)
  specs_max.to_csv(outdir + label + '_specs_maxval.txt', sep='\t', index=True, header=False)
  specs_min = pd.concat([specs.min(axis=1),specs.idxmin(axis=1)], axis=1)
  specs_min = specs_min.sort_values(by=[0,1], axis=0, ascending=False)
  specs_min.to_csv(outdir + label + '_specs_minval.txt', sep='\t', index=True, header=False)

## Read in data to dataframes
count_df = readfile(datafile)
metadata_all = readfile(metafile)

## Define rownames and extract TPM
# Make sure a feature is expressed in something
count_df.index = count_df["gene_transcript"]
tpm = count_df.iloc[:,6:].astype(float)
tpm = tpm.dropna(axis=1)
tpm = tpm[tpm.sum(axis=1) > 0]

## Determine unique tissues, with and without disease
# Filter to qc score 1-3
metadata_all = metadata_all[metadata_all["sample_name"].isin(list(tpm.columns))]
metadata_all = metadata_all[metadata_all["samp_qc_score"].astype(float) < 4]

# Find unique tissues regardless of disease w/more than 5 samples
# (excluding 'multiple' tissue - exclusively iPSCs)
tissues_unique = metadata_all.groupby(["tissue"])
tissues_sample_count = tissues_unique.size()[tissues_unique.size() > 5]
tissues_collapsed_temp = list(tissues_sample_count.index)
tissues_collapsed = [tiss for tiss in tissues_collapsed_temp if tiss != 'multiple']

# Find unique tissues including disease info w/more than 5 samples
# (excluding 'multiple' tissue - exclusively iPSCs)
tissues_unique = metadata_all.groupby(["tissue","disease"])
tissues_sample_count = tissues_unique.size()[tissues_unique.size() > 5]
tissues_bydisease = []
tissues_nondisease = []
for tissue in tissues_sample_count.index:
    if tissue != 'multiple':
        tissues_bydisease.append("_".join(tissue))
        if tissue[1] == "0":
            tissues_nondisease.append("_".join(tissue))

## Do nondisease tissue analysis
tpm_stat_calc(tissues_nondisease, tpm, metadata_all, output, 'filt_qc123_nondisease')
specs_nondisease = specs_calc(tissues_nondisease, tpm, metadata_all)
specs_write(specs_nondisease, output, 'filt_qc123_nondisease')

## Do bydisease tissue analysis
tpm_stat_calc(tissues_bydisease, tpm, metadata_all, output, 'filt_qc123_bydisease')
specs_bydisease = specs_calc(tissues_bydisease, tpm, metadata_all)
specs_write(specs_bydisease, output, 'filt_qc123_bydisease')

## Do collapsed tissue comparison (disease lumped with nondisease)
tpm_stat_calc(tissues_collapsed, tpm, metadata_all, output, 'filt_qc123_all')
specs_collapsed = specs_calc(tissues_collapsed, tpm, metadata_all)
specs_write(specs_collapsed, output, 'filt_qc123_all')
