#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import zarr
import allel
import argparse
from sweep_tools import *

parser = argparse.ArgumentParser(description='Prepare files for Relate to infer whole-genome geneaologies')
parser.add_argument('--pops', type=str, action='store', default=['BFgam'], help='Populations to include')
parser.add_argument('--name', type=str, action='store', default='BFgam', help='Populations name')
parser.add_argument('--samples', type=str, action='store', help='Tab-delimited Samples metadata file')
parser.add_argument('--zarr', type=str, action='store', help='Path to Zarr files')
parser.add_argument('--chrom', type=str, action='store', default=['3L','3R'], help='Which chromosomes to use')
args=parser.parse_args()

### store arguments
chrom = args.chrom
samples = pd.read_csv(args.samples, sep="\t")
popname = args.name

print("-------------------- RelateAg - Sanjay C Nagi --------------------")

pops = [args.pops]
### read in hap array and positions
haps, pop_bool, sweep_region, pos = get_haplos(pops, chrom, 0, 55000000, samples, biallelic=False, zarrpath=args.zarr)

#count alleles and find bi_allelic SNPs and non-singletons
ac = haps.count_alleles()
print(ac.shape)
bial = ac.is_biallelic()
bial_ = np.logical_and(bial[:], ~ac.is_singleton()[:])
ac2 = haps.compress(bial_, axis=0)
print(ac2.shape)

#subset positions
pos_ = pos[bial_]
print("After filtering, there are", len(pos_), f"SNPs remaining for {pops} and {chrom}")

## write file with biallelic SNP positions
pd.Series(pos_).to_csv(f'data/{popname}_{chrom}_filtered_SNPs_pos', index=False)

print("Writing SNPS, .poplabels, ox_codes")

if len(pops[0]) == 1:
    df_samples = samples[samples.population == pops[0]]
else:
    df_samples = samples[samples.population.isin(pops)]

#write ox codes to csv
df_samples.ox_code.to_csv(f"data/{popname}_{chrom}_ox_codes", index=False)

#change male female coding and write .poplabels file
poplabels = df_samples.iloc[:,[0,2,2,9]]
poplabels.columns = ['sample', 'population', 'group', 'sex']
di = {"M":1, "F":2}
poplabels['sex'] = poplabels.sex.map(di)
poplabels['group'] = popname
poplabels.to_csv(f"data/{popname}_{chrom}.poplabels", index=False)

