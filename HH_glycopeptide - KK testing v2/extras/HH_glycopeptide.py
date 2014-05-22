
#! /usr/bin/python

## Some of the codes are modified based on KK's scripts.

## Input:
## Arg 1 -- csv file containing glycopeptide composition
## Arg 2 -- list specifying the candidate glycosylation sites.
## Example: python HH_glycopeptide.py test_data.csv sample_glycosites.txt > report.txt

import sys
import csv
import re

from modification import Modification
from sequencespace import SequenceSpace

site_list = [line.strip() for line in open(sys.argv[2], "r")]
site_list = list(map(int, site_list))

compo_dict = csv.DictReader(open(sys.argv[1], "r"), delimiter=",") 

colnames = compo_dict.fieldnames
match_pattern = re.compile(r'^G:(\w+)')
glycan_identity = []
for x in colnames:
    if match_pattern.match(x):
        glycan_identity.extend(match_pattern.findall(x))

# Each row as a dict.
for row in compo_dict:
    # For each row, read the information of glycan compound, peptide sequence, 
    # ofGlycanAttachmentToPeptide, PeptideModification
    seq_str = row['PeptideSequence']
    num_sites = int(row['#ofGlycanAttachmentToPeptide'])
    if (seq_str == '') | (num_sites == 0): # Not interested
        continue

    # Working on the modification
    mod_str = row['PeptideModification']
    match_pattern = re.compile(r'(\d+)(\w+)')
    items = match_pattern.findall(mod_str)
    mod_list = []
    for i in items:
        mod = Modification(i[1], -1, int(i[0]))
        mod_list.append(mod)

    # Working on the glycosylation sites
    start_pos = int(row['StartAA'])
    end_pos = int(row['EndAA'])
    
    glycan_sites = set(site_list).intersection(range(start_pos, end_pos+1))
    
    if len(glycan_sites) == 0: # No recorded sites.
        continue

    ## Adjust the glycan_sites to relative position
    glycan_sites = [x - start_pos for x in glycan_sites]

    # Working on the glycan composition
    glycan_compo = {}
    for g in glycan_identity:
        glycan_compo[g] = int(row[''.join(['G:',g])])

    ss = SequenceSpace(seq_str, glycan_compo, glycan_sites, mod_list)
    seq_list = ss.getTheoreticalSequence(num_sites)
    for s in seq_list:
        print(seq_str, '\t', row['Compound Key'], mod_str, s.getSequence())
