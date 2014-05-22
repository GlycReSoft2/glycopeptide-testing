
#! /usr/bin/python

## Some of the codes are modified based on KK's scripts.

## Input:
## Arg 1 -- csv file containing glycopeptide composition
## Arg 2 -- list specifying the candidate glycosylation sites.
## Example: python HH_glycopeptide.py test_data.csv sample_glycosites.txt > report.txt

import sys
import csv
import re

from structure.residue import Residue
from structure.modification import Modification
from structure.sequencespace import SequenceSpace
from structure.sequence import Sequence
from structure.stub_glycopeptides import stubs
from structure.fragment import Fragment

def main(result_file, site_file):
    compo_dict = csv.DictReader(open(result_file, "r"), delimiter=",")
    site_list = [line.strip() for line in open(site_file, "r")]
    site_list = list(map(int, site_list))

    colnames = compo_dict.fieldnames
    match_pattern = re.compile(r'^G:(\w+)')

    glycan_identity = []
    for col in colnames:
        if match_pattern.match(col):
            glycan_identity.extend(match_pattern.findall(col))

    fragment_info=[]
# Each row as a dict.
for row in compo_dict:
    # For each row, read the information of glycan compound, peptide sequence,
    # ofGlycanAttachmentToPeptide, PeptideModification
    score = row['Score']
    precur_mass = row['MassSpec MW']
    theo_mass = row['Hypothesis MW']
    glycan_comp = row['Compound Key']
    seq_str = row['PeptideSequence']
    pep_mod = row['PeptideModification']
    pep_mis = row['PeptideMissedCleavage#']
    num_sites = int(row['#ofGlycanAttachmentToPeptide'])
    mass_error = row['PPM Error']
    volume = row['Total Volume']


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
        #print(seq_str, '\t', row['Compound Key'], mod_str, '\n', s.getSequence())
        seq_mod = s.getSequence()

        b_type = s.getFragments('B')
        b_ions = []
        b_ions_HexNAc = []
        all_b_ions = []
        for b in b_type:
            for fm in b:
                key = fm.getFragmentName()
                mass = fm.getMass()
                if "HexNAc" in key:
                    b_ions_HexNAc.append({"key":key, "mass":mass})
                    #all_b_ions.append({"key":key.split("+")[0], "mass":mass})
                else:
                    b_ions.append({"key":key, "mass":mass})
                    #all_b_ions.append({"key":key, "mass":mass})


        y_type = s.getFragments('Y')
        y_ions = []
        y_ions_HexNAc = []
        all_y_ions = []
        for y in y_type:
            for fm in y:
                key = fm.getFragmentName()
                mass = fm.getMass()
                if "HexNAc" in key:
                    y_ions_HexNAc.append({"key":key, "mass":mass})
                    #all_y_ions.append({"key":key.split("+")[0], "mass":mass})
                else:
                    y_ions.append({"key":key, "mass":mass})
                    #all_y_ions.append({"key":key, "mass":mass})


        pep_stubs = stubs(seq_str, pep_mod,num_sites,glycan_comp)
        stub_ions = pep_stubs.getStubs()
        oxonium_ions = pep_stubs.getOxoniums()
        #print('Y ions:\n')
        #print(y_type)

        fragment_info.append({"MS1_Score":score, "Obs_Mass":precur_mass, "Calc_mass":theo_mass, "ppm_error":mass_error,
        "Peptide":seq_str, "Peptide_mod":pep_mod, "Glycan":glycan_comp, "vol":volume, "glyco_sites":num_sites,
        "startAA":start_pos, "endAA":end_pos, "Seq_with_mod":seq_mod, "Glycopeptide_identifier":seq_mod + glycan_comp,
        "Oxonium_ions":oxonium_ions, "pep_stub_ions":stub_ions, "bare_b_ions":b_ions, "bare_y_ions":y_ions,
        "b_ions_with_HexNAc":b_ions_HexNAc,"y_ions_with_HexNAc":y_ions_HexNAc})


print "glycresoft1 output filename"
results = raw_input()
print "glycosylation sites filename"
sites = raw_input()

#site_list = [line.strip() for line in open(sys.argv[2], "r")]
site_list = [line.strip() for line in open(sites, "r")]
site_list = list(map(int, site_list))

#compo_dict = csv.DictReader(open(sys.argv[1], "r"), delimiter=",")
compo_dict = csv.DictReader(open(results, "r"), delimiter=",")

colnames = compo_dict.fieldnames
match_pattern = re.compile(r'^G:(\w+)')
glycan_identity = []
for x in colnames:
    if match_pattern.match(x):
        glycan_identity.extend(match_pattern.findall(x))

fragment_info=[]
# Each row as a dict.
for row in compo_dict:
    # For each row, read the information of glycan compound, peptide sequence,
    # ofGlycanAttachmentToPeptide, PeptideModification
    score = row['Score']
    precur_mass = row['MassSpec MW']
    theo_mass = row['Hypothesis MW']
    glycan_comp = row['Compound Key']
    seq_str = row['PeptideSequence']
    pep_mod = row['PeptideModification']
    pep_mis = row['PeptideMissedCleavage#']
    num_sites = int(row['#ofGlycanAttachmentToPeptide'])
    mass_error = row['PPM Error']
    volume = row['Total Volume']


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
        #print(seq_str, '\t', row['Compound Key'], mod_str, '\n', s.getSequence())
        seq_mod = s.getSequence()

        b_type = s.getFragments('B')
        b_ions = []
        b_ions_HexNAc = []
        all_b_ions = []
        for b in b_type:
            for fm in b:
                key = fm.getFragmentName()
                mass = fm.getMass()
                if "HexNAc" in key:
                    b_ions_HexNAc.append({"key":key, "mass":mass})
                    #all_b_ions.append({"key":key.split("+")[0], "mass":mass})
                else:
                    b_ions.append({"key":key, "mass":mass})
                    #all_b_ions.append({"key":key, "mass":mass})


        y_type = s.getFragments('Y')
        y_ions = []
        y_ions_HexNAc = []
        all_y_ions = []
        for y in y_type:
            for fm in y:
                key = fm.getFragmentName()
                mass = fm.getMass()
                if "HexNAc" in key:
                    y_ions_HexNAc.append({"key":key, "mass":mass})
                    #all_y_ions.append({"key":key.split("+")[0], "mass":mass})
                else:
                    y_ions.append({"key":key, "mass":mass})
                    #all_y_ions.append({"key":key, "mass":mass})


        pep_stubs = stubs(seq_str, pep_mod,num_sites,glycan_comp)
        stub_ions = pep_stubs.getStubs()
        oxonium_ions = pep_stubs.getOxoniums()
        #print('Y ions:\n')
        #print(y_type)

        fragment_info.append({"MS1_Score":score, "Obs_Mass":precur_mass, "Calc_mass":theo_mass, "ppm_error":mass_error, "Peptide":seq_str, "Peptide_mod":pep_mod, "Glycan":glycan_comp, "vol":volume, "glyco_sites":num_sites, "startAA":start_pos, "endAA":end_pos, "Seq_with_mod":seq_mod, "Glycopeptide_identifier":seq_mod + glycan_comp, "Oxonium_ions":oxonium_ions, "pep_stub_ions":stub_ions, "bare_b_ions":b_ions, "bare_y_ions":y_ions, "b_ions_with_HexNAc":b_ions_HexNAc,"y_ions_with_HexNAc":y_ions_HexNAc})

keys = ["MS1_Score","Obs_Mass","Calc_mass","ppm_error","Peptide","Peptide_mod","Glycan","vol","glyco_sites","startAA","endAA","Seq_with_mod","Glycopeptide_identifier","Oxonium_ions","pep_stub_ions","bare_b_ions", "bare_y_ions","b_ions_with_HexNAc","y_ions_with_HexNAc"]
print("enter output filename")
filename = raw_input()
f = open(filename,'wb')
dict_writer = csv.DictWriter(f,keys)
dict_writer.writer.writerow(keys)
dict_writer.writerows(fragment_info)
