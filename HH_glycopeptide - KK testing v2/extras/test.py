#!/usr/bin/python

from residue import Residue
from modification import Modification
from sequence import Sequence
from sequencespace import SequenceSpace

## Test residue module.
aa = Residue('Ala')
print(aa.name, aa.compo, aa.mass)

## Test modification module.
mod = Modification('Deamidated')
print(mod.name, mod.position, mod.number, mod.mass, mod.target)

user_mod = Modification('HexNAc', 3, 2, Residue('HexNAc').mass)
print(user_mod.name, user_mod.position, user_mod.number, user_mod.mass, user_mod.target)

## Test sequence module
str = 'DC(Carbamidomethyl)HLAQVPSHTVVAR'
sequence = Sequence(str)
## Check the sequence information.
print('Mass:', sequence.mass)
print('SEQ:\n')
for res in sequence.seq:
    print(res[0].name, res[0].mass)

print("B ions:\n")
for fragments in sequence.getFragments('B'):
    for fm in fragments:
        print(fm.getFragmentName(), "\t", fm.mass, "\n")
print("Y ions:\n")
for fragments in sequence.getFragments('Y'):
    for fm in fragments:
        print(fm.getFragmentName(), "\t", fm.mass, "\n")

#b_type = sequence.getFragments('B')

#print('B ions:\n')
#print(b_type)
#y_type = sequence.getFragments('Y')
#print('Y ions:\n')
#print(y_type)

str2 = 'DQYELLC(Carbamidomethyl)LDNTR'
new_seq = Sequence(str2)
new_mod = Modification('HexNAc', 9, 1, Residue('HexNAc').mass, 'Asn')
new_seq.appendModification(new_mod)

print('Mass:', new_seq.mass)
#b_type = new_seq.getFragments('B')
#print('B ions:\n')
#print(b_type)
#y_type = new_seq.getFragments('Y')
#print('Y ions:\n')
#print(y_type)
print("B ions:\n")
for fragments in new_seq.getFragments('B'):
    for fm in fragments:
        print(fm.getFragmentName(), "\t", fm.mass, "\n")
print("Y ions:\n")
for fragments in new_seq.getFragments('Y'):
    for fm in fragments:
        print(fm.getFragmentName(), "\t", fm.mass, "\n")

print(new_seq.getSequence())


## Test sequencespace module
seq_str = 'DTIC(Carbamidomethyl)IGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLC(Carbamidomethyl)R'
glycan_compo = {
    'dHex' : 0,
    'HexNAc' : 3,
    'Hex' : 4,
    'NeuAc' : 0,
    'Water' : 0}
num_sites = 2
glycan_sites = (9, 10, 22)

# 1 Deaminated
mod = Modification('Deamidated', -1, 1)
mod_list = [mod]
ss = SequenceSpace(seq_str, glycan_compo, glycan_sites, mod_list)
seq_list = ss.getTheoreticalSequence(num_sites)
for s in seq_list:
    print(s.getSequence())

# 1 Deaminated, 1 Carbamidomethyl
seq_str = 'DTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCR'
mod1 = Modification('Deamidated', -1, 1)
mod2 = Modification('Carbamidomethyl', -1, 1)
mod_list = [mod1, mod2]
ss = SequenceSpace(seq_str, glycan_compo, glycan_sites, mod_list)
seq_list = ss.getTheoreticalSequence(num_sites)
for s in seq_list:
    print(s.getSequence())



