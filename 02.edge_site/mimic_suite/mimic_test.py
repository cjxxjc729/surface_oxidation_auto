#!/usr/bin/python3

#filenamea=input("enter the atom file part A: ")
#filename_ref=input("enter the atom file part ref fragment : ")
#filenameb=input("enter the atom file part B: ")
filenamea="blank.cif"
filename_ref="site_adsorbate_only_biO.cif"
filenameb="site_adsorbate.cif"


#set_cutoffs=[3,4]
similarity_score_threshold=0.9

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry 
import numpy as np
from mimic_functions import *

#---------------------------------begin of pro-------------------------

atomsa=read(filenamea)
atoms_ref=read(filename_ref)
atomsb=read(filenameb)
max_dis_in_ref=get_the_max_dis_in_atoms(atoms_ref)
set_cutoff=max_dis_in_ref+0.1


atoms_for_each_atoma=gen_atoms_for_each_atom(atomsa,set_cutoff)
atoms_for_each_atom_ref=gen_atoms_for_each_atom(atoms_ref,set_cutoff)
atoms_for_each_atomb=gen_atoms_for_each_atom(atomsb,set_cutoff)


print("cutoff is: "+str(set_cutoff) )

for atm_id_b in range(len(atomsb)):
  print(atm_id_b)
  score_b_ref=score_the_similarity(atoms_for_each_atomb[atm_id_b],atoms_for_each_atom_ref[2])
  print(score_b_ref)
  new_atoms=merge_based_on_individual_similar_part(atomsb,atoms_for_each_atomb[atm_id_b],atoms_ref,atoms_for_each_atom_ref[2])
  write("./merged_strs/new_atoms_merge"+str(atm_id_b+1)+"_"+str(3)+"_cutoff_"+str(set_cutoff)+".cif",new_atoms)
