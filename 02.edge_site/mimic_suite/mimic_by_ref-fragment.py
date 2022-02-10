#!/usr/bin/python3

#filenamea=input("enter the atom file part A: ")
#filename_ref=input("enter the atom file part ref fragment : ")
#filenameb=input("enter the atom file part B: ")
filenamea="blank.cif"
filename_ref="site_adsorbate_only_biO.cif"
filenameb="site_adsorbate.cif"


#set_cutoffs=[3,4]
similarity_score_threshold=0.9
adsorbate_number=3

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

center_atm_id=get_center_atm_id(atoms_ref)

print("cutoff is: "+str(set_cutoff) )

if 1>0:
  for atm_id_a in range(len(atomsa)):
    #result=check_whether_surface_site(atoms_for_each_atoma[atm_id_a])
    #if result=="true":
      #print("this is surface site")
    score_a_ref,modea1,modea2=score_the_similarity(atoms_for_each_atoma[atm_id_a],atoms_for_each_atom_ref[center_atm_id])
    #print([score_a_ref,mode1,mode2])
    for atm_id_b in range(len(atomsb)):
      score_b_ref,modeb1,modeb2=score_the_similarity(atoms_for_each_atomb[atm_id_b],atoms_for_each_atom_ref[center_atm_id])
      #print([score_b_ref,mode1,mode2])
      if score_a_ref<similarity_score_threshold and score_b_ref<similarity_score_threshold:
        print(str(atm_id_a+1)+atomsa.get_chemical_symbols()[atm_id_a]+" vs "+str(atm_id_b+1)+atomsb.get_chemical_symbols()[atm_id_b]+": similarity score ="+str(score_a_ref)+" "+str(score_b_ref))
        print("modea1,odea2=")
        print([modea1,modea2])
        print("modeb1,modeb2=")
        print([modeb1,modeb2])
          #print(atomsb)
          #colored_atomsb=switch_symbol_to_HHeli(atomsb)
        merged_atoms=merge_based_on_individual_similar_part(atomsa,atoms_for_each_atoma[atm_id_a],atomsb,atoms_for_each_atomb[atm_id_b],modea1,modeb1,0.4)
        min_dises_for_adsorbate=[]
        for atm_id in range(len(merged_atoms)-adsorbate_number,len(merged_atoms)):
          min_dis=get_the_min_dis_in_atoms_for_atm_id(merged_atoms,atm_id)
          min_dises_for_adsorbate.append(min_dis)
        min_dis=min(min_dises_for_adsorbate)
        print(min_dis)
        if min_dis>=0.95:
          write("./merged_strs/new_atoms_merge"+str(atm_id_a+1)+"_"+str(atm_id_b+1)+"_cutoff_"+str(set_cutoff)+".cif",merged_atoms)
