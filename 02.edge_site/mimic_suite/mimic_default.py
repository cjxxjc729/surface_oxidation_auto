#!/usr/bin/python3

filenamea=input("enter the atom file part A: ")
filenameb=input("enter the atom file part B : ")
#filenamea="../blank.cif"
#filenameb="../CO2_ad.cif"

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry 
import numpy as np
from mimic_functions import *

#---------------------------------begin of atoma-------------------------

atomsa=read(filenamea)
atomsb=read(filenameb)
set_cutoffs=[3,4]
similarity_score_threshold=0.9

for set_cutoff in set_cutoffs:
  atm_nla,M_atomsa_nl=gen_M_atoms_nl(atomsa,set_cutoff)
  atm_nlb,M_atomsb_nl=gen_M_atoms_nl(atomsb,set_cutoff)
  print("for cutoff of: "+str(set_cutoff) )
  for atm_id_a in range(len(atomsa)):
    for atm_id_b in range(len(atomsb)):
      Dija,symbolsa,e_vec_a=gen_Dij(M_atomsa_nl[atm_id_a])
      Dijb,symbolsb,e_vec_b=gen_Dij(M_atomsb_nl[atm_id_b])
      score=score_the_similarity(Dija,symbolsa,Dijb,symbolsb)
      if score<similarity_score_threshold:
        print(str(atm_id_a+1)+atomsa.get_chemical_symbols()[atm_id_a]+" vs "+str(atm_id_b+1)+atomsb.get_chemical_symbols()[atm_id_b]+": similarity score ="+str(score))
        #clored_atomsb_nl=switch_symbol_to_HHeli(M_atomsb_nl[atm_id_b])
        #new_atoms=merge_based_on_individual_similar_part(atomsa,M_atomsa_nl[atm_id_a],e_vec_a,atomsb,clored_atomsb_nl,e_vec_b)
        new_atoms=merge_based_on_individual_similar_part(atomsa,M_atomsa_nl[atm_id_a],e_vec_a,atomsb,M_atomsb_nl[atm_id_b],e_vec_b)
        #clored_atomsb_nl=switch_symbol_to_HHeli(M_atomsb_nl[atm_id_b])
        #new_atoms=merge_based_on_individual_similar_part(M_atomsa_nl[atm_id_a],M_atomsa_nl[atm_id_a],e_vec_a,clored_atomsb_nl,clored_atomsb_nl,e_vec_b)
        write("./merged_strs/new_atoms_merge"+str(atm_id_a+1)+"_"+str(atm_id_b+1)+"_cutoff_"+str(set_cutoff)+".cif",new_atoms)
