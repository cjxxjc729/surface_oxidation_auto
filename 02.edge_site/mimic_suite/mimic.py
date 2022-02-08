#!/usr/bin/python3

#filenamea=input("enter the atom file part A: ")
#filenameb=input("enter the atom file part B : ")
filenamea="NO2-minus.cif"
filenameb="NO2-minus_rot.cif"

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry 
import numpy as np
from mimic_functions import *

#---------------------------------begin of atoma-------------------------

atomsa=read(filenamea)
atomsb=read(filenameb)
set_cutoff=5
atm_nla,M_atomsa_nl=gen_M_atoms_nl(atomsa,set_cutoff)
atm_nlb,M_atomsb_nl=gen_M_atoms_nl(atomsb,set_cutoff)

for atm_id_a in range(len(atomsa)):
  for atm_id_b in range(len(atomsb)):
    Dija,symbolsa,e_veca=gen_Dij(M_atomsa_nl[atm_id_a])
    Dijb,symbolsb,e_vecb=gen_Dij(M_atomsb_nl[atm_id_b])

    score=score_the_similarity(Dija,symbolsa,Dijb,symbolsb)
    print(str(atm_id_a)+" vs "+str(atm_id_b)+"score="+str(score))

