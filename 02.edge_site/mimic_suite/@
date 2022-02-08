#!/usr/bin/python3

#filename=input("enter the originial filename: ")
filename="blank.pw.in"


import numpy as np
from ase.io import write
from ase.io import read
from ase import neighborlist
from ase.neighborlist import NeighborList
from ase import Atoms
from ase import geometry
#from make_box import *


if len(filename.split('.'))>=3 and filename.split('.')[-2]=='pw' and filename.split('.')[-1]=='in':
  atoms=read(filename,format="espresso-in")
else:
  atoms=read(filename)
prefix=filename.split('.')[0]

def find_surf_site(atoms,allowed_element,skin_value,cordi_tol):
  positions=atoms.get_positions()
  cor_center=np.mean(positions, axis=0)
  S=neighborlist.natural_cutoffs(atoms)
  nl=NeighborList(S,skin=skin_value,self_interaction=False)
  nl.update(atoms)
  #----------------------
  atm_id_surf_site=[]
  symbols_hyb=atoms.get_chemical_symbols()
  for atm_id in range(len(atoms)):
    if atoms.get_chemical_symbols()[atm_id] in allowed_element:
      indice, offsets=nl.get_neighbors(atm_id)
      cordi=len(indice)
      if cordi<=cordi_tol and atoms.get_positions()[atm_id][2]>cor_center[2]:
          atm_id_surf_site.append(atm_id)
          symbols_hyb[atm_id]='He'
          
  atoms_surf=atoms[atm_id_surf_site]
  atom_hybrid_surf=atoms      
  atom_hybrid_surf.set_chemical_symbols(symbols_hyb)

  return atoms_surf,atom_hybrid_surf

allowed_element=['O']
skin_value=1
cordi_tol=2

atoms_surf,atom_hybrid_surf=find_surf_site(atoms,allowed_element,skin_value,cordi_tol)

write("atoms_surf.cif",atoms_surf)
write("atoms_hybrid_surf.cif",atom_hybrid_surf)
