#!/usr/bin/python3

filename="blank.pw.in"
filename_ad_mol="CO2.cif"

id_ad_site_vac=[38,39]
id_ad_site_occ=[24]


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

if len(filename_ad_mol.split('.'))>=3 and filename_ad_mol.split('.')[-2]=='pw' and filename_ad_mol.split('.')[-1]=='in':
  atoms_ad_mol=read(filename_ad_mol,format="espresso-in")
else:
  atoms_ad_mol=read(filename_ad_mol)





