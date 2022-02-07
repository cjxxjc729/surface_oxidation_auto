#!/usr/bin/python3

#filename=input("enter the originial filename: ")
filename="blank.pw.in"


import numpy as np
from ase.io import write
from ase.io import read
import numpy as np


if len(filename.split('.'))>=3 and filename.split('.')[-2]=='pw' and filename.split('.')[-1]=='in':
  atoms=read(filename,format="espresso-in")
else:
  atoms=read(filename)
prefix=filename.split('.')[0]

def find_surf_atom(atoms)

positions=atoms.get_positions()
cor_center=np.mean(positions, axis=0)

  return atoms_surf
