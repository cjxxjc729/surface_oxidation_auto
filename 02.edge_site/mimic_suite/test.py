#!/usr/bin/python3

filename="blank.pw.in"

import numpy as np
from ase.io import write
from ase.io import read
from ase import neighborlist
from ase.neighborlist import NeighborList
from ase import Atoms
from ase import geometry

from make_individual_cor import *

atoms=read(filename,format="espresso-in")

atm_nl=gen_atm_nl(atoms,5)
print(atm_nl[38][:,0])

