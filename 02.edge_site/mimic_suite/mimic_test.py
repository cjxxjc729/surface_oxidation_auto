#!/public1/soft/python/3.7.3/bin/python

#filename=input("enter the originial filename: ")
filenamea="NO2.pw.out"
filenameb="NO2-minus.cif"

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry 
import numpy as np
from make_individual_cor import *

#---------------------------------begein of atom a-------------------------
atomsa=read(filenamea)
atomsb=read(filenameb)
set_cutoff=5
x_vecb,y_vecb,z_vecb,atm_clb=gen_xyz_vec_and_atm_cl(atomsb)
print("atomic coordination genrates (x_vec y_vec z_vec) for molecule b")
atm_nlb=gen_atm_nl(atomsb,set_cutoff)
Dijb=gen_Dij(atomsb,atm_nlb,x_vecb,y_vecb,z_vecb)
#-------------------------------end of atom b------------------------------
#------------------------------print properties
#for i in range(1):
#  print("Dijb:"+str(i))
#  print(Dijb[i][2,1])

print("-----------------------")
print(atm_nlb)
#extra_positions,extra_positions_ele=gen_extrapos_from_Dij_assby_atoma(Dija,x_veca,y_veca,z_veca,atomsa)

#extrapos=gen_extrapos_from_Dij_assby_atoma(Dija,x_veca,y_veca,z_veca,atomsa)
extrapos=gen_extrapos_from_Dij_assby_atoma(Dijb,x_vecb,y_vecb,z_vecb,atomsb,atm_nlb)

print(len(extrapos.get_positions()))
#print(extra_positions)
#print(extra_positions_ele)  
geometry.get_duplicate_atoms(extrapos, cutoff=0.1, delete=True)

print(len(extrapos.get_positions()))

write('test.png', extrapos)



#  print(Dijbele[i])
#print("Dij generates. [i,1/Rij,xij/Rij^2,yij/Rij^2,zij/Rij^2,j] ")
#print(Dija[15])
#print("-------")
#print(Dijb[3])

#print("same indx")
#same_indics,diff_indics=same_and_diff(15,3,Dija,Dijb,atomsa,atomsb,0.1)
#print(same_indics)
#print("diff")
#print(diff_indics)

