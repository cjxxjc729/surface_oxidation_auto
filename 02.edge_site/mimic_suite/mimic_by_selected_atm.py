#!/public1/soft/python/3.7.3/bin/python

filenamea=input("enter the originial filename: ")
filenameb=input("enter the goal filename: ")
#filenamea="NO2-minus.cif"
#filenameb="NO2.cif"
ia=input("enter the number of the original atom (ia, no zero): ")
ib=input("enter the number of the goal atom (ib,no zero): ")
ia=eval(ia)-1
ib=eval(ib)-1

from ase.io import write
from ase.io import read
from ase import neighborlist
from ase import geometry 
import numpy as np
from make_individual_cor import *

#---------------------------------begein of atom a-------------------------

if filenamea.split('.')[-2]=='pw':
  atomsa=read(filenamea,format="espresso-in")
else:
  atomsa=read(filenamea)

if filenameb.split('.')[-2]=='pw':
  atomsb=read(filenameb,format="espresso-in")
else:
  atomsb=read(filenameb)
print("atomic number for originial structure is :")
print(len(atomsa))
print("atomic number for goal structure is :")
print(len(atomsb))

set_cutoff=5
atm_nla=gen_atm_nl(atomsa,set_cutoff)
x_veca,y_veca,z_veca,atm_cla=gen_xyz_vec_and_atm_cl(atomsa,atm_nla)
print("atomic coordination genrates (x_vec y_vec z_vec) for molecule a")
print("selected a atom information:")
print(atomsa.get_chemical_symbols()[ia])
Dija=gen_Dij(atomsa,atm_nla,x_veca,y_veca,z_veca)
print("the number of the neighoring a atom selected")
print(len(Dija))
#-------------------------------end of atom a------------------------------
#-------------------------------begin of atom b------------------------------
atm_nlb=gen_atm_nl(atomsb,set_cutoff)
x_vecb,y_vecb,z_vecb,atm_cla=gen_xyz_vec_and_atm_cl(atomsb,atm_nlb)
print("atomic coordination genrates (x_vec y_vec z_vec) for molecule b")
print("selected b atom information:")
print(atomsb.get_chemical_symbols()[ib])
Dijb=gen_Dij(atomsb,atm_nlb,x_vecb,y_vecb,z_vecb)
print("the number of the neighoring b atom selected")
print(len(Dijb))
#same_index,diff_index,Dijb_can_be_used=same_and_diff(ia,ib,Dija,Dijb,0.03)
#print("same_index")
#print(same_index)

#debug
#print("Dija")
#print(Dija[ia])
#print("Dijb")
#print(Dijb[ib])
#same_index,diff_index,Dijb_can_be_used=same_and_diff(ia,ib,Dija,Dijb,1)
#print("same_index")
#print(same_index)
#debug

#-------------------------------end of atom b------------------------------
#------------------------------print properties------------------------
#same_and_diff(ia,ib,Dija,Dijb,0.03)


#------------------------------end print properties-------------------
#-----------------------------begin merge Dijb and Dija -------------
#Dija_merged=merge_Dij(Dija,Dijb,atomsa,atomsb)
Dija_merged=merge_Dij_by_known_id(Dija,Dijb,ia,ib)
#print(len(Dijb_merged[3]))
#print(len(Dijb[3]))
#print(len(Dija[15]))
#-----------------------------end merge Dijb and Dija-----------------
extraposa=gen_extrapos_from_Dij_assby_atoma(Dija_merged,atomsa)
#extraposb=gen_extrapos_from_Dij_assby_atoma(Dijb,atomsb)
print("delete possible repeated stoms")
print("atom unmber before deleted")
print(len(extraposa.get_positions()))
#write('atomsb_before.xyz', extraposb)
print("atom unmber after deleted")
geometry.get_duplicate_atoms(extraposa, cutoff=0.9, delete=True)
print(len(extraposa.get_positions()))

write('atomsb.cif', extraposa)


