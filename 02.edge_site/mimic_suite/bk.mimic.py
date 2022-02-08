#!/public1/soft/python/3.7.3/bin/python

#filename=input("enter the originial filename: ")
filenamea="NO2.pw.out"
filenameb="NO2-minus.cif"

from ase.io import write
from ase.io import read
#from ase.geometry.analysis import Analysis
from ase import neighborlist
import numpy as np

atomsa=read(filenamea)
atomsb=read(filenameb)

from make_individual_cor import *


x_vec,y_vec_dot,first_second_atom_index=gen_xyz_vec_and_atm_connect_list(atomsa)
x_vec_norm=[]
z_vec=[]
y_vec=[]
i=0
print("y_vec_dot1")
print(y_vec_dot[i])
for x_veci in x_vec:
  x_veci=x_veci/np.linalg.norm(x_veci)
  x_vec_norm.append(x_veci)
  z_veci=np.cross(x_veci,y_vec_dot[i])
  z_veci=z_veci/np.linalg.norm(z_veci)
  z_vec.append(z_veci)
  y_veci=np.cross(x_veci,z_veci)
  y_vec.append(y_veci)
  i=i+1
x_vec=x_vec_norm
#print(np.array(x_vec))
#print("--------------")
#print(np.array(y_vec))
#print("--------------")
#print(np.array(z_vec))
print("--------------")
print(np.array(first_second_atom_index))
print("atomic coordination genrates (x_vec y_vec z_vec) for molecule mimic into")
Dija,added_value=generate_neighbor_list(atomsa,5,x_vec,y_vec,z_vec)









print("added value to cut off is " + str(added_value))
x_vec,y_vec_dot,first_second_atom_index=gen_xyz_vec_and_atm_connect_list(atomsb)
i=0
x_vec_norm=[]
z_vec=[]
y_vec=[]

for x_veci in x_vec:
  x_veci=x_veci/np.linalg.norm(x_veci)
  x_vec_norm.append(x_veci)
  z_veci=np.cross(x_veci,y_vec_dot[i])
  z_veci=z_veci/np.linalg.norm(z_veci)
  z_vec.append(z_veci)
  y_veci=np.cross(x_veci,z_veci)
  y_vec.append(y_veci)
  i=i+1
x_vec=x_vec_norm
#print(np.array(x_vec))
#print("--------------")
#print(np.array(y_vec))
#print("--------------")
#print(np.array(z_vec))
print("--------------")
print(np.array(first_second_atom_index))
print("atomic coordination genrates (x_vec y_vec z_vec) for molecule mimic from")
#Dijb=generate_neighbor_list_from_a_settled_cut(atomsb,x_vec,y_vec,z_vec,added_value)
Dijb,added_value=generate_neighbor_list(atomsb,5,x_vec,y_vec,z_vec)

print("Dij generates. [i,1/Rij,xij/Rij^2,yij/Rij^2,zij/Rij^2,j] ")
print(np.around(Dija[15], decimals=4))

print(np.around(Dijb[3],decimals=4))

print("same indx")
same_indics,diff_indics=same_and_diff(15,3,Dija,Dijb,atomsa,atomsb,0.1)
print(same_indics)
print("diff")
print(diff_indics)

