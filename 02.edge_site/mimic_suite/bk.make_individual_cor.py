#!/public1/soft/python/3.7.3/bin/python

import numpy as np
from ase.io import write
from ase.io import read
#from ase.geometry.analysis import Analysis
from ase import neighborlist

from ase import Atoms

#this module will generate the individual coordinate system for each atom. The corrdinate will be built following the role written in deepmd.

def optimiz_cutoff(atoms,number_nl_treshhold):
  cutOff = neighborlist.natural_cutoffs(atoms)
  min_len_indices=1
  a=1
  while (min_len_indices<number_nl_treshhold): #number_nl_treshhold is 2 for calculating x y z vectors
    #print("circle" + str(a))
    nl = neighborlist.NeighborList(cutOff, self_interaction=False,bothways=True)
    nl.update(atoms)
    for i in range(len(atoms)):
      indices, offsets = nl.get_neighbors(i)
      if i==1:
        min_len_indices=len(indices)
      if len(indices)<min_len_indices:
        min_len_indices=len(indices)
      #print(str(indices))
    cutoff=[]
    for j in cutOff:
      cutoff.append(j+0.02)
    cutOff=cutoff
    a=a+1
    added_times=a-1
    added_value=0.02*added_times
  return (cutOff,added_value)

def gen_xyz_vec_and_atm_cl(atoms):
  cutOff,added_value = optimiz_cutoff(atoms,2)
  #use function optimiz_cutoff to get the cutOff
  #cutOff = neighborlist.natural_cutoffs(atoms)
  #--------------------------------------------------------------------------
  #for i in range(len(cutOff)):
  #  cutOff[i]=universe_cutoff_set
  nl = neighborlist.NeighborList(cutOff, self_interaction=False,bothways=True)
  nl.update(atoms)
  #here the nl(neighborlist) is acheved. use nl.get_neighbors(i) to get the atom inside the cutOff


  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  #now begins the vector maker, the first is x_vec  
  x_vec=[]
  atm_nls=[]
  #atm_nls: atoms's neibour list. "s" stands for stort list. Because only 1st and 2n nei is contained. the order is: atomi, atomi's 1st neiindex, atomi's 2nd neiindex
  for i in range(len(atoms)):   #each i represents an atom i.
    nei_indices, offsets = nl.get_neighbors(i)  #nei_indices is the neighbour indices for atom i
    #calculate the neiindex_dis
    neiindex_dis=[]
    for neiindex in nei_indices:
      dis=atoms.get_distance(i,neiindex)
      neiindex_dis.append([neiindex,dis])
    neiindex_dis=np.array(neiindex_dis)

    first_neiindex=neiindex_dis[np.argsort(neiindex_dis[:,1])][0,0]
    first_neiindex=int(first_neiindex)
    second_neiindex=neiindex_dis[np.argsort(neiindex_dis[:,1])][1,0]
    second_neiindex=int(second_neiindex)

    atm_nls.append([i,first_neiindex,second_neiindex])    
    x_vec.append(atoms.get_distance(i,first_neiindex,vector=True))

  x_vec=np.array(x_vec)
  atm_nls=np.array(atm_nls)
  #-----------------------end of generation of x_vec and nei_cl----------------------------------

  #-----------------------begin generation of y_vec_dot------------------------------------
  # the following is to generate y_vec_dot. Then use cross multiply to get z_vec. And finally use z_vec and x_vec to get y_vec.  y_vec_dot is the vect connect first nei and (first nei)^2, it can be find from x_vec. It is important to note, for some atoms (first nei)^2 will return to its origin. This case should be expelled, and we should use first nei's 2nd nei as the replacement of (first nei)^2. 
  y_vec_dot=[]
  atm_cl=[]  #atomic connection list.
  for i in range(len(atoms)):
    index=atm_nls[i,1]  #index is the goal 1st neighor. we here read it from 2nd line in atom_nl list.  It is the start point for y_vec_dot
    j=atm_nls[index,1]  #find what 1st neighbor's 1st neighbor in list. It is the end point for the y_vec_dot
    if j==i:    #but it gets an expectional. If j==i. we should then use the 2nd nei of index as end point
      j=atm_nls[index,2] 
    #print(atoms.get_chemical_symbols()[i]+str(i+1)+'-->'+atoms.get_chemical_symbols()[index]+str(index+1)+'-->'+atoms.get_chemical_symbols()[j]+str(j+1))    #used for check the resutls
    y_vec_dot.append(atoms.get_distance(index,j,vector=True))
    atm_cl.append([i,index,j])
  y_vec_dot=np.array(y_vec_dot)
  #-------------------------end generation of y_vec_dot------------------------------------
  #-------------------------begin generation of z_vec and y_vec----------------------------
  #1. normalize x_vec. z_vec=cross(x_vec.norm,y_vec).  y_vec=cross(x_vec.norm,z_vec.norm)
  z_vec=[]
  y_vec=[]
  x_vec_norm=[]
  i1=0
  for x_veci in x_vec:
    x_veci=x_veci/np.linalg.norm(x_veci)
    x_vec_norm.append(x_veci)
    z_veci=np.cross(x_veci,y_vec_dot[i1])  #get z_vec from x_vec and y_vec_dot
    z_veci=z_veci/np.linalg.norm(z_veci)  
    z_vec.append(z_veci)
    y_veci=np.cross(x_veci,z_veci)  #get z_vec from x_vec and z_vec
    y_vec.append(y_veci) 
    i1=i1+1
  x_vec=np.array(x_vec_norm)  #get x_vec  
  y_vec=np.array(y_vec)
  z_vec=np.array(z_vec)
  atm_cl=np.array(atm_cl) 
  return (x_vec,y_vec,z_vec,atm_cl)

def gen_atm_nl(atoms,set_cutoff):
  cutOff = neighborlist.natural_cutoffs(atoms)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2
  nl = neighborlist.NeighborList(cutOff, self_interaction=False,bothways=True)
  nl.update(atoms)
  atm_nl={}
  for i in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    atm_nl[i]=indices
  return atm_nl

def gen_Dij(atoms,atm_nl,x_vec,y_vec,z_vec):
  Dij={}
  #Dijele={}
  for i in range(len(atoms)):
    indices=atm_nl[i]  
    Dij[i]=[]
    #vec_basis_i=[x_vec[i],y_vec[i],z_vec[i]]
    #Dijele[i]=[]
    for j in indices:
      Rij=atoms.get_distance(i,j)
      Rij_vec=atoms.get_distance(i,j,vector=True)
      vec_basis_j=[x_vec[j],y_vec[j],z_vec[j]]
      #if (i==15 and j==30) or (i==3 and j==12):
      #  print(vec_basis_j)
      Rij_vec=np.mat(Rij_vec)
      #how to measure cor?
      #[x_vec,y_vec,z_vec]*[xij;yij;zij]=[Rij_vec]'  3x3  3x1  3x1
      #=>[xij;yij;zij]=[x_vec,y_vec,z_vec]^-1*[Rij_vec]'
      xij_yij_zij=np.linalg.inv(vec_basis_j)*np.transpose(Rij_vec)
      xij_yij_zij=np.array(xij_yij_zij)
      xij= xij_yij_zij[0][0]
      yij= xij_yij_zij[1][0]
      zij= xij_yij_zij[2][0]
      Dij[i].append([[1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2),j],atoms.get_chemical_symbols()[j]])
      #Dijele[i].append(atoms.get_chemical_symbols()[j])
    Dij[i]=np.array(Dij[i])
    indices2=np.argsort(Dij[i][:,0])
    Dij[i]=Dij[i][indices2]
  return Dij

def gen_extrapos_from_Dij_assby_atoma(Dija,x_veca,y_veca,z_veca,atomsa):
  extra_positions=[]
  extra_positions_ele=''
  for i in range(len(Dija)):
    pos_ref=atomsa.get_positions()[i]
    for j in range(len(Dija[i])):
      xija=Dija[i][j,0][1]/(Dija[i][j,0][0]**2)
      yija=Dija[i][j,0][2]/(Dija[i][j,0][0]**2)
      zija=Dija[i][j,0][3]/(Dija[i][j,0][0]**2)  
      #typical example: Dija[i][j,0]=[0.19027283845980128, 0.11513382090596258, -0.09355007556563862, -0.11914839361404195, 6]
      #typical example2: Dija[i][j,1]=H
      vect=xija*x_veca[i]+yija*y_veca[i]+zija*z_veca[i]
      extra_position=pos_ref+vect
      extra_positions.append((extra_position[0],extra_position[1],extra_position[2]))
      extra_positions_ele+=Dija[i][j,1]
  extrapos=Atoms(extra_positions_ele,extra_positions) 
  return extrapos


#def squeeze_positions(positions,)


def generate_neighbor_list_from_a_settled_cut(atoms,x_vec,y_vec,z_vec,added_value):
  cutOff = neighborlist.natural_cutoffs(atoms)
  cutoff=[]
  for j in cutOff:
    cutoff.append(j+added_value)
  cutOff=cutoff
  nl = neighborlist.NeighborList(cutOff, self_interaction=False,bothways=True)
  nl.update(atoms)
  Dij={}
  for i in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    Dij[i]=[]
    for j in indices:
      Rij=atoms.get_distance(i,j)
      Rij_vec=atoms.get_distance(i,j,vector=True)
      vec_basis_i=[x_vec[i],y_vec[i],z_vec[i]]
      Rij_vec=np.mat(Rij_vec)
      #how to measure cor?
      #[x_vec,y_vec,z_vec]*[x;y;z]=[Rij_vec]'  3x3  3x1  3x1
      #=>[x;y;z]=[x_vec,y_vec,z_vec]-1*[Rij_vec]'
      xij_yij_zij=np.linalg.inv(vec_basis_i)*np.transpose(Rij_vec)
      xij_yij_zij=np.array(xij_yij_zij)
      xij= xij_yij_zij[0][0]
      yij= xij_yij_zij[1][0]
      zij= xij_yij_zij[2][0]
      Dij[i].append([1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2),j])
    Dij[i]=np.array(Dij[i])
    indices2=np.argsort(Dij[i][:,0])
    Dij[i]=Dij[i][indices2]
  return Dij

def same_and_diff(ia,ib,Dija,Dijb,atomsa,atomsb,tol):  #i_from_a,i_from_b,Dij_from_a,Dij_from_b,atomsa,atomsb,tolerance
  same_index=[]
  diff_index=[] 
  ja=1
  for Da_ja in Dija[ia]:
    jb=1
    for Db_jb in Dijb[ib]:
      #print(atomsa.get_chemical_symbols()[int(Da_ja[4]]))
      #print(np.linalg.norm(Da_ja[0:4]-Db_jb[0:4]))
      if atomsa.get_chemical_symbols()[int(Da_ja[4])] == atomsb.get_chemical_symbols()[int(Db_jb[4])] and np.linalg.norm(Da_ja[0:4]-Db_jb[0:4])<tol:
        same_index.append([0.5*Da_ja[0]+0.5*Db_jb[0],int(Da_ja[4]),int(Db_jb[4])])
      else:
        diff_index.append([0.5*Da_ja[0]+0.5*Db_jb[0],int(Da_ja[4]),int(Db_jb[4])])
      jb=jb+1
    ja=ja+1
  return (same_index,diff_index)

