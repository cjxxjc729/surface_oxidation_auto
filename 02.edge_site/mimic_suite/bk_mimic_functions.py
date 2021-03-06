#!/usr/bin/python3

import numpy as np
from ase.io import write
from ase.io import read
#from ase.geometry.analysis import Analysis
from ase import neighborlist

from ase import Atoms



def gen_e_vec(atoms):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  e_vec=[]
  if len(atoms)==1:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==2:
    Ria=atoms.get_distance(0,1,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2)
    ei3=np.cross(ei1,ei2)
  else:
    Ria=atoms.get_distance(0,1,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,2,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  e_vec=[ei1,ei2,ei3]
  return e_vec




def gen_M_atoms_nl(atoms,set_cutoff):
  #this funcitonal will generate lists of the neighbor atoms (within radii of set_cutoff) for all the atoms in file atoms. And the list will be sorted.
  cutOff = neighborlist.natural_cutoffs(atoms)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2   
  nl = neighborlist.NeighborList(cutOff, skin=0)
  nl.update(atoms)
  M_atoms_nl={}
  atm_nl={}
  for atm_id in range(len(atoms)):
    indices, offsets = nl.get_neighbors(i)
    atm_nl[atm_id]=sort_the_indices(atoms,i,indices)
    neibor_id_list_of_atm_id=atm_nl[atm_id][:,0]
    ll=[]
    for i in range(len(neibor_id_list_of_atm_id)):
      ll.append(neibor_id_list_of_atm_id[i].astype('int64'))
    M_atoms_nl[atm_id]=atoms[ll]
  return atm_nl,M_atoms_nl


def sort_the_indices(atoms,i,indices):
  A=[]
  for index in indices:
    dis=atoms.get_distance(i,index,mic=True)
    A=np.append(A,[index,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  for i in range(len(A)):
    A[i,0]=int(A[i,0])
  neibor_indices_sort_by_dis=A
  return neibor_indices_sort_by_dis

def gen_e_vec(atoms):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  e_vec=[]
  if len(atoms)==1:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==2:
    Ria=atoms.get_distance(0,1,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2)
    ei3=np.cross(ei1,ei2)
  else:
    Ria=atoms.get_distance(0,1,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,2,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  e_vec=[ei1,ei2,ei3]
  return e_vec



def gen_Dij(atoms):
  e_vec=gen_e_vec(atoms)
  vec_basis=e_vec
  Dij=[]
  symbols=[]
  for atm_id in range(len(atoms)):
      
    Rij=atoms.get_distance(0,atm_id,mic=True)
    Rij_vec=atoms.get_distance(0,atm_id,mic=True,vector=True)


    Rij_vec=np.mat(Rij_vec)
    #how to measure cor?
    #[ei1,ei2,ei3]*[xij;yij;zij]=[Rij_vec]'  3x3  3x1  3x1
    #=>[xij;yij;zij]=[x_vec,y_vec,z_vec]^-1*[Rij_vec]'
    xij_yij_zij=np.linalg.inv(np.transpose(vec_basis))*np.transpose(Rij_vec)
    xij_yij_zij=np.array(xij_yij_zij)
    xij= xij_yij_zij[0][0]
    yij= xij_yij_zij[1][0]
    zij= xij_yij_zij[2][0]
    Dij=np.append(Dij,[1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2)])
    symbols.append(atoms.get_chemical_symbols()[atm_id])
  Dij=Dij.reshape(-1,4)
  return Dij,symbols,e_vec



def score_the_similarity(Dija,symbolsa,Dijb,symbolsb):
  if len(Dija)!=len(Dijb):
    print("length of Dij is not match, skip")
    score=100
  else:
    speciesa=np.unique(symbolsa)
    speciesb=np.unique(symbolsb)    
    if len(speciesa)!=len(speciesb):
      print("specie number dont match, skip")
      score=100
    else:
      for i in range(len(speciesa)):
        if speciesa[i]!=speciesb[i]:
          print("specie kind dont match, skip")
          score=100
          break
      Dija_ele={}
      Dijb_ele={}
      score=0
      for specie in speciesa:
        Dija_ele[specie]=[]
        Dijb_ele[specie]=[]
        if len(symbolsa)!=len(symbolsb):
          print("specie number dont match,skip")
          score=100
          break
        for i in range(len(symbolsa)):
          index=symbolsa.index(specie,i)
          Dija_ele[specie]=np.append(Dija_ele[specie],Dija[index])
          index=symbolsb.index(specie,i)
          Dijb_ele[specie]=np.append(Dijb_ele[specie],Dijb[index])
        diff_Dijab=Dijb_ele[specie]-Dija_ele[specie]
        print("Dija_ele[specie="+specie)
        print(Dija_ele[specie])
        print("Dijb_ele[specie="+specie)
        print(Dijb_ele[specie])
        score_added=np.dot(diff_Dijab,diff_Dijab)
        print("score_added")
        print(score_added)
        score=score+score_added



  return score





