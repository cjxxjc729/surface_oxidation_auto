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
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  nl.update(atoms)
  M_atoms_nl={}
  atm_nl={}
  for atm_id in range(len(atoms)):
    indices, offsets = nl.get_neighbors(atm_id)
    atm_nl[atm_id]=sort_the_indices(atoms,atm_id,indices)
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
  if len(atoms)==1 or len(atoms)==2:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==3:
    Ria=atoms.get_distance(0,2,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  else:
    Ria=atoms.get_distance(0,2,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,3,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  e_vec=[ei1,ei2,ei3]
  return e_vec


def gen_e_vec_beta(atoms):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  e_vec=[]
  if len(atoms)==1 or len(atoms)==2:
    ei1=[1,0,0]
    ei2=[0,1,0]
    ei3=[0,0,1]
  elif len(atoms)==3:
    Ria=atoms.get_distance(0,2,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
    ei1=Ria/np.linalg.norm(Ria)
    Rib=[0,1,1]
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  else:
    Ria=atoms.get_distance(0,3,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,2,mic=True,vector=True)
    ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
    ei2=ei2s/np.linalg.norm(ei2s)
    ei3=np.cross(ei1,ei2)
  e_vec=[ei1,ei2,ei3]
  return e_vec




def gen_Dij_dot(atoms):
  e_vec=gen_e_vec(atoms)
  vec_basis=e_vec
  Dij=[]
  symbols=[]
  for atm_id in range(len(atoms)):
    if atm_id<=2:
      Dij=np.append(Dij,[0,0,0,0])
      symbols.append(atoms.get_chemical_symbols()[atm_id])
    else:
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
  #if len(Dij)>=1:
  Dij=Dij.reshape(-1,4)
  return Dij,symbols,e_vec


def gen_Dij(atoms):
  e_vec=gen_e_vec_beta(atoms)
  vec_basis=e_vec
  Dij=[]
  symbols=[]
  for atm_id in range(len(atoms)):
    if atm_id<=2:
      Dij=np.append(Dij,[1,0,0,0])
      symbols.append(atoms.get_chemical_symbols()[atm_id])
    else:
      #Rij=atoms.get_distance(0,atm_id,mic=True)
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
      Dij=np.append(Dij,[1,xij,yij,zij])
      symbols.append(atoms.get_chemical_symbols()[atm_id])
  #if len(Dij)>=1:
  Dij=Dij.reshape(-1,4)
  return Dij,symbols,e_vec




def score_the_similarity_dot(Dija,symbolsa,Dijb,symbolsb):

  if len(Dija)!=len(Dijb) or symbolsa[0]!=symbolsb[0] or len(Dija)<=5 or len(Dijb)<=5 :
    #print("length of Dij is not match or central element dont match, or neibour list too shor <2 , skip")
    score=100
  else:
    speciesa=np.unique(symbolsa)
    speciesb=np.unique(symbolsb)    
    if len(speciesa)!=len(speciesb):
      #print("specie number dont match, skip")
      score=100
    else:
      for i in range(len(speciesa)):
        if speciesa[i]!=speciesb[i]:
          #print("specie kind dont match, skip")
          score=100
          break
      Dija_ele={}
      Dijb_ele={}
      score=0
      break2=False
      for specie in speciesa:
        Dija_ele[specie]=[]
        Dijb_ele[specie]=[]
        for i in range(len(symbolsa)):
          #print(symbolsa)
          if symbolsa[i]==specie:
            Dija_ele[specie]=np.append(Dija_ele[specie],Dija[i])
          if symbolsb[i]==specie:
            Dijb_ele[specie]=np.append(Dijb_ele[specie],Dijb[i])
          if len(Dija_ele[specie])!=len(Dijb_ele[specie]):
            #print("specie number dont match,skip")
            score=100
            break2=True
            break
        if (break2):
          score=100
          break
        else:
          #print("Dija_ele[specie="+specie)
          #print(Dija_ele[specie])
          #print("Dijb_ele[specie="+specie)
          #print(Dijb_ele[specie])
          diff_Dijab=Dijb_ele[specie]-Dija_ele[specie]
          #print(diff_Dijab)
          score_added=np.dot(diff_Dijab,diff_Dijab)
          #print("score_added")
          #print(score_added)
          score=score+score_added

  return score

def score_the_similarity(Dija,symbolsa,Dijb,symbolsb):
  
  if len(Dija)!=len(Dijb) or symbolsa[0]!=symbolsb[0] or len(Dija)<=5 or len(Dijb)<=5 :
    #print("length of Dij is not match or central element dont match, or neibour list too shor <2 , skip")
    score=100
  else:
    speciesa=np.unique(symbolsa)
    speciesb=np.unique(symbolsb)
    if len(speciesa)!=len(speciesb):
      #print("specie number dont match, skip")
      score=100
    else:
      for i in range(len(speciesa)):
        if speciesa[i]!=speciesb[i]:
          #print("specie kind dont match, skip")
          score=100
          break
      Dija_ele={}
      Dijb_ele={}
      score=0
      break2=False
      scores=[]
      for specie in speciesa:
        Dija_ele[specie]=[]
        Dijb_ele[specie]=[]
        for i in range(len(symbolsa)):
          #print(symbolsa)
          if symbolsa[i]==specie:
            Dija_ele[specie]=np.append(Dija_ele[specie],Dija[i])
          if symbolsb[i]==specie:
            Dijb_ele[specie]=np.append(Dijb_ele[specie],Dijb[i])
          if len(Dija_ele[specie])!=len(Dijb_ele[specie]):
            #print("specie number dont match,skip")
            score=100
            break2=True
            break
        if (break2):
          score=100
          break
        else:
          #print("Dija_ele[specie="+specie)
          #print(Dija_ele[specie])
          #print("Dijb_ele[specie="+specie)
          #print(Dijb_ele[specie])
          Dija_ele[specie]=Dija_ele[specie].reshape(-1,4)
          Dijb_ele[specie]=Dijb_ele[specie].reshape(-1,4)
          A=[]
          for j in range(len(Dija_ele[specie])):
            dis=np.linalg.norm(Dija_ele[specie][j]-Dijb_ele[specie][j])
            dis=dis**2
            A=np.append(A,[j,dis])
          A=A.reshape(-1,2)
          A=A[np.argsort(A[:,1])]
          #print(specie)
          #print(A)
          
          score=A[len(A)-1,1]
          
          scores.append(score)
      if len(scores)>0:
        #print(scores)
        score=max(scores)
        #print(score)

  return score




def merge_based_on_individual_similar_part(atomsa,atomsa_nl,e_vec_a,atomsb,atomsb_nl,e_vec_b):
  impointb=atomsb_nl.get_positions()[0]
  impoint_ref=atomsa_nl.get_positions()[0]
  e_vec_ref=e_vec_a
  rotrance_atomsb=rotrance_atoms_to_ref(atomsb,impointb,impoint_ref,e_vec_b,e_vec_ref)
  new_atoms=merge_atoms(atomsa,rotrance_atomsb)

  return new_atoms


def merge_atoms(atomsa,atomsb):
  symbolesa=atomsa.get_chemical_symbols()
  symbolesb=atomsb.get_chemical_symbols()
  symboles=symbolesa+symbolesb

  positions=np.append(atomsa.get_positions(),atomsb.get_positions())
  positions=positions.reshape(-1,3)
  atoms=Atoms(symboles,positions,cell=atomsa.get_cell(),pbc=True)

  return atoms



def rotrance_atoms_to_ref(atoms,impoint,impoint_ref,e_vec,e_vec_ref):
  e_vec=np.mat(e_vec)
  e_vec_ref=np.mat(e_vec_ref)
  e_vec=e_vec.reshape(-1,3)
  e_vec_ref=e_vec_ref.reshape(-1,3)
  #e_vec=np.mat(e_vec)
  #e_vec_ref=np.mat(e_vec_ref)
  #base on:
  #(atoms.get_positions()-impoint)*np.linalg.inv(e_vec)*e_vac_ref+impoint_ref
  rotrance_cor=(atoms.get_positions()-impoint)*np.linalg.inv(e_vec)*e_vec_ref+impoint_ref
  symboles=atoms.get_chemical_symbols()
  rotrance_atoms=Atoms(symboles,rotrance_cor,cell=atoms.get_cell(),pbc=True)

  return rotrance_atoms



def cal_the_tranM(e_vec_a,e_vec_b):
  e_vec_a=e_vec_a.reshape(-1,3)
  e_vec_b=e_vec_b.reshape(-1,3)
  e_vec_a=np.mat(e_vec_a)
  e_vec_b=np.mat(e_vec_b)

  e_vec_a_inv=np.linalg.inv(e_vec_a)
  M_trans=e_vec_a_inv*e_vec_b

  return M_trans

def switch_symbol_to_HHeli(atoms):
  symbols=atoms.get_chemical_symbols()
  symbols_switch=[]
  for i in range(len(atoms)):
    if symbols[i]=='H':
      symbols_switch.append('Li')
    if symbols[i]=='O':
      symbols_switch.append('S')
    if symbols[i]=='C':
      symbols_switch.append('Si')
    if symbols[i]=='Co':
      symbols_switch.append('Fe')
    if symbols[i]=='N':
      symbols_switch.append('P')
  new_atoms=atoms
  new_atoms.set_chemical_symbols(symbols_switch)
  return new_atoms



