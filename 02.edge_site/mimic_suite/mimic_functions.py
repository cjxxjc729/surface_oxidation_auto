#!/usr/bin/python3

import numpy as np
from ase.io import write
from ase.io import read
#from ase.geometry.analysis import Analysis
from ase.geometry import get_duplicate_atoms
from ase import neighborlist

from ase import Atoms


def gen_atoms_for_each_atom(atoms,set_cutoff):
  #this funcitonal will generate lists of the neighbor atoms (within radii of set_cutoff) for all the atoms in file atoms. And the list will be sorted.
  cutOff = neighborlist.natural_cutoffs(atoms)
  for i in range(len(cutOff)):
    cutOff[i]=set_cutoff/2   
  nl = neighborlist.NeighborList(cutOff, skin=0, bothways=True)
  nl.update(atoms)
  atoms_for_each_atom={}
  atm_nl={}
  for atm_id in range(len(atoms)):
    indices, offsets = nl.get_neighbors(atm_id)
    atm_nl[atm_id]=sort_the_indices(atoms,atm_id,indices)
    neibor_id_list_of_atm_id=atm_nl[atm_id][:,0]
    ll=[]
    for i in range(len(neibor_id_list_of_atm_id)):
      ll.append(neibor_id_list_of_atm_id[i].astype('int64'))
    atoms_for_each_atom[atm_id]=atoms[ll]
  return atoms_for_each_atom


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

def gen_e_vec(atoms,mode=0):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  mode_dependent_list=[[2,3],[2,4],[3,4],[2,5],[3,5],[4,5]]
  mode_dependent_list=np.array(mode_dependent_list)
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
    id1=min([mode_dependent_list[mode,0],len(atoms)-2])
    id2=min([mode_dependent_list[mode,1],len(atoms)-1])
    Ria=atoms.get_distance(0,id1,mic=True,vector=True)
    ei1=Ria/np.linalg.norm(Ria)
    Rib=atoms.get_distance(0,id2,mic=True,vector=True)
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


def gen_Dij(atoms,mode=0):
  e_vec=gen_e_vec(atoms,mode)
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

def score_the_similarity(atomsa,atomsb):
  mode_scores=[]
  for mode1 in range(3):
    for mode2 in range(3):
      Dija,symbolsa,e_vec_a=gen_Dij(atomsa,mode1)
      Dijb,symbolsb,e_vec_b=gen_Dij(atomsb,mode2)
      score=score_the_similarity_detailed(Dija,symbolsa,Dijb,symbolsb)
      mode_scores=np.append(mode_scores,[mode1,mode2,score])
  mode_scores=mode_scores.reshape(-1,3)
  mode_scores=mode_scores[np.argsort(mode_scores[:,2])]
  score=mode_scores[0,2]
  #print("mode_scores")
  #print(mode_scores)
  mode_of_a=mode_scores[0,0].astype('int64')
  mode_of_b=mode_scores[0,1].astype('int64')

  return score,mode_of_a,mode_of_b


def score_the_similarity_detailed(Dija,symbolsa,Dijb,symbolsb):

  if len(Dija)<len(Dijb) or symbolsa[0]!=symbolsb[0]:
    #print("length of Dij is not match or central element dont match , skip")
    score=100
  else:
    AA=[]
    for atmb_id in range(len(Dijb)):
      specieb=symbolsb[atmb_id]
      A=[]
      for atma_id in range(len(Dija)): 
        speciea=symbolsa[atma_id]
        if speciea==specieb:
          dis=np.linalg.norm(Dija[atma_id]-Dijb[atmb_id])
          dis=dis**2
          A=np.append(A,[atma_id,dis])
      if len(A)==0:
        smallest_dis_of_atma_ids_to_atmb_id=100
      else:
        A=A.reshape(-1,2)
        A=A[np.argsort(A[:,1])]
        smallest_dis_of_atma_ids_to_atmb_id=A[0,1]
      AA=np.append(AA,[atmb_id,smallest_dis_of_atma_ids_to_atmb_id])
    AA=AA.reshape(-1,2)
    AA=AA[np.argsort(AA[:,1])]
    score=AA[len(AA)-1,1]


  return score




def merge_based_on_individual_similar_part_no_del_atomsb_nl(atomsa,atomsa_nl,atomsb,atomsb_nl,mode_of_a=0,mode_of_b=0):
  #score,mode_of_a,mode_of_b=score_the_similarity(atomsa_nl,atomsb_nl)
  print("merging")
  print("mode_of_a,mode_of_b")
  print(mode_of_a,mode_of_b)
  e_vec_a=gen_e_vec(atomsa_nl,mode_of_a)
  e_vec_b=gen_e_vec(atomsb_nl,mode_of_b)

  impointb=atomsb_nl.get_positions()[0]
  impoint_ref=atomsa_nl.get_positions()[0]
  e_vec_ref=e_vec_a
  rotrance_atomsb=rotrance_atoms_to_ref(atomsb,impointb,impoint_ref,e_vec_b,e_vec_ref)
  new_atoms=merge_atoms(atomsa,rotrance_atomsb)

  return new_atoms

def merge_based_on_individual_similar_part(atomsa,atomsa_nl,atomsb,atomsb_nl,mode_of_a=0,mode_of_b=0,del_cutoff=0.1):
  #score,mode_of_a,mode_of_b=score_the_similarity(atomsa_nl,atomsb_nl)
  #print("merging")
  print("mode_of_a,mode_of_b")
  print(mode_of_a,mode_of_b)
  e_vec_a=gen_e_vec(atomsa_nl,mode_of_a)
  e_vec_b=gen_e_vec(atomsb_nl,mode_of_b)

  impointb=atomsb_nl.get_positions()[0]
  impoint_ref=atomsa_nl.get_positions()[0]
  e_vec_ref=e_vec_a
  rotrance_atomsb=rotrance_atoms_to_ref(atomsb,impointb,impoint_ref,e_vec_b,e_vec_ref)
  new_atoms=merge_atoms(atomsa,rotrance_atomsb)
  get_duplicate_atoms(new_atoms, del_cutoff, delete=True)

  return new_atoms


#def substract_atoms(atomsa,atomsb):
#  for atmb_id in range(len(atomsb)):
#    A=[]
#    for atma_id in range(len(atomsa)):
#      dis=np.linalg.norm(atomsa.get_positions()[atma_id]-atomsb.get_positions()[atmb_id])
#      A=np.append(A,[atma_id,dis])
#    A=A.reshape(-1,2)
#    A=A[np.argsort(A[:,1])]
#    del_atm_id=A[0,0]      


  return new_atoms

def merge_atoms(atomsa,atomsb):
  symbolesa=atomsa.get_chemical_symbols()
  symbolesb=atomsb.get_chemical_symbols()
  symboles=symbolesa+symbolesb

  positions=np.append(atomsa.get_positions(),atomsb.get_positions())
  positions=positions.reshape(-1,3)
  new_atoms=Atoms(symboles,positions,cell=atomsa.get_cell(),pbc=True)

  return new_atoms



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
      symbols_switch.append('Li')
    if symbols[i]=='C':
      symbols_switch.append('Li')
    if symbols[i]=='Co':
      symbols_switch.append('Li')
    if symbols[i]=='N':
      symbols_switch.append('Li')
    if symbols[i]=='Bi':
      symbols_switch.append('Li')
  positions=atoms.get_positions()
  new_atoms=Atoms(symbols_switch,positions,cell=atoms.get_cell(),pbc=True)
  #new_atoms.set_chemical_symbols(symbols_switch)
  return new_atoms

def get_the_max_dis_in_atoms(atoms):
  A=[]
  for atm_id in range(len(atoms)):
    for atm_jd in range(atm_id+1,len(atoms)):
      A.append(atoms.get_distance(atm_id,atm_jd,mic=True))
  max_dis=max(A)
  return max_dis

def get_the_min_dis_in_atoms(atoms):
  A=[]
  for atm_id in range(len(atoms)):
    for atm_jd in range(atm_id+1,len(atoms)):
      A.append(atoms.get_distance(atm_id,atm_jd,mic=True))
  min_dis=min(A)
  return min_dis

def get_the_min_dis_in_atoms_for_atm_id(atoms,atm_id):
  A=[]
  for atm_jd in range(len(atoms)):
    if atm_jd!=atm_id:
      A.append(atoms.get_distance(atm_id,atm_jd,mic=True))
  min_dis=min(A)
  return min_dis


def get_center_atm_id(atoms):
  center_cor=np.mean(atoms.get_positions())
  A=[]
  for atm_id in range(len(atoms)):
   dis=np.linalg.norm(atoms.get_positions()[atm_id]-center_cor)
   A=np.append(A,[atm_id,dis])
  A=A.reshape(-1,2)
  A=A[np.argsort(A[:,1])]
  center_atm_id=A[0,0]
  return center_atm_id


def check_whether_surface_site(atoms):
  #center_cor=np.mean(atoms.get_positions()) 
  char_vec=0
  for atm_id in range(2,len(atoms)):
    R_vec=atoms.get_distance(0,atm_id,mic=True,vector=True)
    char_vec=char_vec+R_vec
  char_vec=char_vec/(len(atoms)-2)
  print(char_vec)
  #char_vec=center_cor-atoms.get_positions()[0]
  result="true"
  for atm_id in range(2,len(atoms)):
    test_vec=atoms.get_positions()[atm_id]-atoms.get_positions()[0]
    cos_value=np.dot(char_vec,test_vec)/np.linalg.norm(char_vec)/np.linalg.norm(test_vec)
    if cos_value<-0.5:
      result="false"
      break


  return result
