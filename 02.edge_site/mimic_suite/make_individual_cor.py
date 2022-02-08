#!/usr/bin/python3

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

def gen_e_vec(atoms,atm_nl):
  #------------------------beginiing of generation of x_vec and nei_cl----------------------
  e_vec=[]
  for atm_id in range(len(atoms)):   #each i represents an atom i.
    neibor_id_list_of_atm_id=atm_nl[atm_id][:,0]  #indices is the list of sorted. if neibor_id_list_of_atm_id=[15,7,3,4] means for atm_id, the 1st neibor atom is 15, 2nd neibor atom is 7. 
    if len(neibor_id_list_of_atm_id)==1:
      ei1=[1,0,0]
      ei2=[0,1,0]
      ei3=[0,0,1]
    elif len(neibor_id_list_of_atm_id)==2:
      id_1st_nei=int(neibor_id_list_of_atm_id[1])
      Ria=atoms.get_distance(atm_id,id_1st_nei,mic=True,vector=True)   #notice neibor_id_list_of_atm_id[0] is itself !
      ei1=Ria/np.linalg.norm(Ria)
      Rib=[0,1,1]
      ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
      ei2=ei2s/np.linalg.norm(ei2)
      ei3=np.cross(ei1,ei2)
    else:
      id_1st_nei=int(neibor_id_list_of_atm_id[1])
      Ria=atoms.get_distance(atm_id,id_1st_nei,mic=True,vector=True)
      ei1=Ria/np.linalg.norm(Ria)
      id_2nd_nei=int(neibor_id_list_of_atm_id[2])
      Rib=atoms.get_distance(atm_id,id_2nd_nei,mic=True,vector=True)
      ei2s=Rib-np.dot(np.dot(Rib,ei1),ei1)
      ei2=ei2s/np.linalg.norm(ei2s)
      ei3=np.cross(ei1,ei2)
    e_vec=np.append(e_vec,ei1)
    e_vec=np.append(e_vec,ei2)
    e_vec=np.append(e_vec,ei3)
  e_vec=e_vec.reshape(-1,9)
  return e_vec



def gen_e_vec_by_atm_id0(atoms):
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
  cutOff = neighborlist.natural_cutoffs(atomis)
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
    M_atoms_nl[atm_id]=atoms[neibor_id_list_of_atm_id]

  #atm_nl[i] gives the sorted neighbor list of atm_id=i. 
  #for instance, we have
  #atm_nl[38]=
  #[[38.          0.        ]
  # [40.          3.06549498]
  # [43.          4.17021823]
  # [44.          4.41685152]
  # [39.          4.61817254]]
  # so atm_nl[atm_id][:,0] output the indices for neighbor.
  #    atm_nl[atm_id][:,1] output the associated distances

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


def gen_Dij_large(atoms,atm_nl,e_vec):
  Dij={}
  for atm_id in range(len(atoms)):
    neibor_id_list_of_atm_id=atm_nl[atm_id][:,0]  
    vec_basis=e_vec[atm_id].reshape(-1,3)

    Dij[atm_id]=[]
    j=0
    for atm_jd in neibor_id_list_of_atm_id[1:-1]:             
      atm_jd=int(atm_jd)
      Rij=atoms.get_distance(atm_id,atm_jd,mic=True)
      Rij_vec=atoms.get_distance(atm_id,atm_jd,mic=True,vector=True)


      Rij_vec=np.mat(Rij_vec)
      #how to measure cor?
      #[ei1,ei2,ei3]*[xij;yij;zij]=[Rij_vec]'  3x3  3x1  3x1
      #=>[xij;yij;zij]=[x_vec,y_vec,z_vec]^-1*[Rij_vec]'
      xij_yij_zij=np.linalg.inv(np.transpose(vec_basis))*np.transpose(Rij_vec)
      xij_yij_zij=np.array(xij_yij_zij)
      xij= xij_yij_zij[0][0]
      yij= xij_yij_zij[1][0]
      zij= xij_yij_zij[2][0]
      #Dij[atm_id].append([[1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2),atm_jd],atoms.get_chemical_symbols()[atm_jd],vec_basis])
      Dij[atm_id].append([[1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2),atm_jd]
      #Dijele[i].append(atoms.get_chemical_symbols()[j])
      j=j+1
    Dij[atm_id]=np.array(Dij[atm_id])
    #indices2=np.argsort(Dij[i][:,0])
    #Dij[i]=Dij[i][indices2]
  return Dij

def gen_Dij(atoms):
  e_vec=gen_e_vec_by_atm_id0(atoms)
  vec_basis=e_vec.reshape(-1,3)
  Dij=[]
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
    Dij=np.append(Dij,[[1/Rij,xij/(Rij**2),yij/(Rij**2),zij/(Rij**2),atoms.get_chemical_symbols()[atm_id]])
  Dij=Dij.reshape(-1,5)
  return Dij,e_vec



def score_the_similarity(Dija,Dijb):

  #

  return 






def merge_Dij(Dija,Dijb,atomsa,atomsb):
  #typical example: Dija[i][j][0]=[0.19027283845980128, 0.11513382090596258, -0.09355007556563862, -0.11914839361404195, 6]
  #typical example2: Dija[i][j][1]=H
  #typical example3: Dija[i][j][2]=[x_vec[i][j],y_vec[i][j],z_vec[i][j]]
  #the first is judgement. Whether can Dija[i] and Dijb[j] merge. The criteria is same_index fraction surpasses 0.3
  #The second is for unknown error, i have to create an extra Dija after merged, which is named as Dija_merged. Dija_merged has the same structure with Dija.  
  tol=0.03
  number_treshold=1
  how_many_selected=1   #
  ia_ib_matchscores=[]
  Dija_merged={}
  print("the following is the candidate")
  for ia in range(len(Dija)):
    Dija_merged[ia]=[]
    for i in range(len(Dija[ia])):
      Dija_merged[ia].append(Dija[ia][i])
    for ib in range(len(Dijb)):
      if atomsa.get_chemical_symbols()[ia]==atomsb.get_chemical_symbols()[ib]:
        same_index,diff_index,Dijb_can_be_merged=same_and_diff(ia,ib,Dija,Dijb,tol,2)
        if len(same_index)==0:
          fraction=0
        else:  
          fraction=len(np.unique(same_index[:,0]))/len(Dija[ia])
        if len(same_index)>number_treshold:
          print("the same fraction is "+str(len(same_index))+". merge Dija " + str(ia)+ " and Dijb " + str(ib),", with the fraction of " + str(fraction)+" different index is " + str(len(diff_index)))
          ia_ib_matchscores.append([ia,ib,len(same_index),fraction,len(diff_index)])
  ia_ib_matchscores=np.array(ia_ib_matchscores)
  ia_ib_matchscores=ia_ib_matchscores[np.argsort(ia_ib_matchscores[:,3])[::-1]]     #from highest fraction to the lowest
  print("ia_ib_matchscores matix is: ")
  print(ia_ib_matchscores)
  print("the following is what we chosed")
  for ia_ib_matchscore in ia_ib_matchscores[0:how_many_selected]:
    ia_chose=ia_ib_matchscore[0]
    ib_chose=ia_ib_matchscore[1]
    same_index_chose,diff_index_chose,Dijb_can_be_merged_chose=same_and_diff(ia_chose,ib_chose,Dija,Dijb,tol,2)
    print("same_index_chose")
    print(same_index_chose)
    print(same_index_chose[0][0])
    print(same_index_chose[0][1])
    print("the same number is "+str(int(ia_ib_matchscore[2]))+". merge Dija " + str(int(ia_ib_matchscore[0]))+ " and Dijb " + str(int(ia_ib_matchscore[1]))+", with the fraction of "+ str(ia_ib_matchscore[3]))
    x_vec_y_vec_z_vec_a=Dija[ia_chose][int(same_index_chose[0][0])][2]
    x_vec_y_vec_z_vec_b=Dijb[ib_chose][int(same_index_chose[0][1])][2]
    print("using x_vec_y_vec_z_vec of " + str(x_vec_y_vec_z_vec_a) +"and x_vec_y_vec_z_vec of "+ str(x_vec_y_vec_z_vec_b) + " to generate transform matrix")
    M_trans=cal_the_tranM(x_vec_y_vec_z_vec_b,x_vec_y_vec_z_vec_a)
    print("the M trans is")
    print(M_trans)
    for i in range(len(Dijb[ib_chose])):
      #typical example2: Dija[i][j][2]=[x_vec,y_vec,z_vec]
      x_vec_y_vec_z_vec=[Dijb[ib_chose][i][2][0],Dijb[ib_chose][i][2][1],Dijb[ib_chose][i][2][2]]
      x_vec_y_vec_z_vec=np.array(x_vec_y_vec_z_vec)
      x_vec_y_vec_z_vec=x_vec_y_vec_z_vec*M_trans
      x_vec_y_vec_z_vec=np.array(x_vec_y_vec_z_vec)
      print(i)
      
      Dijb[ib_chose][i][2][0]=x_vec_y_vec_z_vec[0]
      Dijb[ib_chose][i][2][1]=x_vec_y_vec_z_vec[1]
      Dijb[ib_chose][i][2][2]=x_vec_y_vec_z_vec[2]
      Dija_merged[ia_chose].append(Dijb[ib_chose][i])
  return Dija_merged

def merge_Dij_by_known_id(Dija,Dijb,ia,ib):
  tol=0.03
  ia_ib_matchscores=[]
  Dija_merged={}
  for iai in range(len(Dija)):
    Dija_merged[iai]=[]
    for i in range(len(Dija[iai])):
      Dija_merged[iai].append(Dija[iai][i])
  same_index,diff_index,Dijb_can_be_merged=same_and_diff(ia,ib,Dija,Dijb,tol,2)
  if len(same_index)==0:
    print("no same near atoms. cannot generate Dij. try to enlarge Rnearest to 3")
    same_index,diff_index,Dijb_can_be_merged=same_and_diff(ia,ib,Dija,Dijb,tol,3)
  if len(same_index)==0:
    fraction=0 
    print("these two strucures dont have same surrding, keep the original strucurres and merge them directly ")
    for i in range(len(Dijb[ib])):
      x_vec_y_vec_z_vec=[Dijb[ib][i][2][0],Dijb[ib][i][2][1],Dijb[ib][i][2][2]]
      Dijb[ib][i][2][0]=x_vec_y_vec_z_vec[0]
      Dijb[ib][i][2][1]=x_vec_y_vec_z_vec[1]
      Dijb[ib][i][2][2]=x_vec_y_vec_z_vec[2]
      Dija_merged[ia].append(Dijb[ib][i])
  else:
    fraction=len(np.unique(same_index[:,0]))/len(Dija[ia])
    ia_ib_matchscore=[ia,ib,len(same_index),fraction,len(diff_index)]
    print("the same number is "+str(int(ia_ib_matchscore[2]))+". merge Dija " + str(int(ia_ib_matchscore[0]))+ " and Dijb " + str(int(ia_ib_matchscore[1]))+", with the fraction of "+ str(ia_ib_matchscore[3]))
    x_vec_y_vec_z_vec_a=Dija[ia][int(same_index[0][0])][2] 
    x_vec_y_vec_z_vec_b=Dijb[ib][int(same_index[0][1])][2]
    M_trans=cal_the_tranM(x_vec_y_vec_z_vec_b,x_vec_y_vec_z_vec_a)
    print("the M trans (tranverse metrix to merge two different local coordination) is")
    print(M_trans)
    for i in range(len(Dijb[ib])):
      x_vec_y_vec_z_vec=[Dijb[ib][i][2][0],Dijb[ib][i][2][1],Dijb[ib][i][2][2]]
      x_vec_y_vec_z_vec=np.array(x_vec_y_vec_z_vec)
      x_vec_y_vec_z_vec=x_vec_y_vec_z_vec*M_trans
      x_vec_y_vec_z_vec=np.array(x_vec_y_vec_z_vec)
      #print(i)
      Dijb[ib][i][2][0]=x_vec_y_vec_z_vec[0]
      Dijb[ib][i][2][1]=x_vec_y_vec_z_vec[1]
      Dijb[ib][i][2][2]=x_vec_y_vec_z_vec[2]
      Dija_merged[ia].append(Dijb[ib][i])
  return Dija_merged

def cal_the_tranM(x_vec_y_vec_z_vec_a,x_vec_y_vec_z_vec_b):
  #x_vec_y_vec_z_vec_a*M_trans=x_vec_y_vec_z_vec_b
  x1=x_vec_y_vec_z_vec_a[0]
  y1=x_vec_y_vec_z_vec_a[1]
  z1=x_vec_y_vec_z_vec_a[2]
  A=[x1,y1,z1]
  A=np.array(A)
  #print(A)
  b1=[x_vec_y_vec_z_vec_b[0][0],x_vec_y_vec_z_vec_b[1][0],x_vec_y_vec_z_vec_b[2][0]]
  b2=[x_vec_y_vec_z_vec_b[0][1],x_vec_y_vec_z_vec_b[1][1],x_vec_y_vec_z_vec_b[2][1]]
  b3=[x_vec_y_vec_z_vec_b[0][2],x_vec_y_vec_z_vec_b[1][2],x_vec_y_vec_z_vec_b[2][2]]
  M_trans1=np.linalg.solve(A, b1)
  M_trans2=np.linalg.solve(A, b2)
  M_trans3=np.linalg.solve(A, b3)
  M_trans=[M_trans1,M_trans2,M_trans3]
  M_trans=np.mat(M_trans)
  M_trans=np.transpose(M_trans)
  #M_trans=np.array(M_trans)
  return M_trans


def gen_extrapos_from_Dij_assby_atoma(Dija,atomsa):
  extra_positions=[]
  extra_positions_ele=''
  for i in range(len(Dija)):
    pos_ref=atomsa.get_positions()[i]
    for j in range(len(Dija[i])):
      xija=Dija[i][j][0][1]/(Dija[i][j][0][0]**2)
      yija=Dija[i][j][0][2]/(Dija[i][j][0][0]**2)
      zija=Dija[i][j][0][3]/(Dija[i][j][0][0]**2)  
      #typical example: Dija[i][j][0]=[0.19027283845980128, 0.11513382090596258, -0.09355007556563862, -0.11914839361404195, 6]
      #typical example2: Dija[i][j][1]=H
      #typical example2: Dija[i][j][2]=[x_vec,y_vec,z_vec]
      #vect=xija*x_veca[i][j]+yija*y_veca[i][j]+zija*z_veca[i][j]
      vect=xija*Dija[i][j][2][0]+yija*Dija[i][j][2][1]+zija*Dija[i][j][2][2]
      extra_position=pos_ref+vect
      extra_positions.append((extra_position[0],extra_position[1],extra_position[2]))
      extra_positions_ele+=Dija[i][j][1]
  
  extrapos=Atoms(extra_positions_ele,extra_positions,cell=atomsa.get_cell(),pbc=True) 
  return extrapos


def same_and_diff(ia,ib,Dija,Dijb,tol,Rnearest):  #i_from_a,i_from_b,Dij_from_a,Dij_from_b,tolerance,Rnearest means larger than this distance, even if the cordiation is the same , we will also not count.
  same_index=[]
  diff_index=[] 
  ja=0
  #Dija[ia].tolist()
  #Dijb[ib].tolist()
  Dijb_can_be_used=[]
  for Da_ja in Dija[ia]:
    jb=0
    for Db_jb in Dijb[ib]:
      ##typical example: Dija[i][j][0]=[0.19027283845980128, 0.11513382090596258, -0.09355007556563862, -0.11914839361404195, 6]
      #typical example2: Dija[i][j][1]=H
      #typical example2: Dija[i][j][2]=[x_vec,y_vec,z_vec]
      #if jb==1:
        #print("Da_ja[1] ="+Da_ja[1])
        #print("Da_ja[0] ="+str(Da_ja[0]))
      dev=(np.linalg.norm(np.array(Da_ja[0][0:3])-np.array(Db_jb[0][0:3])))/(Da_ja[0][0]+Da_ja[0][0])
      if Da_ja[1] == Db_jb[1] and dev<tol and Da_ja[0][0] > 1/Rnearest:
        
        #print("similiarty"+str(np.linalg.norm(np.array(Da_ja[0][0:3])-np.array(Db_jb[0][0:3]))))
        same_index.append([ja,jb,np.around(dev,decimals=3),Da_ja[1]+str(Da_ja[0][4]+1),Db_jb[1]+str(Db_jb[0][4]+1)])
      else:
        #print(Db_jb)
        Dijb_can_be_used.append(Db_jb)   #bug
        diff_index.append([int(Da_ja[0][4]),int(Db_jb[0][4])])
      jb=jb+1
    ja=ja+1
  same_index=np.array(same_index)
  return (same_index,diff_index,Dijb_can_be_used)

