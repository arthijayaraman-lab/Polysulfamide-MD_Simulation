# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 08:19:02 2023

@author: jayshah
"""

#After the relaxtion
import numpy as np
import matplotlib.pyplot as plt

def read_file(file):
    with open(file) as data:
        lines = data.readlines()
    
    atoms = int(lines[3].split()[0])
    
    ts=200
    xpos = np.zeros([ts, atoms])
    ypos = np.zeros([ts, atoms])
    zpos = np.zeros([ts, atoms])
    
    timesteps = []
    box_lens = np.zeros([ts,3])
    
    atype= np.zeros([ts,atoms], dtype = 'int')
    molecnum= np.zeros([ts,atoms], dtype = 'int')
    
    count =0
    for l_num, line in enumerate(lines):
        line = line.split()
        if (len(line) ==2) and (line[1] == "TIMESTEP"):  
            timesteps.append(int(lines[l_num+1][:-1]))
            xlen= float(lines[l_num+5].split()[1]) - float(lines[l_num+5].split()[0])
            ylen= float(lines[l_num+6].split()[1]) - float(lines[l_num+6].split()[0])
            zlen= float(lines[l_num+7].split()[1]) - float(lines[l_num+7].split()[0])
            box_lens[count] =[xlen,ylen,zlen]
            
            for i in range(atoms):
                molecnum[count][i]=int(lines[l_num+9+i].split()[0])
                atype[count][i] = int(lines[l_num+9+i].split()[1])
                xpos[count][i] = float(lines[l_num+9+i].split()[2])*xlen+ float(lines[l_num+5].split()[0])
                ypos[count][i] = float(lines[l_num+9+i].split()[3])*ylen+ float(lines[l_num+6].split()[0])
                zpos[count][i] = float(lines[l_num+9+i].split()[4])*zlen+ float(lines[l_num+7].split()[0])           
           
            count += 1     
            
    atype=atype[:count,:]
    molecnum=molecnum[:count,:]
    xpos=xpos[:count,:]
    ypos=ypos[:count,:]
    zpos=zpos[:count,:]
    box_lens=box_lens[:count,:]
    return (atype,molecnum,xpos,ypos,zpos,box_lens,atoms,timesteps)
        
     
def pbc_distance(r1,r2, len_box):
    delx=r1[0]-r2[0]
    dely=r1[1]-r2[1]
    delz=r1[2]-r2[2]
    delx=delx-len_box*round(delx/len_box) #Periodic boundary conditions
    dely=dely-len_box*round(dely/len_box)
    delz=delz-len_box*round(delz/len_box)
    r2=delx**2 + dely**2 + delz**2
    return(r2**0.5)  
def pbc_distance_difference(r1,r2, len_box):
    delx=r1[0]-r2[0]
    dely=r1[1]-r2[1]
    delz=r1[2]-r2[2]
    delx=delx-len_box*round(delx/len_box) #Periodic boundary conditions
    dely=dely-len_box*round(dely/len_box)
    delz=delz-len_box*round(delz/len_box)
    r2=delx**2 + dely**2 + delz**2
    return [delx,dely,delz]  
def angle_vector_scale(r1,r2):
    len_box = 1 
    delx=r1[0]-r2[0]
    dely=r1[1]-r2[1]
    delz=r1[2]-r2[2]
    delx=delx-len_box*round(delx/len_box) #Periodic boundary conditions
    dely=dely-len_box*round(dely/len_box)
    delz=delz-len_box*round(delz/len_box)
    
    return [delx,dely,delz]       