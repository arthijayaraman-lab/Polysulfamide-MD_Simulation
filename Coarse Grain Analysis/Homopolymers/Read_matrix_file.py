# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 08:19:02 2023

@author: jayshah
"""

#After the relaxtion
import numpy as np
import matplotlib.pyplot as plt

def get_key_from_value(d, val):
    """
    Retrieve the first key from a dictionary that matches the given value.

    Parameters:
    d (dict): The dictionary to search.
    val: The value to find the corresponding key for.

    Returns:
    The first key that matches the given value, or None if no match is found.
    """
    keys = [k for k, v in d.items() if v == val]
    return keys[0] if keys else None


def read_file(file):
    
    """
    Read the file
    
    Parameter: 
    file is the LAMMPS trajectory file
    
    Returns: 
    dataframe of molecule type ,molecule number,x-position,y-position,z-position,box-length,number of atoms    
    """
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
            timesteps.append(lines[l_num+1])
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
    return (atype,molecnum,xpos,ypos,zpos,box_lens,atoms)
        
        
def pbc_distance(r1,r2, len_box):
    """
    Calculate the minimum distance between two points in a periodic boundary condition (PBC) box.

    Parameters:
    r1 (array-like): A 3D vector representing the coordinates of the first point.
    r2 (array-like): A 3D vector representing the coordinates of the second point.
    len_box (float or array-like): The length of the simulation box in each dimension. 
                                   If float, the box is assumed to be cubic; 
                                   if array-like, it should contain the lengths in x, y, and z.

    Returns:
    float: The minimum distance between the two points, accounting for PBC.
    """
    
    delx=r1[0]-r2[0]
    dely=r1[1]-r2[1]
    delz=r1[2]-r2[2]
    delx=delx-len_box*round(delx/len_box) #Periodic boundary conditions
    dely=dely-len_box*round(dely/len_box)
    delz=delz-len_box*round(delz/len_box)
    r2=delx**2 + dely**2 + delz**2
    return(r2**0.5)  

def angle_vector_scale(r1,r2):
    
    """
    Calculate the scaled distance between two 3D vectors.
    
    Parameters:
    r1 (array-like): A 3D vector representing the first point or direction.
    r2 (array-like): A 3D vector representing the second point or direction.
    
    Returns:
    float: The scaled distance between the two vectors.
   """
    len_box = 1 
    delx=r1[0]-r2[0]
    dely=r1[1]-r2[1]
    delz=r1[2]-r2[2]
    delx=delx-len_box*round(delx/len_box) #Periodic boundary conditions
    dely=dely-len_box*round(dely/len_box)
    delz=delz-len_box*round(delz/len_box)
    
    return [delx,dely,delz]       