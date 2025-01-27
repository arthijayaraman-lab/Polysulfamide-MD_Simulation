# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:10:33 2023

@author: jayshah
"""
def get_keys_from_value(d, val):
    b= [ k for k, v in d.items() if v == val]
    return b[0]


import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from Filename import *
import itertools
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


nRU_dist = [5]
type_pos_dis = ["Configuration_1"]
nChain1_dist = [25,50,75]
aliL1_dist = [2,4,2]
aliL2_dist = [6,6,4]
seed_dist = [0,1500,2000,2500]

## Simulated annealing conditions
eps_start = 6       # Starting of simulated annealing
gap =0.2            # Step rate in simulated annealing
eps_end=12          # Highest interaction strength of simulated annealing

eps = np.linspace(eps_start,eps_end,freq_eps+1)

save_file_directory = '/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Blends/BiModal/Data/'

def Hb_pairs(nRU, nChain1_dist,aliL1_dist,aliL2_dist, model_type, box_len, eps_start, gap, eps_end, save_file_directory):
    
    """
    Calculate hydrogen bond (HB) pairs based on the given parameters.

    Parameters:
    nRU (int): Number of repeat units in the polymer chain.
    nChain1_dist (list or array): Distribution of the number of chains of type 1.
    aliL1_dist (list or array): Distribution of aliphatic units on the left side of chain type 1.
    aliL2_dist (list or array): Distribution of aliphatic units on the left side of chain type 2.
    model_type (str): Type of model being used for the simulation (e.g., coarse-grained or atomistic).
    box_len (float): Length of the simulation box.
    eps_start (float): Starting epsilon value for interaction strength.
    gap (float): Incremental step size for epsilon.
    eps_end (float): Ending epsilon value for interaction strength.
    save_file_directory (str): Directory path to save the output file.

    Returns:
    file: A file containing the following information:
    - Interaction strength.
    - Count of hydrogen bonds (HB) between type A and type A (HB A-A).
    - Count of hydrogen bonds (HB) between type B and type B (HB B-B).
    - Count of hydrogen bonds (HB) between type A and type B (HB A-B).
    """
    
    for type_pos in type_pos_dis:
        for nRU in nRU_dist:
            
            for nChain1 in nChain1_dist:
                nChain2 = int(100-nChain1)
                
                for dist in range(len(aliL1_dist)):
                    aliL1 = aliL1_dist[dist]
                    aliL2 = aliL2_dist[dist]
                    aliR1 = 12-aliL1
                    aliR2 = 12-aliL2
                    
                    
                    for trial in range(1,4):
                        seed = seed_dist[trial]
                    
                    
                    
                        
                        if nRU == 3:
                            len_box =68
                        
                        if nRU == 5:
                            len_box =80
                        if nRU == 10:
                            len_box =100
                        if nRU == 15:
                            len_box =160
                        
                        freq_eps = int(math.ceil((eps_end-eps_start)/gap))
                        
                        step = 125000
                        Timestep = np.linspace(40,80,11)
                        last=int(step*Timestep[-1])                
                        files = load_trajectories(nRU, nChain1,nChain2, aliL1,aliL2, aliR1,aliR2, seed,eps_start,eps_end,freq_eps,trial,last,type_pos)
                        num_files = np.shape(files)[0]
                        print(files)
                        hb_item = np.zeros((num_files,10000,2))
                        normalize_hb = []
                        cutoff = 0.35

                        #import time

                        time_com = []
                        
                        # GRID-based approach for faster CPU calculations:  
                        # 1. Divide simulation box into smaller grid cells.  
                        # 2. Add padding layers for cross-cell interactions.  
                        # 3. Analyze atoms only within each grid cell.  
                        # Reduces interactions and speeds up calculations.  
                        
                        grid_size = 10
                        #start = time.process_time()
                        for timestep in Timestep:
                            time_cal = (timestep)*step
                            file_name = save_file_directory + str(int(last))+'/'+str(type_pos)+'/hbond_count_diff_disulf_nRU'+ str(int(nRU)) + '_nChain_'+str(int(nChain1))+'_'+str(int(nChain2))+'_leftsataliph_'+str(int(aliL1))+'_'+str(int(aliL2))+'_rightsataliph_'+str(int(aliR1))+'_'+str(int(aliR2))+'_timestep'+str(int(time_cal))+'_trial'+str(int(trial))+'.txt'
                     
                            f = open("" +file_name, "w") 
                         
                            for n_file in range(num_files):
                                lines =[]
                                file= files[n_file]
                                with open(file) as data:
                                    lines = data.readlines()
                                   
                                
                                donor_type = 1
                                acceptor_type = 2
                               
                                nChain = nChain1 +nChain2
                                donor_coordinate = np.zeros((2*2*nRU*nChain,3))
                                acceptor_coordinate =np.zeros((2*2*nRU*nChain,3))
                                
                                atoms = int(lines[3])
                        
                                donoritem = np.zeros(2*2*nRU*nChain)
                                acceptoritem = np.zeros(2*2*nRU*nChain)
                        
                                initial = (9 + atoms) * int(timestep)
                                h = lines[initial-atoms:initial]
                                results=[]
                                atom=[]
                                
                                chain_length = (aliL1+aliR1+2*5)*nRU
                               
                                for i in range(atoms):
                                    results.append([x.strip() for x in h[i].split(' ')])
                                    p=[x.strip() for x in h[i].split(' ')]
                                    atom.append([int(p[0]),int(p[1]),float(p[2]),float(p[3]),float(p[4])])
                                
                                atom = np.array(atom)
                                
                                donor =[]
                                acceptor =[]
                                
                                donor = atom[atom[:, 1] == donor_type]
                                acceptor = atom[atom[:, 1] == acceptor_type]
                                
                                donor = donor[donor[:, 1].argsort()]
                                acceptor = acceptor[acceptor[:, 1].argsort()]
                                donoritem = donor[:,0]
                                acceptoritem = acceptor[:,0]
                                
                                donor_coordinate = donor[:,2:]
                                #donor_coordinate = np.float64(donor_coordinate)
                                acceptor_coordinate= acceptor[:,2:]
                                #acceptor_coordinate = np.float64(acceptor_coordinate)
                                donor_count=len(donor)
                                acceptor_count= len(acceptor)   
                                
                               
                                init = np.array([-15,30,0])
                                final = init + len_box
                                
                                donor_coordinate = donor_coordinate[:]*len_box+init[:]
                                acceptor_coordinate= acceptor_coordinate[:]*len_box+init[:] 
                                
                                #df_donor  = pd.DataFrame((donoritem,donor_coordinate), columns=['a', 'b'])
                                
                                #df_acceptor = pd.DataFrame((acceptoritem,acceptor_coordinate), columns=['a', 'b'])
                                donor ={donoritem[i]:donor_coordinate[i] for i in range(donor_count)}
                                acceptor ={acceptoritem[i]:acceptor_coordinate[i] for i in range(acceptor_count)}
                                
                                
                                padding_dis = 1
                            
                            #Grid
                           
                                xgrids = np.array([init[0]+ i *((final[0]-init[0])/grid_size) for i in range(grid_size+1)])
                                ygrids = np.array([init[1]+ i *((final[1]-init[1])/grid_size) for i in range(grid_size+1)])
                                zgrids = np.array([init[2]+ i *((final[2]-init[2])/grid_size) for i in range(grid_size+1)])
                                grids  = [[xgrids[i],ygrids[j],zgrids[k]] for k in range(len(zgrids)) for j in range(len(ygrids)) for i in range(len(xgrids))]
                                
                                grid_number = {}
                                
                                
                                grid_ = (final[0]-init[0])/grid_size
                                #print(grid_size - 1)
                                
                                count = 0
                                for i in range(grid_size):
                                    for j in range(grid_size):
                                        for k in range(grid_size):
                                            grid_number[count] = [i,j,k]
                                            count += 1
                         
                                keys = np.linspace(0,(grid_size)**3-1,(grid_size)**3,dtype =int)
                                donor_grid ={key: [] for key in keys}
                                
                                for cor in range(3):
                                    grids_cord = [init[cor]+ i *((final[cor]-init[cor])/grid_size) for i in range(grid_size+1)]
                                    donor_coordinate[:,cor]=np.where(donor_coordinate[:,cor] < max(grids_cord), donor_coordinate[:,cor], max(grids_cord))
                                    acceptor_coordinate[:,cor] = np.where(acceptor_coordinate[:,cor] < max(grids_cord), acceptor_coordinate[:,cor], max(grids_cord))
                                    donor_coordinate[:,cor]=np.where(donor_coordinate[:,cor] > min(grids_cord), donor_coordinate[:,cor], min(grids_cord))
                                    acceptor_coordinate[:,cor] = np.where(acceptor_coordinate[:,cor] > min(grids_cord), acceptor_coordinate[:,cor], min(grids_cord))
                                
                                count_multi = 0
                                for i in range(donor_count):
                                    position = []
                                    for cor in range(3):
                                        grids_cord = [init[cor]+ i *((final[cor]-init[cor])/grid_size) for i in range(grid_size+1)]
                                        
                                        
                                        for j in range(grid_size):
                                            if donor_coordinate[i][cor] <= grids_cord[j+1] and donor_coordinate[i][cor] >= grids_cord[j]:  
                                                if donor_coordinate[i][cor] > grids_cord[j] + padding_dis and  donor_coordinate[i][cor] <  grids_cord[j+1] - padding_dis:
                                                    position.append([j])
                                                    break
                                                if donor_coordinate[i][cor] <= grids_cord[j] + padding_dis:
                                                    if j == 0:
                                                        position.append([0,grid_size-1])
                                                    else:
                                                        position.append([j-1,j])
                                                    break
                                                if donor_coordinate[i][cor] >=  grids_cord[j+1] - padding_dis:
                                                
                                                    if j<grid_size:
                                                        position.append([j,0])
                                                    else:
                                                        position.append([j,j+1])
                                                    #print(donor_coordinate[i][cor],position)
                                                    break
                                   
                                    pos = list(itertools.product(position[0], position[1],position[2] ))
                                    for n_pos in range(len(pos)):
                                        position_it = list(pos[n_pos])
                                        a = get_keys_from_value(grid_number, position_it)
                                        
                                        donor_grid[a].append(donoritem[i])
                                
                                
                                acceptor_grid ={key: [] for key in keys}
                                
                                for i in range(acceptor_count):
                                    position = []
                                    for cor in range(3):
                                        
                                        grids_cord = [init[cor]+ i *((final[cor]-init[cor])/grid_size) for i in range(grid_size+1)]
                                        
                                        for j in range(grid_size):
                                            if acceptor_coordinate[i][cor] <= grids_cord[j+1] and acceptor_coordinate[i][cor] >= grids_cord[j]:  
                                                if acceptor_coordinate[i][cor] > grids_cord[j] + padding_dis and  acceptor_coordinate[i][cor] <  grids_cord[j+1] - padding_dis:
                                                    position.append([j])
                                                    break
                                                if acceptor_coordinate[i][cor] <= grids_cord[j] + padding_dis:
                                                    if j == 0:
                                                        position.append([0,grid_size-1])
                                                    else:
                                                        position.append([j-1,j])
                                                    break
                                                if acceptor_coordinate[i][cor] >=  grids_cord[j+1] - padding_dis:
                                                    if j<grid_size:
                                                        position.append([j,0])
                                                    else:
                                                        position.append([j,j+1])
                                                    break
                                    pos = list(itertools.product(position[0], position[1],position[2] ))
                                    for n_pos in range(len(pos)):
                                        position_it = list(pos[n_pos])
                                        a = get_keys_from_value(grid_number, position_it)
                                        
                                        acceptor_grid[a].append(acceptoritem[i])
                                        
                                count_hb = 0
                                count_hb_aa = 0
                                count_hb_bb = 0
                                count_hb_ab = 0
                                hb =[]
                                hb_s =[]
                                hb_d =[]
                                a__ = []
                                for i in range(grid_size**3):
                                    donor_in_grid = donor_grid[i]
                                    acceptor_in_grid =acceptor_grid[i]
                                    
                                    for d in range(np.shape(donor_in_grid)[0]):
                                        for a in range(np.shape(acceptor_in_grid)[0]):
                                            item_donor = int(donor_in_grid[d])
                                            item_acceptor = int(acceptor_in_grid[a])
                                            check = [item_donor, item_acceptor]
                                            # Not double counting
                                            if hb.count(check) != 1 :
                                                if(item_donor//chain_length == item_acceptor//chain_length and abs(item_donor - item_acceptor) <= 5):
                                                    pass
                                                
                                                else:
                                                    x=donor[item_donor][0]-acceptor[item_acceptor][0]
                                                    y=donor[item_donor][1]-acceptor[item_acceptor][1]
                                                    z=donor[item_donor][2]-acceptor[item_acceptor][2]
                                                    
                                                    
                                                    x=x-len_box*round(x/len_box) #Periodic boundary conditions
                                                    y=y-len_box*round(y/len_box)
                                                    z=z-len_box*round(z/len_box)
                                                    r = math.sqrt(x**2 + y**2 + z**2)
                                                    
                                                    if r < cutoff:
                                                        count_hb +=1
                                                        hb.append([item_donor,item_acceptor])
                                                        if item_donor > nChain1* 22*nRU and item_acceptor > nChain1* 22*nRU :
                                                            hb_s.append([item_donor,item_acceptor])
                                                            count_hb_bb += 1
                                                        elif item_donor < nChain1* 22*nRU and item_acceptor < nChain1* 22*nRU :
                                                            hb_s.append([item_donor,item_acceptor])
                                                            count_hb_aa += 1
                                                        else:
                                                            hb_d.append([item_donor,item_acceptor])
                                                            count_hb_ab += 1
                                                        break
                                
                                            
                                normalize_s = (count_hb_bb+count_hb_aa)/(2*2*nRU*nChain)
                                normalization = min(2*2*nRU*nChain1,2*2*nRU*nChain2)
                                normalize_d = count_hb_ab/normalization
                                normalize = count_hb/(2*2*nRU*nChain)
                                normalize_hb.append(normalize)
                                print("Hydrogen Bonds in Eps %s is: %s" %(eps[n_file],count_hb))
                                #print(answers_should[n_file]-count_hb )
                        
                                T_current_eps = round(eps[n_file],1)
                                f.write(str(float(T_current_eps)) + "," +str(float(count_hb_aa))+ "," +str(float(count_hb_bb))+ "," +str(float(count_hb_ab)))
                                f.write('\n')
                            f.close()
