"""
@author: jayshah
"""

import numpy as np
import math
from Filename import *
import itertools
from Read_matrix_file import * 

nRU = 5
nChain = 100

aliL_dist = [3]
aliR_dist = [3]
type_pos_dis = ["Configuration_1","Configuration_2"]

save_file_directory = '/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Aliphatic/Data/'

## Simulated annealing conditions
eps_start = 6       # Starting of simulated annealing
gap =0.2            # Step rate in simulated annealing
eps_end=12          # Highest interaction strength of simulated annealing


eps = np.linspace(eps_start,eps_end,freq_eps+1)

def angles_between_HBond_chains(nRU, nChain,left_distribution, right_distribution, model_type, eps_start, gap, eps_end, save_file_directory):
    
    """    
    Calculating the angle between chains with at least one H-bonding interaction.

    
    Parameters:

    nRU (int): Number of repeat units.
    nChain (int): Number of chains in the simulation box.
    left_distribution (list or array): Array representing the distribution of beads on the left side of monomer.
    right_distribution (list or array): Array representing the distribution of beads on the right side of monomer.
    model_type (list or array): Array showing the various configurations used in the simulation.
    
    box_len (float): Length of the simulation box.
    
    # Simulated annealing conditions:
    
    eps_start (float): Initial interaction strength for simulated annealing.
    gap (float): Step size or increment for interaction strength during simulated annealing.
    eps_end (float): Maximum or highest interaction strength during simulated annealing.
    
    Returns:
    file: A file containing the distribution of angles, with the following data:
        - Bin edges representing the angle intervals.
        - Counts of angles falling into each bin.
    """
    
    for type_pos in type_pos_dis:
        for l in range(len(aliL_dist)):
            aliL = aliL_dist[l]
            aliR = aliR_dist[l]
            seed_dist = [0,1500,2000,2500]      # Various seeds that I used (change accordingly)

            for trial in range(1,4):
                seed = seed_dist[trial]              
                
                
                if nRU == 5:
                    len_box =80     # Specify Box length of the simulaiton in the production run
           
                freq_eps = int(math.ceil((eps_end-eps_start)/gap))    
                files = load_trajectories(nRU, nChain, aliL, aliR, seed,eps_start,eps_end,freq_eps,trial,type_pos)
                num_files = np.shape(files)[0]
                
                hb_item = np.zeros((num_files,10000,2))
                normalize_hb = []
                cutoff = 0.35
                
                step = 125000 # Step size of collecting simulation snapshot
                
                
                #import time
                #start = time.process_time()
                time_com = []
                
                Timestep = np.linspace(40,80,11) # Collecting last 11 snapshots
                
                # GRID-based approach for faster CPU calculations:  
                # 1. Divide simulation box into smaller grid cells.  
                # 2. Add padding layers for cross-cell interactions.  
                # 3. Analyze atoms only within each grid cell.  
                # Reduces interactions and speeds up calculations.  
                
                
                grid_size = 10
                #start = time.process_time()
                for n_file in range(num_files):
                    T_current_eps = round(eps[n_file],1)
                    
                    file_name = save_file_directory+str(type_pos)+'/Angle_disulf_nRU'+ str(int(nRU)) + '_nChain'+str(int(nChain))+'_leftsataliph'+str(int(aliL))+'_rightsataliph'+str(int(aliR))+'_eps'+str(float(T_current_eps))+'_trial'+str(int(trial))+'.txt'
                    f = open("" +file_name, "w") 
                     
                    for timestep in Timestep:
                        time_cal = (timestep)*step
                    
                 
                        f.write(str(int(time_cal)))
                        f.write('\n')
                        lines =[]
                        file= files[n_file]
                        print(file)
                        print("yes")
                        with open(file) as data:
                            lines = data.readlines()
                           
                        
                        donor_type = 1
                        acceptor_type = 2
                       
                        
                        donor_coordinate = np.zeros((2*2*nRU*nChain,3))
                        acceptor_coordinate =np.zeros((2*2*nRU*nChain,3))
                        
                        atoms = int(lines[3])
                
                        donoritem = np.zeros(2*2*nRU*nChain)
                        acceptoritem = np.zeros(2*2*nRU*nChain)
                
                        initial = (9 + atoms) * int(timestep)
                        h = lines[initial-atoms:initial]
                        results=[]
                        atom=[]
                        
                        chain_length = (aliL+aliR+2*5)*nRU
                       
                        for i in range(atoms):
                            results.append([x.strip() for x in h[i].split(' ')])
                            p=[x.strip() for x in h[i].split(' ')]
                            atom.append([int(p[0]),int(p[1]),float(p[2]),float(p[3]),float(p[4])])
                        
                        atom = np.array(atom)
                        atom = atom[atom[:, 0].argsort()]
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
                        hb =[]
                        hb_item =[]
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
                                                hb_item.append([item_donor-1,item_acceptor-1])
                                                hb.append([item_donor,item_acceptor])
                                                break
                                    
                        normalize = count_hb/(2*2*nRU*nChain)
                        normalize_hb.append(normalize)
                        print("Hydrogen Bonds in Eps %s is: %s" %(eps[n_file],count_hb))
                        
                        angle=[]
                        sulfamide_type = 3
                        sulfamide_h_coordinate = []
                        donor_sulfamide=[]
                        donor_side=[]
                        acceptor_sulfamide=[]
                        acceptor_side=[]
                        side_1=[]
                        side_2=[]
                        vector=[]
                        
                        V1_Coordinate = []
                        V2_Coordinate = []
                        
                        sulfamide_list =[]
                        
                       
                       
                        for i in range(np.shape(hb_item)[0]):
                            if int(atom[int(hb_item[i][0])-1][1]) == sulfamide_type:
                                a = atom[int(hb_item[i][0])-1][-4:]
                                sulfamide = np.array([a[1],a[2],a[3]], dtype=float)
                                
                                d_sulfamide_item= int(hb_item[i][0])-1
                                
                                b = atom[int(hb_item[i][0])+4][-4:]
                                s1 = np.array([b[1],b[2],b[3]], dtype=float)
                                
                                vector_1=angle_vector_scale(sulfamide,s1)
                

                            elif int(atom[int(hb_item[i][0])-3][1]) == sulfamide_type:
                                a =atom[int(hb_item[i][0])-3][-4:]
                                sulfamide = np.array([a[1],a[2],a[3]], dtype=float)
                                
                                d_sulfamide_item= int(hb_item[i][0])-3
                                hb = atom[int(hb_item[i][0])+2][-4:]
                                s1 = np.array([hb[1],hb[2],hb[3]], dtype=float)

                               
                                vector_1=angle_vector_scale(sulfamide,s1)
                                
                                
                            
                            acceptor_sulfamide =[]
                            acceptor_side_chain=[]
                            
                           
                            if int(atom[int(hb_item[i][1])-2][1]) == sulfamide_type:
                                l = atom[int(hb_item[i][1])-2][-4:]
                                psulfamide = np.array([l[1],l[2],l[3]], dtype=float)
                                acceptor_sulfamide.append(l[0])
                                
                                
                                a_sulfamide_item= int(hb_item[i][1])-2
                                
                                pb = atom[int(hb_item[i][1])+3][-4:]
                                ps1 = np.array([pb[1],pb[2],pb[3]], dtype=float)

                                vector_2=angle_vector_scale(psulfamide,ps1)
                                
                               
                            
                            elif int(atom[int(hb_item[i][1])-4][1]) == sulfamide_type:
                                l = atom[int(hb_item[i][1])-4][-4:]
                                psulfamide = np.array([l[1],l[2],l[3]], dtype=float)

                                
                                a_sulfamide_item= int(hb_item[i][1])-4
                                
                                pb = atom[int(hb_item[i][1])+1][-4:]
                                ps1 = np.array([pb[1],pb[2],pb[3]], dtype=float)
                                acceptor_side_chain.append(pb[0])
                                
                                vector_2=angle_vector_scale(psulfamide,ps1)
                            
                            
                            if set([d_sulfamide_item,a_sulfamide_item]) not in sulfamide_list[:]:
                                sulfamide_list.append(set([d_sulfamide_item,a_sulfamide_item]))

                                vector_1 = np.transpose(vector_1)
                                vector_2 = np.transpose(vector_2)
                                unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
                                unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
                                dot_product = np.dot(unit_vector_1, unit_vector_2)
                                angles = np.arccos(dot_product)
                                angle_vector = math.degrees(angles)
                                
                                # Append the angle to the list, ensuring it stays within 0-90 degrees.
                                if angle_vector < 90:
                                    angle.append(angle_vector)
                                else:
                                    angle.append(180-angle_vector)
                            
     
                        div=5                
                        num = int(90/div)  +1
                        bins_count = np.zeros((num_files,num))
                        bins = np.linspace(0,90,num)
                        for i in range(np.shape(angle)[0]):
                            bins_count[n_file][int(angle[i]/div)] += 1
                        
                        
                        for i in range(num):
                            
                            f.write(str(int(bins[i])) + "," +str(int(bins_count[n_file][i])))
                            f.write('\n')
                    
                    f.close()
            
            
                #end = time.process_time()
                #print("Time taken = " + str(1000*(end - start)) + " ms")   
          