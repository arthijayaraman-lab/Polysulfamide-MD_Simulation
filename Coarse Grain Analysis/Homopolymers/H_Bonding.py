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
seed_dist = [0,1500,2000,2500]

aliL_dist = [3]
aliR_dist = [3]

type_pos_dis = ["Configuration_1","Configuration_2"]

save_file_directory = '/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Aliphatic/Data/'

if nChain == 100:
    len_box =80
if nChain == 98:
    len_box =100
if nChain == 267:
    len_box =160
    
## Simulated annealing conditions
eps_start = 6       # Starting of simulated annealing
gap =0.2            # Step rate in simulated annealing
eps_end=12          # Highest interaction strength of simulated annealing

eps = np.linspace(eps_start,eps_end,freq_eps+1)



def analyze_HBond(nRU, nChain,left_distribution, right_distribution, model_type, box_len, eps_start, gap, eps_end, save_file_directory):
    """    
    Calculating the propensity of H-bonding

    Parameters:
        
    nRU (int): Number of repeat units.
    nChain (int): Number of chains in the simulation box.
    left_distribution (list or array): Array representing the distribution of beads on the left side of monomer.
    right_distribution (list or array): Array representing the distribution of beads on the right side of monomer.
    model_type (list or array): Array showing the various configurations used in the simulation.
    
    box_len (float): Length of the simulation box.

    # Simulated annealing conditions:

    eps_start (float): Initial interaction strength for simulated annealing.
    gap (float): Step size for interaction strength increments during simulated annealing.
    eps_end (float): Maximum interaction strength for simulated annealing.
    
    Returns:
    file: A file containing interaction strengths and the corresponding number of hydrogen bonds.
    
    """
    for type_pos in type_pos_dis:
        for l in range(len(aliL_dist)):
            aliL = aliL_dist[l]
            aliR = aliR_dist[l]
            for trial in range(1,4):
        	
                seed = seed_dist[trial]
                freq_eps = int(math.ceil((eps_end-eps_start)/gap))
   
                files = load_trajectories(nRU, nChain, aliL, aliR, seed,eps_start,eps_end,freq_eps,trial,type_pos)
                num_files = np.shape(files)[0]
                
                hb_item = np.zeros((num_files,10000,2))
                normalize_hb = []
                cutoff = 0.35
                
                step = 125000
                
                #import time
                #start = time.process_time()
                
                time_com = []
                
                Timestep = np.linspace(40,80,11) 
                
                # GRID-based approach for faster CPU calculations:  
                # 1. Divide simulation box into smaller grid cells.  
                # 2. Add padding layers for cross-cell interactions.  
                # 3. Analyze atoms only within each grid cell.  
                # Reduces interactions and speeds up calculations.  
                
                grid_size = 30
                
                for timestep in Timestep:
                    time_cal = (timestep)*step
                    file_name = save_file_directory + str(type_pos)+'/hbond_disulf_nRU'+ str(int(nRU)) + '_nChain'+str(int(nChain))+'_leftsataliph'+str(int(aliL))+'_rightsataliph'+str(int(aliR))+'_timestep'+str(int(time_cal))+'_trial'+str(int(trial))+'.txt'
                    f = open("" +file_name, "w") 
                 
                    for n_file in range(num_files):
                        lines =[]
                        file= files[n_file]
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
                        
                        donor =[]
                        acceptor =[]
                        
                        donor = atom[atom[:, 1] == donor_type]
                        acceptor = atom[atom[:, 1] == acceptor_type]
                        
                        donor = donor[donor[:, 1].argsort()]
                        acceptor = acceptor[acceptor[:, 1].argsort()]
                        donoritem = donor[:,0]
                        acceptoritem = acceptor[:,0]
                        
                        donor_coordinate = donor[:,2:]
                        acceptor_coordinate= acceptor[:,2:]

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
                                                break
                                    
                        normalize = count_hb/(2*2*nRU*nChain)
                        normalize_hb.append(normalize)
                        print("Hydrogen Bonds in Eps %s is: %s" %(eps[n_file],count_hb))
                
                        T_current_eps = round(eps[n_file],1)
                        f.write(str(float(T_current_eps)) + "," +str(float(normalize)))
                        f.write('\n')
                    f.close()

