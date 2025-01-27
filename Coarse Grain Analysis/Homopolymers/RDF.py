
"""
@author: jayshah

Radial Distribution function
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from Read_matrix_file import *
from Filename import *
nRU = 5
nChain = 100

aliL_dist = [3]
aliR_dist = [3]
type_pos_dis = ["Configuration_1","Configuration_2"]

eps_range = [6,8,9,10,11,12]


def rdf(nRU, nChain,left_distribution, right_distribution, model_type, box_len, eps_range, save_file_directory):
    """    
    
    Radial distribution function between Sulfamide bead - Sulfmaide bead

    Parameters:
        
    nRU (int): Number of repeat units.
    nChain (int): Number of chains in the simulation box.
    left_distribution (list or array): Array representing the distribution of beads on the left side of monomer.
    right_distribution (list or array): Array representing the distribution of beads on the right side of monomer.
    model_type (list or array): Array showing the various configurations used in the simulation.
    
    
    box_len (float): Length of the simulation box.
    
    # Simulated annealing conditions
    
    eps_range (array): Array of interaction strengths used to evaluate the radial distribution function (RDF).
    
    Returns:
    file: A file containing bin values and corresponding counts for the radial distribution function (RDF).
    
    """
    for type_pos in type_pos_dis:
        for l in range(len(aliL_dist)):
            aliR = aliR_dist[l]
            aliL = aliL_dist[l]
            seed_dist = [0,1500,2000,2500]
            #seed = 1500
            for trial in range(1,4):
                seed = seed_dist[trial] 
                
                
                for current_eps in eps_range:
        
                    file = file_name(nRU, nChain, aliL, aliR, seed, current_eps,trial,type_pos)
                    
                    #input - Nru, Nchain, left, right, seeed, eps
                    y=read_file(file)
                    
                    #rdf of type 33
                    step = 125000
                    Timestep = np.linspace(40,80,11)
                    end_timestep = 10
                    
                    
                    time_cal = []
                    atoms_types =[]
                    moleculnum=[]
                    x_pos=[]
                    y_pos=[]
                    zpos=[]
                    box=[]
                    
                    for i in range(len(Timestep)):
                        a = int(Timestep[i])
                        time_cal.append((Timestep[i])*step)
                        #print(y[0][a])
                        atoms_types.append(y[0][a])
                        moleculnum.append(y[1][a])
                        x_pos.append(y[2][a])
                        y_pos.append(y[3][a])
                        zpos.append(y[4][a])
                        box.append(y[5][a])
                    
                    for t in range(len(Timestep)):
                        file_name1 = save_file_directory +str(type_pos)+'/rdf_disulf_nRU'+ str(int(nRU)) + '_nChain'+str(int(nChain))+'_leftsataliph'+str(int(aliL))+'_rightsataliph'+str(int(aliR))+'_timestep'+str(int(time_cal[t]))+'_eps'+str(float(current_eps))+'_trial'+str(int(trial))+'.txt'
                        f = open("" +file_name1, "w")
                        
                        rmax = 20
                        dr = 0.05
                        nbins = int(rmax/dr)
                        volumes = [0]*nbins
                        counts = [0]*nbins
                        rdf = [0]*(nbins)
                        Natoms = y[6]
                        
                        bin_radius = np.arange(0, rmax+dr, dr)
                        volumes = [4*math.pi*(((r+dr)**3 - (r**3)))/3 for r in bin_radius]
                        
                        type_3_atom = np.zeros((10000,3))
                        num_type3 = 0
                        item_num_3=np.zeros(10000)
                        chain_length = (aliL+aliR+2*5)*nRU
                        
                        for n in range(Natoms):
                            if (atoms_types[t][n] == 3):
                                type_3_atom[num_type3][:] = [x_pos[t][n],y_pos[t][n],zpos[t][n]]
                                item_num_3[num_type3] = moleculnum[t][n]
                                num_type3 += 1
                                
                        type_3_atom = type_3_atom[:num_type3]
                        len_box = box[t][0]
                        
                        for i in range(num_type3):
                            for j in range(i+1,num_type3):
                                if(item_num_3[i]//chain_length != item_num_3[j]//chain_length ):
                                    r1 = type_3_atom[i]
                                    r2 = type_3_atom[j]
                                    r = pbc_distance(r1, r2, len_box)
                                    if r < rmax:
                                        bin_i = int(r / dr )# round down
                                        counts[bin_i] += 2
                                
                         #Normalize the bin counts and store rdf
                        for h in range(nbins):       
                            rdfshell = (counts[h])/volumes[h]
                            particles= (2*nRU*nChain)
                            
                            
                           
                            avg_dens = particles
                            vol_box = len_box**3
                            avg_dens = particles / vol_box
                            rdf[h] = rdfshell/particles/avg_dens
                        
                        r_vals = (bin_radius[1:]+bin_radius[:-1])/2
                        
                        for p in range(np.shape(r_vals)[0]):
                            f.write(str(float(r_vals[p])) + "," +str(float(rdf[p])))
                            f.write('\n')
                            
                        plt.plot(r_vals, rdf,label="Timestep %s" %t)
                        plt.legend()
                        plt.title("L_%s R_%s" %(aliL,aliR))
                    
                    
                    f.close()  