# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 12:39:03 2023

@author: jayshah
"""


import numpy as np
import math
import matplotlib.pyplot as plt
from Read_matrix_file import *
from Filename import *
nRU = 5
nChain1_dist = [25,50,75]
type_pos_dis = ["Configuration_1"]

## Simulated annealing conditions
eps_start = 6       # Starting of simulated annealing
gap =0.2            # Step rate in simulated annealing
eps_end=12          # Highest interaction strength of simulated annealing
eps = np.linspace(eps_start,eps_end,freq_eps+1)

current_eps_range = [6,7,12]

save_file_directory = '/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Blends/BiModal/Data/'
def cross_rdf(nRU, nChain1_dist, model_type, box_len, eps_range, save_file_directory):
    """    
    
    Radial distribution function between Sulfamide bead - Sulfmaide bead

    Parameters:
        
    nRU (int): Number of repeat units.
    nChain1_dist(int) - number of chains of type A in simulation box
    model_type (list or array): Array showing the various configurations used in the simulation.
    
    box_len (float) - length of the simulation box
    
    # Simulated annealing conditions
    
    eps_range (array): Array of interaction strengths used to evaluate the radial distribution function (RDF).
    
    Returns:
    file: 3 file containing bin values and corresponding counts for the radial distribution function (RDF) for A-A, B-B and A-B chain types.
    
    """
    for type_pos in type_pos_dis:
        for nChain1 in nChain1_dist:
            nChain1 = nChain1
            nChain2 = int(100-nChain1)
            
            # Define a list of possible values for left side lengths of chain type 1 and iterate over them
            aliL1_dist = [2, 4]
            
            for aliL1 in aliL1_dist:
                aliL2 = 6  # Set fixed value for left side lengths of chain type 2
                aliR1 = 12 - aliL1  # Calculate right side  based on left side for chain type 1
                aliR2 = 12 - aliL2  # Calculate right side  based on left side for chain type 2
    
                for current_eps in current_eps_range:
                   
                    trial_range = [1,2,3]
                    seed_dist = [0,1500,2000,2500]
                    for trial in trial_range:
                        seed = seed_dist[trial]
                        step = 125000
                        Timestep = np.linspace(40,80,11) 
                        last=int(step*Timestep[-1])
        
                        file = file_name(nRU, nChain1,nChain2, aliL1,aliL2, aliR1,aliR2, seed, current_eps,trial,last,type_pos)
                        
                        y=read_file(file)
                        
    
                        
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
                            file2 = save_file_directory+str(int(last))+'/'+str(type_pos)+'/rdf11_disulf_nRU'+ str(int(nRU)) + '_nChain_'+str(int(nChain1))+'_'+str(int(nChain2))+'_leftsataliph_'+str(int(aliL1))+'_'+str(int(aliL2))+'_rightsataliph_'+str(int(aliR1))+'_'+str(int(aliR2))+'_eps'+str(float(current_eps))+'_timestep'+str(int(time_cal[t]))+'_trial'+str(int(trial))+'.txt'
                            f1 = open("" +file2, "w")
                    
                            file3 = save_file_directory+str(int(last))+'/'+str(type_pos)+'/rdf22_disulf_nRU'+ str(int(nRU)) + '_nChain_'+str(int(nChain1))+'_'+str(int(nChain2))+'_leftsataliph_'+str(int(aliL1))+'_'+str(int(aliL2))+'_rightsataliph_'+str(int(aliR1))+'_'+str(int(aliR2))+'_eps'+str(float(current_eps))+'_timestep'+str(int(time_cal[t]))+'_trial'+str(int(trial))+'.txt'
                            g1 = open("" +file3, "w")
                    
                            file4 = save_file_directory+str(int(last))+'/'+str(type_pos)+'/rdfdiff_disulf_nRU'+ str(int(nRU)) + '_nChain_'+str(int(nChain1))+'_'+str(int(nChain2))+'_leftsataliph_'+str(int(aliL1))+'_'+str(int(aliL2))+'_rightsataliph_'+str(int(aliR1))+'_'+str(int(aliR2))+'_eps'+str(float(current_eps))+'_timestep'+str(int(time_cal[t]))+'_trial'+str(int(trial))+'.txt'
                            h1 = open("" +file4, "w")
                            
                            rmax = 40
                            dr = 0.05
                            nbins = int(rmax/dr)
                            volumes = [0]*nbins
                            counts_1_1 = [0]*nbins
                            counts_2_2 = [0]*nbins
                            counts_diff = [0]*nbins
                            rdf_1_1 = [0]*(nbins)
                            rdf_2_2 = [0]*(nbins)
                            rdf_diff = [0]*(nbins)
                            Natoms = y[6]
                            
                            bin_radius = np.arange(0, rmax+dr, dr)
                            volumes = [4*math.pi*(((r+dr)**3 - (r**3)))/3 for r in bin_radius]
                            
                            type_3_atom = np.zeros((5000,3))
                            num_type3 = 0
                            item_num_3=np.zeros(5000)
                            chain_length = (aliL1+aliR1+2*5)*nRU
                            
                            nChain = nChain1 +nChain2
                            for n in range(Natoms):
                                if (atoms_types[t][n] == 3):
                                    type_3_atom[num_type3][:] = [x_pos[t][n],y_pos[t][n],zpos[t][n]]
                                    item_num_3[num_type3] = moleculnum[t][n]
                                    num_type3 += 1
                                    
                            type_3_atom = type_3_atom[:num_type3]
                            len_box = box[t][0]
                            
                            print(len_box)
                            for i in range(num_type3):
                                for j in range(i+1,num_type3):
                                    if(item_num_3[i]//chain_length != item_num_3[j]//chain_length ):
                                        if (item_num_3[i]//chain_length <nChain1 ):          # FOr crossed RDF
                                            if (item_num_3[j]//chain_length <nChain1 ):
                                                r1 = type_3_atom[i]
                                                r2 = type_3_atom[j]
                                                r = pbc_distance(r1, r2, len_box)
                                                if r < rmax:
                                                    bin_i = int(r / dr )# round down
                                                    counts_1_1[bin_i] += 2
                                            else:
                                                r1 = type_3_atom[i]
                                                r2 = type_3_atom[j]
                                                r = pbc_distance(r1, r2, len_box)
                                                if r < rmax:
                                                    bin_i = int(r / dr )# round down
                                                    counts_diff[bin_i] += 1
                                        else:          
                                            if (item_num_3[j]//chain_length < nChain1 ):
                                                r1 = type_3_atom[i]
                                                r2 = type_3_atom[j]
                                                r = pbc_distance(r1, r2, len_box)
                                                if r < rmax:
                                                    bin_i = int(r / dr )# round down
                                                    counts_diff[bin_i] += 1
                                            else:
                                                r1 = type_3_atom[i]
                                                r2 = type_3_atom[j]
                                                r = pbc_distance(r1, r2, len_box)
                                                if r < rmax:
                                                    bin_i = int(r / dr )# round down
                                                    counts_2_2[bin_i] += 2
                            #Normalize the bin counts and store rdf
                            for h in range(nbins):       
                                rdfshell = (counts_1_1[h])/volumes[h]
                                particles= (2*nRU*nChain)*(nChain1/nChain)
                                vol_box = len_box**3
                                avg_dens = particles / vol_box
                                rdf_1_1[h] = rdfshell/particles/avg_dens
                                
                                rdfshell = (counts_2_2[h])/volumes[h]
                                particles= (2*nRU*nChain)*(nChain2/nChain)
                                vol_box = len_box**3
                                avg_dens = particles / vol_box
                                rdf_2_2[h] = rdfshell/particles/avg_dens
                                
                                rdfshell = (counts_diff[h])/volumes[h]
                                particles1 = (2*nRU*nChain)*(nChain1/nChain)
                                particles2 = (2*nRU*nChain)*(nChain2/nChain)
                                vol_box = len_box**3
                                avg_dens = particles1 / vol_box
                                rdf_diff[h] = rdfshell/particles2/avg_dens
                            
                            r_vals = (bin_radius[1:]+bin_radius[:-1])/2
                            
                            for p in range(np.shape(r_vals)[0]):
                                f1.write(str(float(r_vals[p])) + "," +str(float(rdf_1_1[p])))
                                f1.write('\n')
                                g1.write(str(float(r_vals[p])) + "," +str(float(rdf_2_2[p])))
                                g1.write('\n')
                                h1.write(str(float(r_vals[p])) + "," +str(float(rdf_diff[p])))
                                h1.write('\n')
                       
                        
                        f1.close()  
                        g1.close()
                        h1.close()