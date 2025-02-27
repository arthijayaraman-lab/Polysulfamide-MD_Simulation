# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 10:19:08 2024

@author: jayshah
"""

import numpy as np
import math
import pandas as pd
from collections import Counter
#%% initializing the positions of the CH2/CH3 and sulfamide group

CH2 = {"C": [ "C", -0.12 ,0.000000  ,   0.000000  ,  0.000000],
       "H1":[ "H", 0.06, 0.000000  ,  0.631044  ,   0.892431],
       "H2":[ "H",0.06, 0.000000   ,  0.631044   ,  -0.892431]}

CH3 = {"C": [ "C", -0.18 ,0.000000  ,   0.000000  ,  0.000000],
       "H1":["H",  0.06 , 0.000000  ,  0.631044  ,   0.892431],
       "H2":["H", 0.06 , 0.000000   ,  0.631044   ,  -0.892431],
       "H3":["H", 0.06  ,-0.892431  ,  -0.631044   ,   0.000000]}

nh = 1.04
Sulfamide = {"N1": [ "N", -0.18 ,0.000000  ,   0.000000  ,  0.000000],
           "NH1":["H",  0.06 , 0.000000  ,  0.000000 ,   -nh],
           "S":["S", 0.06 , 0.000000   ,  1.398   , 0.8075],
           "N2":["N", 0.06 , 0.000000   ,  2.796   ,  0.0000],
           "NH2":["H", 0.06 , 0.000000   ,    2.796   ,  -nh],
           "O1":["O", 0.06 ,-1.227   , 1.398   ,  1.5202],
           "O2":["O", 0.06 , 1.227   , 1.398   ,  1.5202],
           }

# CH2_second = {"CN": ["C", -0.12 ,0.000000  ,   0.000000  ,  0.000000],
#        "H_1N":["H", 0.06, 0.000000  ,  0.631044  ,   0.892431],
#        "H_2N":["H",0.06, 0.000000   ,  0.631044   ,  -0.892431]}

#%% For the moving the atoms positions

L_cc = 1.5350 
L_cn = 1.442 
angle_cc =  109.5
angle_cn = 120
move_cc_y  = 0.5*L_cc * math.cos(math.radians(angle_cc/2))
move_cn_y  = 0.5*L_cn * math.cos(math.radians(angle_cn/2))


# %%

def rotate(theta, pos , direction):
    
    cos_angle = math.cos(math.radians(theta))
    sin_angle = math.sin(math.radians(theta))
    
    if direction == "x":
        rotate_matrix = [[1,0,0],[0,cos_angle,-sin_angle] ,[0,sin_angle, cos_angle]]
    if direction == "y":
        rotate_matrix = [[cos_angle,0,sin_angle],[0,1,0] ,[-sin_angle,0, cos_angle]]
    if direction == "z":
        rotate_matrix = [[cos_angle,-sin_angle,0],[sin_angle,cos_angle,0] ,[0,0, 1]]
    #new_coordinates = np.matmul(rotate_matrix, np.transpose(pos)) 
    new_coordinates =np.dot(rotate_matrix, np.transpose(pos))
    return np.transpose(new_coordinates)

#%%

move_cc_x = L_cc * math.sin(math.radians(angle_cc/2))
move_cn_x = L_cn * math.sin(math.radians(angle_cc/2))
def after_rotation_move(pos, move,direction):
    if direction == "x":
        new = [row + [move,0,0] for row in pos]
        #new = pos+ [move,0,0]
    if direction == "y":
        new = pos+ [0,move,0]
    if direction == "z":
        new = pos+ [0,0,move]    
    return new  

#%% Making -CH2-CH2-

coordinates_CH2 = np.zeros((6,3))
moving_coordinate = [0,0,0]
for i in range(len(CH2)):
    coordinates_CH2[i,:] = CH2[list(CH2.keys())[i]][-3:]
    
coord_atoms =  coordinates_CH2[:3]
print(coord_atoms)
coord_atoms =  coord_atoms + [0,move_cc_y,0] # Moving coordinates first in y direction equals to half of C-C bond length along with taking angle into consideration
coordinates_CH2[:3] = coord_atoms # First CH2 groups positions

# Next set of CH2 group
coord_atoms = [value[-3:] for value in CH2.values()]

rotating = rotate(180, coord_atoms , "x") # rotate along x -axis
moving = after_rotation_move(rotating,move_cc_x,"x" ) # moving along x -axis
coordinates_CH2[3:] = moving

atom_type = ["C1", "H1","H2","C2","H3","H4"]
element_type = ["C", "H","H","C","H","H"]
g = open("Residue_structures/CCA.xyz", "w")
g.write(str(6) +"\n")
g.write("Sulfamide \n")
for i in range(6):
    g.write(str(atom_type[i]) +"    " + str(float(coordinates_CH2[i][0])) +"    " + str(float(coordinates_CH2[i][1]))+"    " + str(float(coordinates_CH2[i][2])) +"\n")
g.close()


CCA_coordinates =pd.DataFrame({'Atom_type': atom_type, 'Element_type': element_type,'x': list(coordinates_CH2[:,0]),'y': list(coordinates_CH2[:,1]),'z': list(coordinates_CH2[:,2])})

#%% CCF Making CH3-CH2-

coordinates_CH3 = np.zeros((7,3))
moving_coordinate = [0,0,0]
for i in range(len(CH3)):
    coordinates_CH3[i,:] = CH3[list(CH3.keys())[i]][-3:]
    
coord_atoms =  coordinates_CH3[:4]
print(coord_atoms)
coord_atoms =  coord_atoms + [0,move_cc_y,0] # Moving coordinates first in y direction equals to half of C-C bond length along with taking angle into consideration
coordinates_CH3[:4] = coord_atoms # First CH2 groups positions

# Next set of CH3 group
coord_atoms = [value[-3:] for value in CH2.values()]
rotating = rotate(180, coord_atoms , "x") # rotate along x -axis
moving = after_rotation_move(rotating,move_cc_x,"x" ) # moving along x -axis
coordinates_CH3[4:] = moving

atom_type = ["C1", "H1","H2","H3","C2","H4","H5"]
element_type = ["C", "H","H","H","C","H","H"]
p = open("Residue_structures/CCF.xyz", "w")
p.write(str(7) +"\n")
p.write("Sulfamide \n")
for i in range(7):
    p.write(str(atom_type[i]) +"    " + str(float(coordinates_CH3[i][0])) +"    " + str(float(coordinates_CH3[i][1]))+"    " + str(float(coordinates_CH3[i][2])) +"\n")
p.close()

CCF_coordinates =pd.DataFrame({'Atom_type': atom_type, 'Element_type': element_type,'x': list(coordinates_CH3[:,0]),'y': list(coordinates_CH3[:,1]),'z': list(coordinates_CH3[:,2])})

#%% CCL Making -CH2-CH3

# Reverse -CH3-CH2
coordinates_CH3[:,0] = -coordinates_CH3[:,0]
coordinates_CH3 = coordinates_CH3+ [move_cc_x,0,0]

atom_type = ["C2", "H3","H4","H5","C1","H1","H2"]
element_type = ["C", "H","H","H","C","H","H"]

p = open("Residue_structures/CCL.xyz", "w")
p.write(str(7) +"\n")
p.write("Sulfamide \n")
for i in range(7):
    p.write(str(atom_type[i]) +"    " + str(float(coordinates_CH3[i][0])) +"    " + str(float(coordinates_CH3[i][1]))+"    " + str(float(coordinates_CH3[i][2])) +"\n")
p.close()

CCL_coordinates_CH2Ch3 =pd.DataFrame({'Atom_type': atom_type,'Element_type': element_type, 'x': list(coordinates_CH3[:,0]),'y': list(coordinates_CH3[:,1]),'z': list(coordinates_CH3[:,2])})

#%% Making - SUS

coordinates_sus = np.zeros((13,3))
coordinates_sus[:3]=rotate(180, coordinates_CH2[:3] , "x")
last_C = coordinates_CH2[3]

# Adding sulfmide group
coord_atoms = [value[-3:] for value in Sulfamide.values()]
rotating = rotate(90, coord_atoms , "z") # rotate along y -axis
rotating = rotate(90, rotating , "x") # rotate along y -axis
# coordinates_sul[6:13] = moving+last_C # First CH2 groups positions
coordinates_sus[3:10] = rotating-last_C # First CH2 groups positions
last_N = coordinates_sus[7]

# Adding CH2 GROUP again

# Shifting Coordinate CH2 bond - Reflecting the CH2 - cH2 along S

Sulfur= coordinates_sus[5,0]

coordinates_CH2_moving = np.zeros((3,3))
coordinates_CH2_moving[:,0] = 2*Sulfur- coordinates_sus[:3,0]
coordinates_CH2_moving[:,1:] = coordinates_sus[:3,1:]
coordinates_sus[10:] = coordinates_CH2_moving 

#Shifting to origin
coordinates_sus = coordinates_sus - coordinates_sus[10]

######## NOTE --------------------------- CHANGING THE CA1 --> CB1 [ but using the oplsa parameter of the alpha carbon only], the naming is because to not change other carbon residues. for bonds


atom_type = ["CB2", "HA3","HA4","N2", "NH2","S","N1","NH1","O1","O2","CB1", "HA1","HA2"]
element_type = ["C","H","H","N", "H","S","N","H","O","O","C", "H","H"]


h = open("Residue_structures/SUS.xyz", "w")
h.write(str(7) +"\n")
h.write("Sulfamide \n")
for i in range(3,10):
    h.write(str(atom_type[i]) +"    " + str(float(coordinates_sus[i][0])) +"    " + str(float(coordinates_sus[i][1]))+"    " + str(float(coordinates_sus[i][2])) +"\n")
h.close()

h = open("Residue_structures/SUS_full.xyz", "w")
h.write(str(13) +"\n")
h.write("Sulfamide \n")
for i in range(13):
    h.write(str(atom_type[i]) +"    " + str(float(coordinates_sus[i][0])) +"    " + str(float(coordinates_sus[i][1]))+"    " + str(float(coordinates_sus[i][2])) +"\n")
h.close()

h = open("Residue_structures/SUS_CH22.xyz", "w")
h.write(str(3) +"\n")
h.write("Sulfamide \n")
for i in range(10,13):
    h.write(str(atom_type[i]) +"    " + str(float(coordinates_sus[i][0])) +"    " + str(float(coordinates_sus[i][1]))+"    " + str(float(coordinates_sus[i][2])) +"\n")
h.close()

SUS_coordinates =pd.DataFrame({'Atom_type': atom_type,'Element_type': element_type, 'x': list(coordinates_sus[:,0]),'y': list(coordinates_sus[:,1]),'z': list(coordinates_sus[:,2])})

#%% Defining monomer - Residue

def monomer (molecules_residue):
    CCF_count = 1
    CCL_count = 1
    SUL_count = 1
    CCA_count = 1
    Residue_number = []
    Residue_atom_wise =[]
    #molecules_in_monomer = ["CH3","CH2","CH2","CH3"]
    if molecules_residue[0] != "CCF":
        print("Change residue!!")
        
    
    coordinates = []
    atom_types = []
    element_types = []
    moving_coordinate =[]
    atom_number = 0 
    for i in range(len(molecules_residue)):
        if i == 0:
            # coordinates.append(np.array(CCL_coordinates_CH3Ch2.iloc[:,-3:]))
            coordinates.append(rotate(180,(np.array(CCF_coordinates.iloc[:,-3:])),"x"))
            atom_types.append(np.array(CCF_coordinates.iloc[:,0]))
            element_types.append(np.array(CCF_coordinates.iloc[:,1]))
            moving_coordinate = np.array(CCF_coordinates.iloc[4,-3:])
            residue_atom= (molecules_residue[i]+',')*7
            
            Residue_atom_wise.append(residue_atom.split(',')[:-1])
            Residue_number.append([CCF_count for j in range(7)])
            CCF_count += 1
            atom_number += 6
        if i == len(molecules_residue) - 1  :
            coordinates.append(np.array(CCL_coordinates_CH2Ch3.iloc[:,-3:])+[moving_coordinate[0],0,0]+[move_cc_x,0,0]-[0,1.2*move_cc_y,0])
            atom_types.append(np.array(CCL_coordinates_CH2Ch3.iloc[:,0]))
            element_types.append(np.array(CCL_coordinates_CH2Ch3.iloc[:,1]))
            residue_atom= (molecules_residue[i]+',')*7
            Residue_atom_wise.append(residue_atom.split(',')[:-1])
            Residue_number.append([CCL_count for j in range(7)])
            atom_number += 6
        else:
            if molecules_residue[i] == "SSF" or molecules_residue[i] == "SSL":
                coordinates.append(np.array(SUS_coordinates.iloc[:,-3:])+[moving_coordinate[0],0,0]+[1.2*move_cc_x,move_cc_y,0])
                #moving_coordinate = np.array(CCL_coordinates_CH3Ch2.iloc[4,-3:])
                moving_coordinate = np.array(coordinates[i][0])                     #Make sure to change this
                residue_atom= (molecules_residue[i]+',')*13
                Residue_atom_wise.append(residue_atom.split(',')[:-1])
                atom_types.append(np.array(SUS_coordinates.iloc[:,0]))
                element_types.append(np.array(SUS_coordinates.iloc[:,1]))
                Residue_number.append([SUL_count for j in range(13)])
                SUL_count += 1
                
            if molecules_residue[i] == "CCA":
                # coordinates.append(np.array(CCA_coordinates.iloc[:,-3:])+[moving_coordinate[0],0,0]+[1.3*move_cc_x,0,0])
                coordinates.append(rotate(180,(np.array(CCA_coordinates.iloc[:,-3:])+[moving_coordinate[0],0,0]+[move_cc_x,0,0]),"x"))
                #moving_coordinate = np.array(CCL_coordinates_CH3Ch2.iloc[4,-3:])
                moving_coordinate = np.array(coordinates[i][3]) 
                residue_atom= (molecules_residue[i]+',')*6
                Residue_atom_wise.append(residue_atom.split(',')[:-1])
                atom_types.append(np.array(CCA_coordinates.iloc[:,0]))
                element_types.append(np.array(CCA_coordinates.iloc[:,1]))
                Residue_number.append([CCA_count for j in range(6)])
                CCA_count += 1
                   
    atom_types=np.concatenate(atom_types)
    element_types=np.concatenate(element_types)
    coordinates = np.concatenate(coordinates)
    Residue_atom_wise = np.concatenate(Residue_atom_wise) 
    
    # # Changing residue type --> To convert the
    #     # CCA - SSF --> CSA - SSF 
        
    number_SSF = int(np.sum((Residue_atom_wise == "SSF")) // 13)
    for i in range (number_SSF):
        ssf_index = np.where(Residue_atom_wise == "SSF")[0]
        Residue_atom_wise[:ssf_index[i*13]][-6:] = "CSA"
        
    number_SSL = int(np.sum((Residue_atom_wise == "SSL")) // 13)
    for i in range (number_SSL):
        ssl_index = np.where(Residue_atom_wise == "SSL")[0]
        Residue_atom_wise[ssl_index[(i+1)*13-1]:][1:7] = "SCA" 
        
    # ## Adding a new residue if the SUL are sepearted by just one -C2H4 group --> SCS new residue    
    # first_sulfamide_index = sul_index[::13]  
    # indices = np.where(np.diff(first_sulfamide_index) == 25)[0]
    # for index in range(len(indices)):
        
    #     last_sulfamide_index = first_sulfamide_index[indices[index] ]+12 # index already count the first element of the SUL residue
    #     # print(last_sulfamide_index)
    #     # print(Residue_atom_wise[last_sulfamide_index])
    #     # print(Residue_atom_wise[last_sulfamide_index+1:last_sulfamide_index+7])
    #     Residue_atom_wise[last_sulfamide_index+1:last_sulfamide_index+7] ="SCS"
    #     # print("-------------------")
    
    Residue_number = np.concatenate(Residue_number)  
    
    # h = open("molecule.xyz", "w")
    # h.write(str(len(atom_types)) +"\n")
    # h.write("Sulfamide \n")
    # for i in range(len(atom_types)):
    #     h.write(str(atom_types[i]) +"    " + str(float(coordinates[i][0])) +"    " + str(float(coordinates[i][1]))+"    " + str(float(coordinates[i][2])) +"\n")
    # h.close()
    
    monomer_dataframe =pd.DataFrame({'Atom_type': atom_types,'Element_type': element_types,'Residue':Residue_atom_wise ,'Residue_Number':Residue_number , 'x': list(coordinates[:,0]),'y': list(coordinates[:,1]),'z': list(coordinates[:,2])})
    
     
    monomer_dataframe =pd.DataFrame({'Atom_type': atom_types,'Element_type': element_types,'Residue':Residue_atom_wise ,'Residue_Number':Residue_number , 'x': list(coordinates[:,0]),'y': list(coordinates[:,1]),'z': list(coordinates[:,2])})
    monomer_data = np.array(monomer_dataframe)
    return monomer_dataframe


# %% Structure ---> more chains

shift_x = 0
nchains = 1
adding_coordinates = np.zeros((nchains, 3))
count = 0
break_outer_loop = False
for i in range(math.ceil(math.sqrt(nchains))):
    adding_coordinates[count] = [0,5*i,0]
    for j in range(nchains//int(math.sqrt(nchains))):
        adding_coordinates[count] = [shift_x*j,5*i,5*j]
        count += 1
        if count >= nchains:
            break_outer_loop = True
            break
    if break_outer_loop:
        break 
    
#%%        

#molecules_residue = ["CCF","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","SUL","CCA","CCL"]
molecules_residue = ["CCF","CCA"]+["SSF","SSL","CCA","CCA","CCA","CCA"]*20+["CCA","CCL"]
# molecules_residue[-1] = "CCL"
molecules_residue = np.array(molecules_residue)
positions_of_SUS = np.where(np.logical_or(molecules_residue == "SSF", molecules_residue == "SSL"))[0]
count_SUS = len(positions_of_SUS)


#%%
system = pd.DataFrame(columns = ['Atom_type', 'Element_type', 'Residue','Residue_Number', 'x', 'y','z'])
for i in range(nchains):
    chain=monomer(molecules_residue)
    
    chain.iloc[:,-3:] = chain.iloc[:,-3:] + adding_coordinates[i]
    # system.append(chain)
    system = pd.concat([system,chain])
system.reset_index(drop=True, inplace=True)

print("DONE-------------------")

# %% XYZ file for the system

f = open("System.xyz", "w")
f.write(str(len(system)) +"\n")
f.write("Sulfamide \n")
for i in range(len(system)):
    f.write(str(system["Atom_type"][i]) +"    " + str(float(system["x"][i])) +"    " + str(float(system["y"][i]))+"    " + str(float(system["z"][i])) +"\n")
f.close()
#%% Chains

import string
# Generate alphabet characters from 'A' to 'Z'
alphabet_upper = string.ascii_uppercase

# Create array starting from 'A' to 'Z'
array_from_A_to_Z = list(alphabet_upper)

alphabet_lower = string.ascii_lowercase

# Create array starting from 'A' to 'Z'
array_from_a_to_z = list(alphabet_lower)


number = np.linspace(1,9,9,dtype = int)

chain_numbers =np.concatenate((array_from_A_to_Z,array_from_a_to_z,number))


num_of_atoms_per_chain = int(len(system)/nchains)
chain_array = [num for num in chain_numbers for _ in range(num_of_atoms_per_chain)]    
# %%
system.iloc[:,-3:] = np.array(system.iloc[:,-3:], dtype=float)
coordinates = np.round(np.array(system.iloc[:,-3:]), 7)

Temperaturefactor = 0.00
g = open("2_10_chains_"+str(int(nchains))+"_sulf_"+str(count_SUS)+".pdb", "w")
# g = open("shift_"+str(int(shift_x))+"_2_10_chains_"+str(int(nchains))+"_sulf_"+str(count_SUS)+".pdb", "w")
g.write(str("COMPOUD") +"\n")

for i in range(len(system)):
    filled_string = list(" " * 80)
        
    # filled_string[0:6] ="HETATM"
    # filled_string[6:11] =str(i+1).rjust(4)
    # filled_string[12:16] =str(system["Atom_type"][i]).ljust(3)
    # filled_string[17:20] =str(system["Residue"][i]).rjust(3)
    # filled_string[21] = str(chain_array[i])
    # filled_string[22:26] =  str(system["Residue_Number"][i]).rjust(3)
    # filled_string[30:38] =  str(float(coordinates[i][0]))[:7:].rjust(7)
    # filled_string[38:46] =  str(float(coordinates[i][1]))[:7:].rjust(7)
    # filled_string[46:54] =  str(float(coordinates[i][2]))[:7:].rjust(7)
    # filled_string[54:60] =  str(1.00).rjust(5)
    # filled_string[60:66] =  str(Temperaturefactor).rjust(5)
    # filled_string[76:77] =  str(system["Element_type"][i]).rjust(1)  #Segment Identifier
    #filled_string[77:78] =  str(atom_type[i]).rjust(1)  # Element
    #filled_string[79:80] =  str(0).rjust(1)  # Charge
    
    
    filled_string[0:6] ="HETATM"
    filled_string[7:11] =str(i+1).rjust(4)
    filled_string[13:16] =str(system["Atom_type"][i]).ljust(3)
    filled_string[17:20] =str(system["Residue"][i]).rjust(3)
    filled_string[21] = str(chain_array[i])
    filled_string[22:26] =  str(system["Residue_Number"][i]).rjust(3)
    filled_string[31:38] =  str(float(coordinates[i][0]))[:7:].rjust(7)
    filled_string[39:46] =  str(float(coordinates[i][1]))[:7:].rjust(7)
    filled_string[47:54] =  str(float(coordinates[i][2]))[:7:].rjust(7)
    filled_string[55:60] =  str(1.00).rjust(5)
    filled_string[61:66] =  str(Temperaturefactor).rjust(5)
    filled_string[76:77] =  str(system["Element_type"][i]).rjust(1)
    # filled_string[73:76] =  str(element_types[i]).ljust(3)  #Segment Identifier
    # filled_string[77:78] =  str(element_types[i]).rjust(1)  # Element

    
    line = "".join(filled_string)
    g.write(str(line) +"\n")
g.close()

'''
#%% Connect

connect_atoms =pd.DataFrame({'Atom_type': atom_types,'Residue':Residue_atom_wise ,'Residue_Number':Residue_number,'Atom_Number':np.linspace(1,len(atom_types), len(atom_types))})

## Within each residue
# CCL

connect_CCL = {"C2": ["H3","H4","H5","C1"],
               "C1": ["H1","H2","C2"],
               "H1": ["C1"],
               "H2": ["C1"],
               "H3": ["C2"],
               "H4": ["C2"],
               "H5": ["C2"],
               }
connect_CCA = {"C2": ["H3","H4","C1"],
               "C1": ["H1","H2","C2"],
               "H1": ["C1"],
               "H2": ["C1"],
               "H3": ["C2"],
               "H4": ["C2"],
               }

connect_SUL = {"CB1": ["HB1","HB2","CA1"],
               "HB1":["CB1"],
               "HB2":["CB1"],
               "CB2": ["HB3","HB4","CA2"],
               "HB3":["CB2"],
               "HB4":["CB2"],
               "CA1": ["HA1","HA2","CB1","N1"],
               "HA1":["CA1"],
               "HA2":["CA1"],
               "CA2": ["HA3","HA4","CB2","N2"],
               "HA3":["CA2"],
               "HA4":["CA2"],
               "N1": ["NH1","S","CA1"],
               "NH1":["N1"],
               "N2": ["NH2","S","CA2"],
               "NH2":["N2"],
               "S":["N1","N2","O1","O2"],
               "O1":["S"],
               "O2":["S"]
               }

Interconnection = [["C1","C2"],
             ["C2","C1"],
             ["CB2","C1"],
             ["C1","CB2"],
             ["CB1","C2"],
             
             ["C2","CB1"]
             ]
connect = {i:[] for i in range(1,len(atom_types)+1)}
connects = []
num_atom = 0
# Connecting within the resiude
for i in range(len(molecules_residue)):
    connects.append(globals()["connect_"+molecules_residue[i]])
    k = globals()["connect_"+molecules_residue[i]]
    #for j in range(len(k)):
    array = connect_atoms[num_atom:num_atom+len(k)]
    num_atom += len(k)
    
    for j in range(len(list(k.keys()))):
        #y = array.loc[array['Atom_type'] == list(k.keys())[j], 'Atom_Number'].values[0]
        #print(y)
        for l in range(len(k[list(k.keys())[j]])):
            #y = array.loc[array['Atom_type'] == k[list(k.keys())[j]][l], 'Atom_Number'].values[0]
            connect[array.loc[array['Atom_type'] == list(k.keys())[j], 'Atom_Number'].values[0]].append(array.loc[array['Atom_type'] == k[list(k.keys())[j]][l], 'Atom_Number'].values[0])

num_atom = 0
connected_numbers =[]
#Between the residue
for i in range(len(molecules_residue)-1):
    current_residue = globals()["connect_"+molecules_residue[i]]
    current_residue_array = connect_atoms[num_atom:num_atom+len(current_residue)]
    num_atom += len(current_residue_array)
   
    next_residue = globals()["connect_"+molecules_residue[i+1]]
    next_residue_array = connect_atoms[num_atom:num_atom+len(next_residue)]
    #num_atom += len(next_residue_array)
    if (molecules_residue[i] == "CCL" and molecules_residue[i+1] == "CCA") or molecules_residue[i] == "CCL" and molecules_residue[i+1] == "CCL":
        first = current_residue_array.loc[current_residue_array['Atom_type'] == "C1", 'Atom_Number'].values[0]
        second =next_residue_array.loc[next_residue_array['Atom_type'] == "C1", 'Atom_Number'].values[0]
        print(first)
        print(current_residue_array)
        print(second)
        print(next_residue_array)
        
    if molecules_residue[i] == "CCA" and molecules_residue[i+1] == "CCA":
        first = current_residue_array.loc[current_residue_array['Atom_type'] == "C2", 'Atom_Number'].values[0]
        second = next_residue_array.loc[next_residue_array['Atom_type'] == "C1", 'Atom_Number'].values[0]
        print(first)
        print(current_residue_array)
        print(second)
        print(next_residue_array)  
        
    if molecules_residue[i] == "CCA" and molecules_residue[i+1] == "SUL":
        first = current_residue_array.loc[current_residue_array['Atom_type'] == "C2", 'Atom_Number'].values[0]
        second = next_residue_array.loc[next_residue_array['Atom_type'] == "CB1", 'Atom_Number'].values[0]

    if (molecules_residue[i] == "SUL" and molecules_residue[i+1] == "CCA") or (molecules_residue[i] == "SUL" and molecules_residue[i+1] == "CCL"):
        first = current_residue_array.loc[current_residue_array['Atom_type'] == "CB2", 'Atom_Number'].values[0]
        second = next_residue_array.loc[next_residue_array['Atom_type'] == "C1", 'Atom_Number'].values[0]
    
    connected_numbers.append([first,second])
    connect[first].append(second)
    connect[second].append(first)

number_bonds = [len(value) for value in connect.values()]
monomer["bonds"] = number_bonds

#%%

g = open("Structure_connect.pdb", "w")
g.write(str("COMPOUD") +"\n")

for i in range(len(atom_types)):
    filled_string = list(" " * 80)
    filled_string[0:6] ="HETATM"
    filled_string[7:11] =str(i+1).rjust(4)
    filled_string[13:16] =str(atom_types[i]).ljust(3)
    filled_string[18:20] =str(Residue_atom_wise[i]).rjust(3)
    filled_string[22] = str(1)
    filled_string[23:26] =  str(Residue_number[i]).rjust(3)
    filled_string[31:38] =  str(float(coordinates[i][0]))[:7:].rjust(7)
    filled_string[39:46] =  str(float(coordinates[i][1]))[:7:].rjust(7)
    filled_string[47:54] =  str(float(coordinates[i][2]))[:7:].rjust(7)
    filled_string[55:60] =  str(1.00).rjust(5)
    filled_string[61:66] =  str(Temperaturefactor).rjust(5)
    filled_string[73:76] =  str(element_types[i]).ljust(3)  #Segment Identifier
    filled_string[77:78] =  str(element_types[i]).rjust(1)  # Element
    #filled_string[79:80] =  str(0).rjust(1)  # Charge
    line = "".join(filled_string)
    g.write(str(line) +"\n")
    
for i in range(len(atom_types)):    
    connect_string = list(" " * 80)
    connect_string[0:6] ="CONECT"
    connect_string[7:11] =str(int((list(connect.keys())[i]))).rjust(4)
    for l in range(len(connect[list(connect.keys())[i]])):
        print(l)
        connect_string[12+(5*l):16+(5*l)] =str(int(connect[list(connect.keys())[i]][l])).rjust(4)
    line2 = "".join(connect_string)
    g.write(str(line2) +"\n")
g.close() 
'''# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 09:49:52 2024

@author: jayshah
"""

