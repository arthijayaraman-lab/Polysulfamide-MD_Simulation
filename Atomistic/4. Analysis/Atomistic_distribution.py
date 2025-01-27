
"""
@author: jayshah
"""

import mdtraj as md
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
import matplotlib.ticker as ticker
import seaborn as sns

def pbc_distance(r1,r2, len_box):
    """
    Calculate the periodic boundary condition (PBC) distance between two sets of points in a simulation box.
    
    Parameters:
        r1 (array): Array of coordinates for the first set of points (shape: [frames, particles, 3]).
        r2 (array): Array of coordinates for the second set of points (shape: [frames, particles, 3]).
        len_box (array): Array of box dimensions (shape: [frames, 3]) representing the size of the simulation box.
    
    Returns:
        list: List of distances between corresponding points from r1 and r2, accounting for periodic boundary conditions.
    """
    frame = np.shape(r1)[0]
    distance =[]
    for i in range(frame):
        box_len = len_box[i,:]
        
        r1_cor = np.array(r1[i,:])
        r2_cor = np.array(r2[i,:])
        
        delta = r1_cor - r2_cor
        print(np.shape(r1_cor))
        print(np.shape(r2_cor))
        print(np.shape(delta))
        distance.append(np.linalg.norm((delta)-box_len*np.round((delta)/box_len),axis = 1))
    return distance   

def find_middle_number(arr):
    """
    Find the middle number (median) in an array.
    
    Parameters:
        arr (array-like): Input array of numbers.
    
    Returns:
        float: The middle number in the array. If the array has an odd length, returns the middle element.
                If the array has an even length, returns the average of the two middle elements.
    """
    sorted_arr = np.sort(arr)
    n = len(sorted_arr)
    
    if n % 2 == 1:
        # Array has odd length
        middle_number = sorted_arr[n // 2]
    else:
        # Array has even length
        middle_number = (sorted_arr[n // 2 - 1] + sorted_arr[n // 2]) / 2
    
    return middle_number
#%% Loading the trajectories

temp = 350
trial = 3
nchains = 1
num_sulfamide_chain =10 # manually depending on the initial structure 
shift = 0
Solvent = "Water"

if Solvent == "Water":  
    files_dir = 'F:/Jay/Training_MD/Atomistic/Gromacs/4_8/T_'+str(temp)+"/"
    
    
    if shift == 0:
        file=files_dir +str(nchains) +"_chains_"+str(num_sulfamide_chain) +"_sulf_"+str(trial)+"/"
    else:
        file=files_dir +str(nchains) +"_chains_"+str(num_sulfamide_chain) +"_sulf_"+str(trial)+"_shift_"+str(int(shift))+"/"    
    xtc =  "md_prod_"+"com_unwrap.xtc" 
    top =  "md_prod_unwrapped"+"_com_unwrap.gro" 
    traj = md.load(file+xtc, top=file+top)
else:
    files_dir = 'F:/Jay/Training_MD/Atomistic/Gromacs/Solvent/'
    solvent_file = Solvent +'_Solvent/T_'+str(temp)+'/'
    file = files_dir + solvent_file + "4_8_chains_" +str(nchains) +"_sulf_"+str(num_sulfamide_chain) +"_"+str(trial)+"/"
    traj = md.load(file +"md_prod_unwrapped_"+"com_unwrap.gro")   

#%% Loading the unwrap co-ordinate and elements


topology = traj.topology
nframes = traj.n_frames
positions = traj.xyz
time = traj.time
len_box = traj.unitcell_lengths

S = topology.select("name == S")
O1 = topology.select("name == O1")
O2 = topology.select("name == O2")
N1 = topology.select("name == N1")
N2 = topology.select("name == N2")
NH1 = topology.select("name == NH1")
NH2 = topology.select("name == NH2")
CA1 = topology.select("name == CA1")
CA2 = topology.select("name == CA2")
CB1 = topology.select("name == CB1")
CB2 = topology.select("name == CB2")
COM = topology.select("name == COM")
CN1 = topology.select("name == CN1")
CN2 = topology.select("name == CN2")
C1 = topology.select("name == C1")
C2 = topology.select("name == C2")
print(traj.topology.n_atoms)

#%%

water = topology.select("resname == HOH")
number_atoms = traj.topology.n_atoms - len(water)
chain_length = number_atoms // nchains

print("number of atoms/species in a chain %s" %chain_length)

#%% Configurations of NH groups in sulfamide --> Dihedral between two NH groups
# This will help in deciding the positions of donors

atom_pairs = list(zip(N1, NH1))
atom_pairs2 = list(zip(N2,NH2))


distance = md.compute_displacements(traj, atom_pairs)   # Compute the displacement vector between pairs of atoms in each frame of a trajectory.
distance2 = md.compute_displacements(traj, atom_pairs2) # Compute the displacement vector between pairs of atoms in each frame of a trajectory.

vectors_array1 = np.array(distance)
vectors_array2 = np.array(distance2)

# Calculate the dot product of each pair of vectors
dot_products = np.sum(vectors_array1 * vectors_array2, axis=2)

# Calculate the magnitudes of each vector
magnitude1 = np.linalg.norm(vectors_array1, axis=2)
magnitude2 = np.linalg.norm(vectors_array2, axis=2)

# Calculate the cosine of the angle between each pair of vectors
cosine_angles = dot_products / (magnitude1 * magnitude2)

# Calculate the angles in radians
angles_radians = np.arccos(np.clip(cosine_angles, -1.0, 1.0))

# Convert angles to degrees if needed
angles_degrees = np.degrees(angles_radians)

# Display the result
print("Angles between vectors (degrees):", angles_degrees)


fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()
# hist, bin_edges = np.histogram(angles_degrees.flatten(), bins=50, density=True)   
# plt.plot(bin_edges[1:], hist/np.sum(hist), label="Configuration", color='blue', marker='o',markersize=2)
# plt.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)

sns.histplot(data=angles_degrees.flatten(), stat = "density", kde=True,label="dihedral angle",color = "#4687AC",edgecolor='#224153',bins = 30)

plt.legend(fontsize="18",loc ="upper left")

ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Angle',fontsize=18) 

file_path = file + 'angles_data.csv'

import csv
with open(file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for angle in angles_degrees.flatten():
        writer.writerow([angle])
csvfile.close()

plt.title("Angle distribution between NH groups",fontsize=18)
#%% End to end distance --> as it's in wrap coordinated, not using the PBC conditions. 

C1 = topology.select("name == C1 and resname == CCF")
C2 = topology.select("name == C2 and resname == CCL")

split = len(C1)/nchains

first_atom = np.array_split(C1, split)[0]
last_atom = np.array_split(C2, split)[-1]

first_atom_positions = traj.xyz[:, first_atom, :]
last_atom_positions = traj.xyz[:, last_atom, :]


Ree = pbc_distance(first_atom_positions,last_atom_positions, len_box)

fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()
ree_chains= []
for i in range(np.shape(Ree)[0]):
    ree_chains.append(Ree[i])

ax.set_ylabel('End to end distance',fontsize=18)
ax.set_xlabel('Timesteps',fontsize=18) 

plt.plot(ree_chains) 
if shift == 0: 
    plt.title('End to end distance')
else: 
    plt.title(shift)
    
data = {
    'Ree': np.concatenate(ree_chains)
}    
df = pd.DataFrame.from_dict(data, orient='index').transpose()
df.to_csv(file+"end_end.csv", index=False)
#%% Radius of gyration

Radius_of_gyration = md.compute_rg(traj, masses=None)
plt.plot(Radius_of_gyration)
 
if shift == 0: 
    plt.title('Radius of Gyration')
else: 
    plt.title(shift)

data = {
    'Rg': Radius_of_gyration
}
df = pd.DataFrame.from_dict(data, orient='index').transpose()
df.to_csv(file+"radius_gyration.csv", index=False)
    

#%% MSE MEAN SQUARE END TO END DISTANCE R^2 - for the entire chain 
C1 = topology.select("name == C1")
C2 = topology.select("name == C2")
## Need to correct !!!

backbone = np.hstack((S,CA1,CB1,CB2,CA2,C1,C2,N1,N2))
backbone = np.sort(backbone)
backbone[ [0, 1]] = backbone[ [1, 0]]


Cn = np.zeros((len(backbone)-1))
Rn = np.zeros((len(backbone)-1))
length =  np.zeros((len(backbone)-1))

for i in range((len(backbone)-1)):
    # print(int(backbone[i+1]))
    li = traj.xyz[:,int(backbone[i+1]), :] - traj.xyz[:,int(backbone[i]), :]
    length[i] = np.mean(np.linalg.norm(li,axis = 1))

l_avg = np.mean(length)    

for n in range(1,(len(backbone))):  
    r_2 =  np.zeros(len(backbone)-n)
    for i in range(len(backbone)-n): 

        backbone_atom_range = backbone[i:i+n+1]
        rn_ith = traj.xyz[:,int(backbone_atom_range[-1]), :] - traj.xyz[:,int(backbone_atom_range[0]), :]
        r_2[i] += np.mean(np.linalg.norm(rn_ith,axis = 1)**2)
             
    Rn[n-1] = np.mean(r_2)   
    Cn[n-1] = Rn[n-1]/(n*(l_avg**2))
    
plt.plot(Cn)
plt.xscale("log")         
plt.show()  
  
#%% RDF 

# COM = topology.select("name == COM")
COM_pairs = traj.top.select_pairs("name COM", "name COM")
# Compute the RDF
bins = 300
r_max = 4
r_min = 0.01
rdf_da, r_da = md.compute_rdf(traj, COM_pairs, (r_min, r_max), n_bins=bins)

rdf_da_df = pd.DataFrame({'Distance (nm)': rdf_da, 'g(r)': r_da})
rdf_da_df.to_csv(file+'rdf_da.csv', index=False)

# Plot the RDF
import matplotlib.pyplot as plt
plt.plot(rdf_da,r_da)
plt.plot(rdf_da,r_da, "o", label="mdtraj", alpha=0.1)
plt.xlabel('Distance (nm)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function between COMs')
plt.show()

donor_acceptor_pair = traj.top.select_pairs("name NH1" or "name NH2", "name O1" or "name O2")
filtered_pairs = [pair for pair in donor_acceptor_pair if abs(pair[0] - pair[1]) >= 5]

bins = 300
r_max = 4
r_min = 0.01
rdf, r = md.compute_rdf(traj, filtered_pairs, (r_min, r_max), n_bins=bins)


# Plot the RDF
import matplotlib.pyplot as plt
plt.plot(rdf,r)
plt.plot(rdf,r, "o", label="mdtraj", alpha=0.1)
plt.xlabel('Distance (nm)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function between donor-acceptor')
plt.show()


#%% Probability distirbution of distance
# Distance between the Adjacent Carbon's in the backbon
C1 = topology.select("name == C1")
C2 = topology.select("name == C2")
CA1 = topology.select("name == CA1")
CA2 = topology.select("name == CA2")
CB1 = topology.select("name == CB1")
CB2 = topology.select("name == CB2")
S = topology.select("name == S")
index = np.hstack((CA1,CB1,CB2,CA2,C1,C2))
index = np.sort(index) 
sulfur_index = np.sort(np.hstack((S)))
atom_pairs2 =[]
for i in range(len(index)-1):
    a = np.where(sulfur_index<index[i])
    b = np.where(sulfur_index<index[i+1])
    if (len(a[0]) == len(b[0])):
        atom_pairs2.append([index[i],index[i+1]])
        
atom_pairs2 = np.array(atom_pairs2, dtype=float)           
atom_pairs2[1]= [0.0,7.0] # CHANING BECAUSE THE GRO FILE HAVE THE FIRST C1 AND C2 in the revese way.  

distance_cc = md.compute_distances(traj, atom_pairs2)
distance_cc = np.array(distance_cc)


fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()

sns.histplot(data=distance_cc.flatten(), stat = "density", kde=True,label="Adj. Carbon in backbone",color = "#4687AC",edgecolor='#224153')

print("CG ratio of sulfamide bead: carbon bead %s" %(0.5/0.3))
#%%
atom_pairs3 = list(zip(COM,CA1))+list(zip(COM,CA2)) # Distance between the COM and C in the adjacent side
distance_sc = md.compute_distances(traj, atom_pairs3)
distance_sc = np.array(distance_sc)

sns.histplot(data=distance_sc.flatten(), stat = "density", kde=True,label="Sulfamide and carbon",color = "#CCAE06",edgecolor='#575216')

plt.legend(fontsize="18",loc ="upper left",bbox_to_anchor=(1.00, 1.05))
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02)) 
ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Distance',fontsize=18) 

plt.title("Distance distribution",fontsize=18)

#%% Distance - Bond Stretching potential

unit_in_kcal_mol = 6.9477*10**-21 #kcal/mol from Joule/molecule

fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()

kb = 1.380649* 10**-23 #joule per kelvin (K)

data = {
    'cc': distance_cc.flatten(),
    'sc': distance_sc.flatten()
}
df = pd.DataFrame.from_dict(data, orient='index').transpose()
df.to_csv(file+"distance_bond.csv", index=False)


#%% Probability distirbution of Angles  - configurations 

angle_indices=list(zip(CN2,COM,O1))
angle_indices_2=list(zip(CN2,COM,O2))

angle = md.compute_angles(traj, angle_indices)
angle2 = md.compute_angles(traj, angle_indices_2)

angle = np.array(np.rad2deg(angle)) # converting the angles to degress
angle2 = np.array(np.rad2deg(angle2))

fig = plt.figure()
fig.set_size_inches(6.00, 5.00, forward=True)
ax = fig.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) 


hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
hist2, bin_edges2 = np.histogram(angle2.flatten(), bins=50, density=True)

ax.plot(bin_edges[1:], hist/np.sum(hist), label="D2-COM-A1", color='purple', marker='o',markersize=2)
ax.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)

ax.plot(bin_edges2[1:], hist2/np.sum(hist2), label="D2-COM-A2", color='skyblue', marker='o',markersize=2)
ax.fill_between(bin_edges2[1:], 0, hist2/np.sum(hist2), color='skyblue', alpha=0.5)

plt.legend(fontsize="18",loc ="upper left")
print("number of angles %s" %(len(angle.flatten())))
plt.title("Angle distribution",fontsize=18)
ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) 
#ax.set_xlim(106,125)
ax.set_ylim(0, 0.15)
ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Angle',fontsize=18) 
fig.savefig('Angles.svg')


# Need to look at the configurations possible

angle_indices_d1_a1=list(zip(CN1,S,O1))
angle_indices_d1_a2=list(zip(CN1,S,O2))
angle_indices_d2_a1=list(zip(CN2,S,O1))
angle_indices_d2_a2=list(zip(CN2,S,O2))

angle_d1_a1 = np.array(np.rad2deg(md.compute_angles(traj, angle_indices_d1_a1)))
angle_d1_a2 = np.array(np.rad2deg(md.compute_angles(traj, angle_indices_d1_a2)))
angle_d2_a1 = np.array(np.rad2deg(md.compute_angles(traj, angle_indices_d2_a1)))
angle_d2_a2 = np.array(np.rad2deg(md.compute_angles(traj, angle_indices_d2_a2)))

data = {
    'angle_d2_a1': angle_d2_a1.flatten(),
    'angle_d2_a2': angle_d2_a2.flatten(),
    'angle_d1_a1': angle_d1_a1.flatten(),
    'angle_d1_a2': angle_d1_a2.flatten()
}

df = pd.DataFrame.from_dict(data, orient='index').transpose()
df.to_csv(file+"angle_d_com_a.csv", index=False)

angle_sum = angle_d1_a1+angle_d1_a2

fig = plt.figure()
fig.set_size_inches(6.00, 5.00, forward=True)
ax = fig.gca()

ax.xaxis.set_major_locator(ticker.MultipleLocator(10)) 
hist, bin_edges = np.histogram(angle_sum.flatten(), bins=50, density=True)
plt.plot(bin_edges[1:], hist/np.sum(hist), label="Sum of angles", color='purple', marker='o',markersize=2)
plt.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)
plt.legend(fontsize="18",loc ="upper left")
ax.xaxis.set_major_locator(ticker.MultipleLocator(10)) 
ax.set_ylim(0, 0.15)
ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Angle',fontsize=18)
fig.savefig('Angles_Sum.svg')
#

angle_d1_a1 = angle_d1_a1.flatten()
angle_d1_a2 = angle_d1_a2.flatten()
angle_d2_a1 = angle_d2_a1.flatten()
angle_d2_a2 = angle_d2_a2.flatten()

count_1 = 0
count_0 = 0
other_angles_d1_a1 = []
other_angles_d1_a2 = []
other_angles_d2_a1 = []
other_angles_d2_a2 = []
for i in range(len(angle_d1_a1)):
    if angle_d1_a1[i] < 100:
        if angle_d1_a2[i] > 220 - 100:
            other_angles_d1_a1.append(angle_d1_a1[i])
            other_angles_d1_a2.append(angle_d1_a1[i])
            other_angles_d2_a1.append(angle_d2_a1[i])
            other_angles_d2_a2.append(angle_d2_a2[i])
            count_1 += 1
        else:
            count_0 +=1
other_angles_d1_a1 = np.array(other_angles_d1_a1)
other_angles_d1_a2 = np.array(other_angles_d1_a2)
other_angles_d2_a1 = np.array(other_angles_d2_a1)
other_angles_d2_a2 = np.array(other_angles_d2_a2)  

hist, bin_edges = np.histogram(other_angles_d1_a1+other_angles_d2_a1, bins=50, density=True)
plt.plot(bin_edges[1:], hist/np.sum(hist), label="Approximation of the angle", color='purple', marker='o',markersize=2)
plt.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)    
    
#%% Backbone Angle Potential

index_angle  = np.hstack((CB1,CB2,C1,C2,CA1,CA2))
index_angle = np.sort(index_angle)
angle_indices_backbone = np.zeros((np.shape(index_angle)[0]-2,3))

sulfur_index = np.sort(np.hstack((S)))

C_C_C =[]

for i in range(np.shape(index_angle)[0]-2):
    a = np.where(sulfur_index<index_angle[i])
    b = np.where(sulfur_index<index_angle[i+1])
    c = np.where(sulfur_index<index_angle[i+2])
    if ((len(a[0]) == len(b[0])) and (len(c[0]) == len(b[0]))):
        C_C_C.append([index_angle[i],index_angle[i+1],index_angle[i+2]])
    
# need to change the first carbon C2 = 0, C1 = 3, C2 = 7 so need [3,0,7] INSTEAD of [0,3,7]        
C_C_C[0] = [3.0,0.0,7.0]
C_C_C[1] = [0.0,7.0,10.0]

print("number of C-C-C %s" %len(C_C_C))
fig = plt.figure()
fig.set_size_inches(6.00, 5.00, forward=True)
ax = fig.gca()

angle_C_C_C  = md.compute_angles(traj, C_C_C)
angle_C_C_C = np.array(np.rad2deg(angle_C_C_C))
hist_C_C_C, bin_edges_C_C_C = np.histogram(angle_C_C_C.flatten(), bins=100, density=True)

bin_center_C_C_C = 0.5 * (bin_edges_C_C_C[1:] + bin_edges_C_C_C[:-1])
hist_C_C_C = hist_C_C_C/np.sum(hist_C_C_C)
V_B_C_C_C = -kb * temp *np.log(hist_C_C_C/(np.sin(np.deg2rad((bin_center_C_C_C)))))/unit_in_kcal_mol
ax.scatter(bin_center_C_C_C,V_B_C_C_C,s=4,marker="d",color = "#714D8E", label = "C-C-C")

C_S_C = list(zip(CB1,COM,CB2))
print("number of C-S-C %s" %len(C_S_C))

angle_C_S_C   = md.compute_angles(traj, C_S_C )
angle_C_S_C  = np.array(np.rad2deg(angle_C_S_C ))
hist_C_S_C, bin_edges_C_S_C = np.histogram(angle_C_S_C.flatten(), bins=25, density=True)
bin_center_C_S_C = 0.5 * (bin_edges_C_S_C[1:] + bin_edges_C_S_C[:-1])
hist_C_S_C = hist_C_S_C/np.sum(hist_C_S_C)
V_B_C_S_C = -kb * temp *np.log(hist_C_S_C/(np.sin(np.deg2rad(bin_center_C_S_C))))/unit_in_kcal_mol
ax.scatter(bin_center_C_S_C,V_B_C_S_C,s=4,marker="d",color = "red", label = "C-S-C")


S_C_C= list(zip(COM,CB2,CB2+3))+list(zip(CB1-3,CB1,COM))
def remove_duplicates(input_list):
    return list(set(input_list))

S_C_C = remove_duplicates(S_C_C)

S_C_C_smaller_Carbon = []
S_C_C_more_Carbon =[]
for i in range(len(S_C_C)):
    if traj.topology.atom(S_C_C[i][1]).name == "CB2":
        print("yes")
        atom = traj.topology.atom(S_C_C[i][2])
        element = atom.name
        if element == "CB1":
            S_C_C_smaller_Carbon.append(S_C_C[i])
        if element == "C1":
            S_C_C_more_Carbon.append(S_C_C[i])
     
    if traj.topology.atom(S_C_C[i][1]).name == "CB1":
        print("yes")
        atom = traj.topology.atom(S_C_C[i][0])
        element = atom.name
        if element == "CB2":
            S_C_C_smaller_Carbon.append(S_C_C[i])
        if element == "C2":
            S_C_C_more_Carbon.append(S_C_C[i])

print("number of S-C-C - less carbon  %s" %len(S_C_C_smaller_Carbon)) 
print("number of S-C-C - bigger carbon chain  %s" %len(S_C_C_more_Carbon)) 
 
angle_S_C_C   = md.compute_angles(traj, S_C_C )
angle_S_C_C  = np.array(np.rad2deg(angle_S_C_C ))


hist_S_C_C, bin_edges_S_C_C = np.histogram(angle_S_C_C.flatten(), bins=50, density=True)
bin_center_S_C_C = 0.5 * (bin_edges_S_C_C[1:] + bin_edges_S_C_C[:-1])
hist_S_C_C = hist_S_C_C/np.sum(hist_S_C_C)
V_B_S_C_C = -kb * temp *np.log(hist_S_C_C/(np.sin(np.deg2rad(bin_center_S_C_C))))/unit_in_kcal_mol
ax.scatter(bin_center_S_C_C,V_B_S_C_C,s=5,marker="d",color = "black", label = "S-C-C")


print("C-C-C angle %s" %bin_center_C_C_C[np.argmax(hist_C_C_C)])
print("C-S-C angle %s" %bin_center_C_S_C[np.argmax(hist_C_S_C)])
print("S-C-C angle %s" %bin_center_S_C_C[np.argmax(hist_S_C_C)])

angle_S_C_C_smaller   = md.compute_angles(traj, S_C_C_smaller_Carbon )
angle_S_C_C_smaller  = np.array(np.rad2deg(angle_S_C_C_smaller ))

angle_S_C_C_more_chain   = md.compute_angles(traj, S_C_C_more_Carbon )
angle_S_C_C_more_chain  = np.array(np.rad2deg(angle_S_C_C_more_chain ))

data = {
    'ccc': angle_C_C_C.flatten(),
    'csc': angle_C_S_C.flatten(),
    'scc_s': angle_S_C_C_smaller.flatten(),
    'scc_l': angle_S_C_C_more_chain.flatten()

}
df = pd.DataFrame.from_dict(data, orient='index').transpose()
df.to_csv(file+"angle.csv", index=False)

plt.legend(fontsize="18",loc ="upper left")
plt.title("Angle Boltzmann Inversion",fontsize=18)
plt.show()
#%% S-C-C angle 
index = np.where(angle_S_C_C<120)
for i in range(np.shape(index)[1]):
    trajectory_number =index[0][i] 
    SCC_index = index[1][i] 
    print(S_C_C[SCC_index]) 
#%% Angle between acceptor - COM and the next carbon

angle_indices=list(zip(O2,COM,CA2))

angle = md.compute_angles(traj, angle_indices)

angle = np.array(np.rad2deg(angle)) # converting the angles to degress

fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50)) 


hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
hist2, bin_edges2 = np.histogram(angle2.flatten(), bins=50, density=True)

plt.plot(bin_edges[1:], hist/np.sum(hist), label="O1-COM-CA2", color='purple', marker='o',markersize=2)
plt.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)

plt.legend(fontsize="18",loc ="upper left")
print("number of angles %s" %(len(angle.flatten())))
plt.title("Angle distribution",fontsize=18)
ax.xaxis.set_major_locator(ticker.MultipleLocator(50)) 
ax.set_ylim(0, 0.15)
ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Angle',fontsize=18) 

#%% Probability distirbution of Torsion
 
dihedra_indices=list(zip(O1,COM,N1,NH1))
dihedral = md.compute_dihedrals(traj, dihedra_indices)

dihedral = np.array(np.rad2deg(dihedral)) +180 # converting the angles to degress
hist, bin_edges = np.histogram(dihedral.flatten(), bins=50, density=True)
plt.plot(bin_edges[1:], hist/np.sum(hist), label="CN1-S-CN2", color='purple', marker='o',markersize=2)
plt.title("Dihedral Distribution",fontsize=18)

#%% Backbone Angle
index = np.hstack((S,N1,N2,CA1,CB1,CB2,CA2,C1,C2))
index = np.sort(index)
angle_indices_backbone = np.zeros((np.shape(index)[0]-2,3))

sulfur_index = np.sort(np.hstack((S,N1,N2)))
sulfur_backbone = []
carbon_sulfur =[]
carbon_backbone = []

for i in range(np.shape(index)[0]-2):
    angle_indices_backbone[i]=[index[i],index[i+1],index[i+2]]
    a = np.in1d([index[i],index[i+1],index[i+2]], sulfur_index)

    if np.sum(a)== 3:
        sulfur_backbone.append([index[i],index[i+1],index[i+2]])
    elif np.sum(a)== 0:
        carbon_backbone.append([index[i],index[i+1],index[i+2]])
    else:
        carbon_sulfur.append([index[i],index[i+1],index[i+2]])
        

angle = md.compute_angles(traj, angle_indices_backbone)
angle = np.array(np.rad2deg(angle))

fig = plt.figure()
fig.set_size_inches(6.00, 5.00, forward=True)
ax = fig.gca()
hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
ax.plot(bin_edges[1:], hist/np.sum(hist), color='purple', marker='o',markersize=2)
ax.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)
plt.title("Backbone - Angle")
plt.show()

angle = md.compute_angles(traj, sulfur_backbone)
angle = np.array(np.rad2deg(angle))

fig = plt.figure()
fig.set_size_inches(6.00, 5.00, forward=True)
ax = fig.gca()
hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
ax.plot(bin_edges[1:], hist/np.sum(hist), label="sulfur_backbone",color='purple', marker='o',markersize=2)
ax.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)

angle = md.compute_angles(traj, carbon_backbone)
angle = np.array(np.rad2deg(angle))
hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
ax.plot(bin_edges[1:], hist/np.sum(hist), label="carbon_backbone",color='green', marker='o',markersize=2)
ax.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='green', alpha=0.3)

angle = md.compute_angles(traj, carbon_sulfur)
angle = np.array(np.rad2deg(angle))
hist, bin_edges = np.histogram(angle.flatten(), bins=50, density=True)
ax.plot(bin_edges[1:], hist/np.sum(hist), label="carbon_sulfur",color='blue', marker='o',markersize=2)
ax.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='blue', alpha=0.3)

plt.title("Backbone - Angle")
plt.legend()
plt.show()
#%% Probability distirbution of distance

atom_pairs = list(zip(COM,O1))
atom_pairs2 = list(zip(CA1,CB1))

distance = md.compute_distances(traj, atom_pairs)
distance2 = md.compute_distances(traj, atom_pairs2)
distance = np.array(distance)*2
distance2 = np.array(distance2)

fig = plt.figure()
fig.set_size_inches(5.00, 5.00, forward=True)
ax = fig.gca()

hist, bin_edges = np.histogram(distance.flatten(), bins=50, density=True)
hist2, bin_edges2 = np.histogram(distance2.flatten(), bins=50, density=True)

plt.plot(bin_edges[1:], hist/np.sum(hist), label="2*(COM-O1)", color='purple', marker='o',markersize=2)
plt.fill_between(bin_edges[1:], 0, hist/np.sum(hist), color='purple', alpha=0.3)

#sns.histplot(distance.flatten(),bins=100,stat="density",kde=True, color='skyblue', fill=True)

plt.plot(bin_edges2[1:], hist2/np.sum(hist2), label="(CA1-CB1)", color='skyblue', marker='o',markersize=2)
plt.fill_between(bin_edges2[1:], 0, hist2/np.sum(hist2), color='skyblue', alpha=0.5)

print("CG ratio of sulfamide bead: carbon bead %s" %(0.5/0.3))

print("Atomistic ratio of sulfamide bead: carbon bead %s" %np.average(find_middle_number(distance)/find_middle_number(distance2)))

atom_pairs3 = list(zip(COM,CA1))
distance3 = md.compute_distances(traj, atom_pairs3)
distance3 = np.array(distance3)

hist3, bin_edges3 = np.histogram(distance3.flatten(), bins=50, density=True)

plt.plot(bin_edges3[1:], hist3/np.sum(hist3), label="COM - CA1", color='green', marker='o',markersize=2)
plt.fill_between(bin_edges3[1:], 0, hist3/np.sum(hist3), color='green', alpha=0.5)

print("r1:r2 %s" %np.average(find_middle_number(distance)/find_middle_number(distance2)))
print("r2:r3 %s" %np.average(find_middle_number(distance2)/find_middle_number(distance3)))


plt.legend(fontsize="18",loc ="upper left")
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.02)) 
ax.set_ylabel('Probability',fontsize=18)
ax.set_xlabel('Distance',fontsize=18) 
print("number of distance %s" %(len(distance.flatten())))
print("number of distance %s" %(len(distance2.flatten())))
plt.title("Distance distribution",fontsize=18)
#%%  H-bonding 

print(md.baker_hubbard(traj,  freq=1, exclude_water=True, periodic=True, sidechain_only=False))
hydrogen = md.wernet_nilsson(traj)
for i in range(len(hydrogen)):
    if len(hydrogen[i])> 0:
        print("Yay")
        
array_acceptor = np.sort(np.append(O1,O2))
array_donor = np.sort(np.append(NH1,NH2))

valid_pairs = []
for i in array_acceptor:
    for j in array_donor:
        if abs(i-j) > 5:  
            valid_pairs.append((i, j))

distance = md.compute_distances(traj, valid_pairs)
n_hb = []

for i in range(0,np.shape(distance)[0],2):
    index = np.where(distance[i] < 0.35)
    if len(index) < 0:
        print("Yes!")
        hbonding_angle = np.zeros((len(index[0]),3))
        pairs_acceptor = []
        pairs_donor = []
        count_hb = 0
        for j in index[0]:
            #print(valid_pairs[j])
            pairs_acceptor.append(valid_pairs[j][0])
            pairs_donor.append(valid_pairs[j][1])
            hbonding_angle[count_hb] = [int(valid_pairs[j][0]),int(valid_pairs[j][1]),int(valid_pairs[j][1]-1)]
            count_hb += 1
        angles = md.compute_angles(traj[i:i+1], list(hbonding_angle))    
        n_hb.append(len(np.where(np.rad2deg(angles[0]) > 120)[0]))
    else:
        print("No!")
    #print("------")
    plt.plot(np.rad2deg(angles[0]))
    #plt.scatter(pairs_acceptor,pairs_donor)
plt.show()
sns.histplot(data=distance.flatten(), stat = "density", kde=True)