"""
@author: jayshah

"""

import math 

def file_name(nRU, nChain, aliL, aliR, seed, current_eps,trial,type_pos):
    file_path = ""
    file_pos_type = type_pos
    file = "/%s_%s/matrix_disulf_rampAuto_nRU%s_nChain%s_leftsataliph%s_rightsataliph%s__HBeps%s_%s.cube" %(aliL, aliR,nRU, nChain, aliL, aliR, current_eps,seed)
    return file_path+file_pos_type+file

def load_trajectories(nRU, nChain, aliL, aliR, seed,eps_start,eps_end,freq_eps,trial,type_pos):
    total = aliL+ aliR
    files=[]
    file_path = ""
    file_pos_type = type_pos
    for i in range(freq_eps+1):
        current_eps = round(eps_start + ((eps_end-eps_start)/freq_eps)*i,1)
        if (math.floor(current_eps) == current_eps):
            current_eps = int(current_eps)
            files.append(file_path+file_pos_type+"/%s_%s/matrix_disulf_rampAuto_nRU%s_nChain%s_leftsataliph%s_rightsataliph%s__HBeps%s_%s.cube" %(aliL, aliR,nRU, nChain, aliL, aliR, current_eps,seed))
        
        else:
            files.append(file_path+file_pos_type+"/%s_%s/matrix_disulf_rampAuto_nRU%s_nChain%s_leftsataliph%s_rightsataliph%s__HBeps%s_%s.cube" %(aliL, aliR,nRU, nChain, aliL, aliR, current_eps,seed))
    return files
