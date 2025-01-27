# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:27:16 2023

@author: jayshah
"""
import math 

def file_name(nRU, nChain1,nChain2, aliL1,aliL2, aliR1,aliR2, seed, current_eps,trial,last,type_pos):
    total = aliL1+ aliR1
    file = "/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Blends/BiModal/%s/%s_mers/%s/matrix_disulf_rampAuto_nChain_%s_%s_leftsataliph_%s_%s_rightsataliph_%s_%s_HBeps%s_%s.cube" %(type_pos,nRU,last,nChain1,nChain2, aliL1, aliL2,aliR1, aliR2, current_eps,seed)
    return file

def load_trajectories(nRU, nChain1,nChain2, aliL1,aliL2, aliR1,aliR2, seed,eps_start,eps_end,freq_eps,trial,last,type_pos):
    total = aliL1+ aliR1
    files=[]
    for i in range(freq_eps+1):
        current_eps = round(eps_start + ((eps_end-eps_start)/freq_eps)*i,1)
        if (math.floor(current_eps) == current_eps):
            current_eps = int(current_eps)
            files.append("/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Blends/BiModal/%s/%s_mers/%s/matrix_disulf_rampAuto_nChain_%s_%s_leftsataliph_%s_%s_rightsataliph_%s_%s_HBeps%s_%s.cube" %(type_pos,nRU,last,nChain1,nChain2, aliL1, aliL2,aliR1, aliR2, current_eps,seed))
        
        else:
            files.append("/ocean/projects/dmr200035p/jshah3/polysulfamide/Modifying_CG/Donor_Acceptor_position/Blends/BiModal/%s/%s_mers/%s/matrix_disulf_rampAuto_nChain_%s_%s_leftsataliph_%s_%s_rightsataliph_%s_%s_HBeps%s_%s.cube" %(type_pos,nRU,last,nChain1,nChain2, aliL1, aliL2,aliR1, aliR2, current_eps,seed))
    return files
