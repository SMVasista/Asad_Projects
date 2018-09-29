from __future__ import division
import os, re, sys, math, random, copy
import numpy as NP
from glob import glob
import geometrix as gx
import bioprop as bp
import pdbutils as pdb
import cwgutils

def initializeDataDump(name):
    #Initializing project simulation dump location
    dataLoc = os.path.join(os.getcwd(), str(name))
    #cwgutils.mkpdir(dataLoc)
    return dataLoc
    
def initializeEDF(EDFdict):
    if len(EDFdict.keys()) == 0:
        EDFdict_ = {'Total': 0.0, 'Stearic': 0.0, 'Elec_v': 0.0, 'Hyd_phi': 0.0, 'LJP': 0.0, 'Bondconst': 0.0, 'Hbond': 0.0, 'DSbond': 0.0, 'dx': 0.0, 'dy': 0.0, 'dz': 0.0}
        return EDFdict_
    else:
        return EDFdict

def minimizeEnergyStruc(projName, coOrdData):
    #Capturing data
    COORD_DATA_HOLDER = []
    
    #dataLoc = initializeDataDump(projName)
    
    COORD_DATA_HOLDER.append(coOrdData)
    
    #Initializing optimization
    t = 100
    k_temp = 300
    pH = 6.5
    rTemp = (k_temp)/300
    h = 0.05
    k_stearic = 2
    k_elec = 2
    k_hydp = 2
    k_ljp = 0.01
    k_bondcnst = -2
    k_hbond = 0.1
    k_dsbond = 0.75
    
    EDF = {}
    
    for i in range(t):
        print "Time", i,"ns"
        
        EDF[i] = {}
        print EDF[i]
    
        data = copy.deepcopy(COORD_DATA_HOLDER[-1])
        
        #Capturing AA sequence
        aaSeq = [(j[0], j[1]) for j in data]
        
        for n, aa in enumerate(aaSeq):
            
            aaTag = str(aa[0])+'_'+str(aa[1])
            EDF[i][aaTag] = {}
            EDF[i][aaTag] = initializeEDF(EDF[i][aaTag])
            
            aa_ = bp.map_aa_names(aa[0], '1')
            vol = bp.VDW_vol[aa_]
            eCord = gx.getCoOrdinates(data, aa[1])
            
            #Identifying amino acids in vicinity
            sub_data = []
            alpha = 0.02
            while len(sub_data) < 3:
                sub_data = gx.extractSphere(data, str(aa[1]), alpha*50)
                alpha += 0.01
            
            #Iterating over region to optimize locations of each AA in region
            for _aa in sub_data:
                if _aa[0] != aa[0] and _aa[1] != aa[1]:
                    _aa_ = bp.map_aa_names(_aa[0], '1')
                    _vol = bp.VDW_vol[_aa_]
                    _eCord = gx.getCoOrdinates(data, _aa[1])
                    
                    #Getting distance of sub-region entity from main AA
                    dist = gx.eucDist(eCord, _eCord)
                    
                    #Calculating stearic force
                    vol_append_dist = (vol**0.33)*0.2388 + (_vol**0.33)*0.2388
                    F_stearic = k_stearic*(vol_append_dist/dist)**2
                    
                    #Calculating LJ Potential Derived force
                    F_LJP = k_ljp*((3.4/dist)**12 - (3.4/dist)**6)
                    
                    #Calculating charge attraction force
                    F_elec = k_elec*((bp.Charge[aa_]*(bp.pKa[aa_] - pH))+0.001)*((bp.Charge[_aa_]*((bp.pKa[_aa_] - pH)))+0.001)/(dist**2)
                    
                    #Calculating hydrophobic interaction force
                    F_hydp = k_hydp*(bp.Hydrophobicity[aa_] - NP.mean(bp.Hydrophobicity.values()))*(bp.Hydrophobicity[_aa_] - NP.mean(bp.Hydrophobicity.values()))/NP.mean(bp.Hydrophobicity.values())**2
                    
                    #Calculating bond-constraint on Ca
                    F_bondconst = 0
                    try:
                        if _aa[1] == aaSeq[n-1][1] or _aa[1] == ssSeq[n+1][1]:
                            F_bondconst = k_bondcnst*(dist - 1.5)**2
                    except:
                        pass
                        
                    #Calculating Hydrogen Bonding potential
                    F_hbond = (dist - 1.5)**2 * k_hbond*0.5*(bp.HBond_donor_prop[aa_]*bp.HBond_acceptor_prop[_aa_] + bp.HBond_donor_prop[_aa_]*bp.HBond_acceptor_prop[aa_])
                    
                    #Calculating di-sulphide bonding
                    F_dsbond = 0
                    if aa_ == 'C' and _aa_ == 'C':
                        F_dsbond = k_dsbond*(dist - 1.7)**2

                    #T is the total force experienced by sub-region element. This represents the magnitude of force vector
                    T = (F_stearic + F_elec + F_hydp + F_LJP + F_bondconst + F_hbond + F_dsbond)
                    
                    #Resolving direction vector
                    O = eCord
                    A = _eCord
                                       
                    ##Force experienced by regional-element
                    OA = [A[0] - O[0], A[1] - O[1], A[2] - O[2]]
                    unit_OA = float(1/(math.sqrt(OA[0]**2 + OA[1]**2 + OA[2]**2)))
                    alphaA = math.acos(OA[0]*unit_OA)
                    betaA = math.acos(OA[1]*unit_OA)
                    gammaA = math.acos(OA[2]*unit_OA)
                    dFx = T*math.cos(alphaA)*h
                    dFy = T*math.cos(betaA)*h
                    dFz = T*math.cos(gammaA)*h
                    print aa, "-> (", _aa[0], ',', _aa[1], ") :", dFx, dFy, dFz
                    
                    #Recording data
                    EDF[i][aaTag]['Total'] += T
                    EDF[i][aaTag]['Stearic'] += F_stearic
                    EDF[i][aaTag]['Elec_v'] += F_elec
                    EDF[i][aaTag]['Hyd_phi'] += F_hydp
                    EDF[i][aaTag]['LJP'] += F_LJP
                    EDF[i][aaTag]['Bondconst'] += F_bondconst
                    EDF[i][aaTag]['Hbond'] += F_hbond
                    EDF[i][aaTag]['DSbond'] += F_dsbond
                    EDF[i][aaTag]['dx'] += (bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFx
                    EDF[i][aaTag]['dy'] += (bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFy
	            EDF[i][aaTag]['dz'] += (bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFz
	                
	                #Modifying co-ordinates of regional-element and original element
	            for elem in data:
	                if elem[0] == bp.map_aa_names(_aa_,'3') and elem[1] == _aa[1]:
	                    elem[2] = elem[2] + 0.5*(bp.Flexibility[_aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFx + (rTemp**5)*0.01*(random.random() - 0.5)
	                    elem[3] = elem[3] + 0.5*(bp.Flexibility[_aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFy + (rTemp**5)*0.01*(random.random() - 0.5)
	                    elem[4] = elem[4] + 0.5*(bp.Flexibility[_aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFz + (rTemp**5)*0.01*(random.random() - 0.5)
	                if elem[0] == bp.map_aa_names(aa_,3) and elem[1] == aa[1]:
	                    elem[2] = elem[2] - 0.5*(bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFx + (rTemp**5)*0.01*(random.random() - 0.5)
	                    elem[3] = elem[3] - 0.5*(bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFy + (rTemp**5)*0.01*(random.random() - 0.5)
	                    elem[4] = elem[4] - 0.5*(bp.Flexibility[aa_]/(bp.Flexibility[_aa_]+bp.Flexibility[aa_]))*dFz + (rTemp**5)*0.01*(random.random() - 0.5)

        COORD_DATA_HOLDER.append(data)
    return COORD_DATA_HOLDER, EDF
    
#def tbdSimulation():

def calculateEnergyFunc(COORD_DATA_HOLDER, EDF, pockets):
    
    
