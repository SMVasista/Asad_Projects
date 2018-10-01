from __future__ import division
import os, re, sys, math, random, copy
import numpy as NP
from glob import glob
import geometrix as gx
import bioprop as bp
import pdbutils as pdb
import cwgutils

cwl = os.getcwd()
k_temp = 300
pH = 6.5

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
    t = 25
    rTemp = (k_temp)/300
    h = 0.07
    k_stearic = 2
    k_elec = 2
    k_hydp = -2
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
                sub_data = gx.extractSphere(data, aa[1], alpha*50)
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
                    F_elec = k_elec*((bp.Charge[aa_]+(bp.pKa[aa_] - pH))+0.001)*((bp.Charge[_aa_]+((bp.pKa[_aa_] - pH)))+0.001)/(dist**2)
                    
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

def calculateEnergyFunc(pockets, COORD_DATA_HOLDER, EDF, outLoc):
    '''This program computes the time dependent energy functions of whole body (protein) and pockets'''
    
    dataC = copy.deepcopy(COORD_DATA_HOLDER)
    dataE = copy.deepcopy(EDF)
    
    #Calculating stability parameters
    inst = len(dataC)
    seq = [x[0] for x in dataC[0]]
    
    #1) Calculating RMS profile    
    RMSF = []
    init = dataC[0]
    
    for i in range(inst):
        rsq = 0
        for j in range(len(seq)):
            rsq += ((dataC[i][j][2] - dataC[0][j][2])**2 + (dataC[i][j][3] - dataC[0][j][3])**2 + (dataC[i][j][4] - dataC[0][j][4])**2)*(1/len(dataC))
        RMSF.append(math.sqrt(rsq))
        
    #2) Calculating energy minimization profile
    Emin = [0]
    for i in dataE.keys():
        etotal = 0
        for j in dataE[i].keys():
            etotal += dataE[i][j]['Total']*(1/len(dataC))
        Emin.append(etotal)
        
    #3) Calculating Gyration profile
    rGyration = []
    for i in range(inst):
        x_min = min([x[2] for x in dataC[i]])
        x_max = max([x[2] for x in dataC[i]])
        y_min = min([y[3] for y in dataC[i]])
        y_max = max([y[3] for y in dataC[i]])
        z_min = min([z[4] for z in dataC[i]])
        z_max = max([z[4] for z in dataC[i]])
        p_center = [(x_min+x_max)/2, (y_min+y_max)/2, (z_min+z_max)/2]
        rGyration.append(p_center)
        
    #4) Calculating H-Bond formation profile
    HBondCount = [0]
    for i in dataE:
        hbp = 0
        for j in dataE[i]:
            if abs(dataE[i][j]['Hbond']) > 0:
                hbp += 1
        HBondCount.append(hbp)
        
    HBondStrength = [0]
    for i in dataE:
        hbs = 0
        for j in dataE[i]:            
            hbs += abs(dataE[i][j]['Hbond'])
        HBondStrength.append(hbs)
    
    #Finished calculating stability parameters. Proceeding to calculate pocket-wise activity profile
    
    for n, loc in enumerate(pockets):
        pkID = 'PKT_'+str(n+1)
        
        apVol = (loc[0][1] - loc[0][0])*(loc[1][1] - loc[1][0])*(loc[2][1] - loc[2][0])
        centroid = [(loc[0][1] + loc[0][0])/2, (loc[1][1] + loc[1][0])/2, (loc[2][1] + loc[2][0])/2]
        
        _seq = []
        alpha = 0.02
        while len(_seq) < 3:
            _seq = gx.extractSphereCentroid(dataC[-1], centroid, alpha*50)
            alpha += 0.01
            
        #Calculating volume encroachment
        PKTencr = {}
        udataC = dataC[-1]
        udataE = dataE[sorted(dataE.keys())[len(dataE.keys())-1]]
        
        dvol = 0
        for aa in _seq:
            aa_ = bp.map_aa_names(aa[0], '1')
            ld = gx.eucDist(centroid, gx.getCoOrdinates(udataC, aa[1])) - ((bp.VDW_vol[aa_])**0.33)*0.2388
            sd = ld - gx.eucDist(centroid, [loc[0][0], loc[1][0], loc[2][0]])
            dvol -= sd*0.5*bp.VDW_vol[aa_]
        PKTencr[pkID] = apVol - dvol
        
        #Calculating Charge Density
        PKTchrd = {}
        
        chard = 0.0
        for aa in _seq:
            aa_ = bp.map_aa_names(aa[0], '1')
            chard += ((bp.Charge[aa_] + (bp.pKa[aa_] - pH))+0.001)* (1/gx.eucDist([aa[2], aa[3], aa[4]], centroid))**2 *(1/apVol)
        PKTchrd[pkID] = chard
        
        #Calculating Hydrophobicity
        PKThydp = {}
        
        hyphi = 0.0
        for aa in _seq:
            aa_ = bp.map_aa_names(aa[0], '1')
            hyphi += ((bp.Hydrophobicity[aa_] - NP.mean(bp.Hydrophobicity.values())))/NP.mean(bp.Hydrophobicity.values())*(1/apVol)*(1/len(_seq))
        PKThydp[pkID] = hyphi
        
        #Calculating arbitrary internal bonding energy
        PKTbeny = {}
        
        en = 0.0
        for aa in _seq:
            en += udataE[str(aa[0])+'_'+str(aa[1])]['Total'] * (1/gx.eucDist([aa[2], aa[3], aa[4]], centroid))**2 * (1/len(_seq)) * (1/apVol)
        PKTbeny[pkID] = en
        
        #Tabulating report
        
        with open(os.path.join(outLoc, 'catalytic_pocket_parameters.csv'), 'a') as f:
            for ids in PKTencr.keys():
                f.write(str(ids)+'\n')
                f.write('Centroid Location:'+str(centroid)+'\n')
                f.write('Volume Encroachment:'+str(PKTencr[ids])+'\n')
                f.write('Charge Density:'+str(PKTchrd[ids])+'\n')
                f.write('Pocket Hydrophobicity:'+str(PKThydp[ids])+'\n')
                f.write('Internal Energy:'+str(PKTbeny[ids])+'\n')
                f.write('\n')

    print "Finished calculating energy parameters, tabulating..."

    with open(os.path.join(outLoc, 'stability_parameters.csv'), 'a') as f:
    
        f.write('Time RMS Energy_Minimization HBond_Count HBond_profile\n')
        for n in range(len(RMSF)):
            f.write(str(n)+' '+str(RMSF[n])+' '+str(Emin[n])+' '+str(HBondCount[n])+' '+str(HBondStrength[n])+'\n')
    
    
    
    
