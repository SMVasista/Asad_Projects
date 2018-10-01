from __future__ import division
import os, re, sys, math, random, copy, pickle
import numpy as NP
from glob import glob
import cwgutils
import pdbutils as pdb
import datutils as dutils
import bioprop as bp
import geometrix as gx
import energetix as ex
import xyzviz as viz

cwl = os.getcwd()

def parseInput(fileLoc):
    #Init
    projName = None
    pdbLoc = None
    mutations = []
    
    data = cwgutils.readLines(fileLoc)
    
    for line in data:
        if list(line)[0] != '#': #Skip commented lines
            if 'gene_name' in line:
                projName = str(line.split(':')[1])
            elif 'pdbLoc' in line or 'PDBloc' in line:
                pdbLoc = str(line.split(':')[1])
            elif 'mutation' in line:
                mutations.append(parseMutation(line.split(':')[1]))
    return projName, pdbLoc, mutations

def parseMutation(mutaString):
    try:
        a = list(mutaString)
        start = a.pop(0)
        end = a.pop(-1)
        num = int(''.join(a))
        return (start, num, end)
    except:
        raise IOError
        return None
        
#def queryMutation(gene, mutation):
    #TODO
    
def generateProteinSequence(coOrdData, mutList):
    M_SEQUENCES = {}
    for sig in mutList:
        tdata = copy.deepcopy(coOrdData)
        start = bp.map_aa_names(sig[0], '3')
        end = bp.map_aa_names(sig[2], '3')
        for line in tdata:
            print line, start, end, sig[1]
            if line[1] == sig[1] and line[0] == start:
                print "Applying mutation:", sig, "on protein"
                line[0] = end
        M_SEQUENCES[sig] = tdata
    return M_SEQUENCES

if __name__=="__main__":
	
	script, inputfile = sys.argv
	
	projName, pdbLoc, mutList = parseInput(inputfile)
	
	H = pdb.readPDB(pdbLoc)
	uModel = pdb.stringPDB(H)
	aaCordData = gx.influenceSurface(uModel, uModel[0][3], cutOff=10)
	
	DBLoc = []
	
	DBLoc.append(str(os.path.join(cwl, projName)))
	
	##################################################
	#Applying Energy Minimization
	if os.path.exists(os.path.join(cwl, projName)) != True:
	    os.makedirs(os.path.join(cwl, projName))
	if os.path.exists(os.path.join(cwl, projName, 'traj')) != True:
	    os.makedirs(os.path.join(cwl, projName, 'traj'))
	
	C, E = ex.minimizeEnergyStruc(projName, aaCordData)
	
	with open(os.path.join(cwl, projName, 'sout.p'), 'w') as f:
	    pickle.dump(C,f)
	with open(os.path.join(cwl, projName, 'eout.p'), 'w') as f:
	    pickle.dump(E,f)
	
	pockets = gx.identifyCatalyticPockets(C[-1])
	#The pockets identified here are putative.
	
	ex.calculateEnergyFunc(pockets, C, E, os.path.join(cwl, projName))
	
	d, l = viz.extractSimulData(os.path.join(cwl, projName, 'sout.p'))
	viz.writeTrajectory(d, l, os.path.join(cwl, projName, 'traj'))
    
	##################################################
	#Applying mutation 
	#mutList = queryMutation(mutation)
	seq = generateProteinSequence(aaCordData, mutList)
    
	##################################################
	#Applying Energy Minimization (on mutations)
	for mut in mutList:
	    nName = str(projName)+'_'+str(mut[0])+str(mut[1])+str(mut[2])
	    DBLoc.append(str(os.path.join(cwl, nName)))
	    
	    if os.path.exists(os.path.join(cwl, nName)) != True:
	        os.makedirs(os.path.join(cwl, nName))
	    if os.path.exists(os.path.join(cwl, nName, 'traj')) != True:
	        os.makedirs(os.path.join(cwl, nName, 'traj'))
		
		C, E = ex.minimizeEnergyStruc(nName, seq[mut])
		
	    with open(os.path.join(cwl, nName, 'sout.p'), 'w') as f:
                pickle.dump(C,f)
	    with open(os.path.join(cwl, nName, 'eout.p'), 'w') as f:
                pickle.dump(E,f)
		
            ex.calculateEnergyFunc(pockets, C, E, os.path.join(cwl, nName))
        
            d, l = viz.extractSimulData(os.path.join(cwl, nName, 'sout.p'))
            viz.writeTrajectory(d, l, os.path.join(cwl, nName, 'traj'))
            
	with open(os.path.join(cwl, 'dbloc'), 'w') as f:
		pickle.dump(DBLoc, f)
    
	##################################################
	#Collating data
	sample, train, test = dutils.collateReports(inputfile, DBLoc)
	dutils.writeLFiles(sample, train, test)



