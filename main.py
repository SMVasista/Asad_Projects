from __future__ import division
import os, re, sys, math, random, copy
import numpy as NP
from glob import glob
import cwgutils
import pdbutils as pdb
import geometrix as gx
import energetix as ex

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


if __name__=="__main__":
	script, fileLoc, projName = sys.argv
	H = pdb.readPDB(fileLoc)
	uModel = pdb.stringPDB(H)
	#mutation = raw_input('Enter mutation : ')
	#mutation = parseMutation(mutation)
	aaCordData = gx.influenceSurface(uModel, '10', cutOff=1)
	print aaCordData
	
	##################################################
	#Applying Energy Minimization
	#C, E = ex.minimizeEnergyStruc(aaCordData)
	#pockets = gx.identifyCatalyticPockets(C[-1])
	#The pockets identified here are putative.
	#ex.calculateEnergyFunc(pockets, C, E)
    
    ##################################################
    #Applying mutation 
    #mutList = queryMutation(mutation)
    #seq = generateProteinSequence(aaCordData, mutList)
    
    ##################################################
    #Applying Energy Minimization (on mutations)
    #for mut in mutList:
    #    C, E = ex.minimizeEnergyStruc(seq[mut])
    #    ex.calculateEnergyFunc(seq[mut])
    
    ##################################################
    #Collating data
    
    
    



