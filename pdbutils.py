from __future__ import division
import os, re, sys, math, random, copy
import numpy as NP
from glob import glob
import geometrix as gx
import cwgutils

def readPDB(fileLoc):
	'''typically reads a location which contain one or more files with .pdb extention. Assumes all pdb for same protein'''
	HOLDER = {}
	files = glob(os.path.join(str(fileLoc),'*.pdb'))
	for fl in files:
		HOLDER[fl] = []
		#data = cwgutils.readLinesAndSplit(fl, ' ')
		#print data
		with open(fl, 'r') as fo:
			for line in fo:
				'''ATOM, AACID, CHAIN,AAPOS,X,Y,Z'''
				if line.startswith('ATOM'):
					line = list(line)
					line = [' ']+line
					ATOM = str(line[13:16]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', '')
					AAcid = str(line[18:21]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', '')
					Chain = line[22]
					AApos = int(str(line[23:27]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', ''))	
					X = float(str(line[31:39]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', ''))
					Y = float(str(line[39:47]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', ''))
					Z = float(str(line[47:55]).replace(',', '').replace('[', '').replace(']', '').replace("'", '').replace(' ', ''))
					#print ATOM,AAcid, Chain, AApos, X, Y, Z			
					HOLDER[fl].append((ATOM,AAcid, Chain, AApos, X, Y, Z))
	
	return HOLDER

def stringPDB(HOLDER):
	FRAG = {}
	if len(HOLDER.keys()) > 1:
	    for i, elem in enumerate(HOLDER):
		    #print HOLDER[elem]
		    cov = sorted(list(set([int(k[3]) for k in HOLDER[elem]])))
		    #print '###cov%s'%cov
		    FRAG[len(cov)] = (cov, elem)
	    #print FRAG

	    newFRAG = []
	    elemKey = []

	    for i in sorted(FRAG.keys(), reverse=True):
		    if len(sorted(list(set(newFRAG).union(FRAG[i][0])))) > len(sorted(newFRAG))+2:
			    newFRAG = newFRAG + FRAG[i][0]
			    elemKey.append((FRAG[i][1], FRAG[i][0]))
			    print "adding peptide", FRAG[i][1], " of length:", i
			
		    else:
			    print "ignoring peptide", FRAG[i][1], " of length:", i
	    #print elemKey
	    model = gx.mergePDBS(elemKey, HOLDER)
	    return model
	    
	else:
	    return HOLDER[HOLDER.keys()[0]]

if __name__ == '__main__':
	pass


