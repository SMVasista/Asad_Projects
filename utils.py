from __future__ import division
import os, re, sys, math, random, copy, pickle
import numpy as NP
from glob import glob
import cwgutils
import pdbutils as pdb
import datutils as dutils
import bioprop as bp

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
                mutations.append(dutils.parseMutation(line.split(':')[1]))
    return projName, pdbLoc, mutations

def extractFromSimulation(inputfile, DBLoc):
    projName, pdbloc, mutList = parseInput(inputfile)
    NSFL = dutils.parseAllCombinations([2000], projName, mutList, inputfile)
    with open(DBLoc, 'r') as f:
    	location = pickle.load(f)
    sample, train, test = dutils.collateReports(inputfile, location)
    dutils.writeLFiles(sample, train, test, NSFL)

if __name__=="__main__":
    script, inputfile, arg = sys.argv
    if arg == 'source' or arg == 'Source' or arg == "SOURCE":
        dbLoc = raw_input('DB-Loc file: ')
        extractFromSimulation(inputfile, dbLoc)

    #Add other functions for processing data
