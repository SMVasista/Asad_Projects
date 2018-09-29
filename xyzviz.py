from __future__ import division
import os, re, sys, math, random, copy, pickle
from glob import glob
import cwgutils

def extractSimulData(fileLoc):
    with open(fileLoc, 'r') as f:
        data = pickle.load(f)
    
    #Capature data-labels
    LABELS = []
    for telem in data[0]:
        LABELS.append(str(telem[1])+'_'+str(telem[0]))
        
    LABELS = list(set(LABELS))
            
    return data, LABELS
    
def writeTrajectory(data, LABELS, outputLoc):

    print "Writing trajectories"
    
    for aa in LABELS:
        with open(os.path.join(outputLoc, str(aa)), 'a') as f:
            for elem in data:
                for key in elem:
                    if key[1] == int(aa.split('_')[0]):
                        f.write(str(key[0])+' '+str(key[1])+' '+str(key[2])+' '+str(key[3])+' '+str(key[4])+'\n')
                        
    print "Finished writing individual trajectories"
    
    for i, elem in enumerate(data):
        with open(os.path.join(outputLoc, str(i)), 'a') as f:
            for key in elem:
                f.write(str(key[0])+' '+str(key[1])+' '+str(key[2])+' '+str(key[3])+' '+str(key[4])+'\n')
    print "Done"

if __name__=="__main__":
    script, fileLoc, outputLoc = sys.argv
    data, labels = extractSimulData(fileLoc)
    writeTrajectory(data, labels, outputLoc)
    

                
