from __future__ import division
import os, re, sys, math, random, copy, pickle
from glob import glob
import cwgutils

cwl = os.getcwd()

Features = ['rmsf', 'emin', 'HBondTotal', 'HBondEnergy', 'pkt1_dvol', 'pkt1_chr', 'pkt1_hydp', 'pkt1_ienergy', 'pkt2_dvol', 'pkt2_chr', 'pkt2_hydp', 'pkt2_ienergy', 'pkt3_dvol', 'pkt3_chr', 'pkt3_hydp', 'pkt3_ienergy', 'pkt4_dvol', 'pkt4_chr', 'pkt4_hydp', 'pkt4_ienergy', 'outcome']

def collateReports(mutData, DBLoc):
	'''Collecting data'''
	sample = {}
        for loc in DBLoc:
            mName = os.path.basename(loc)
            sample[mName] = {'outcome': None, 'rmsf': 0.0, 'emin': 0.0, 'HBondTotal': 0.0, 'HBondEnergy': 0.0, 'pkt1_dvol': 0.0, 'pkt1_chr': 0.0, 'pkt1_hydp': 0.0, 'pkt1_ienergy': 0.0, 'pkt2_dvol': 0.0, 'pkt2_chr': 0.0, 'pkt2_hydp': 0.0, 'pkt2_ienergy': 0.0, 'pkt3_dvol': 0.0, 'pkt3_chr': 0.0, 'pkt3_hydp': 0.0, 'pkt3_ienergy': 0.0, 'pkt4_dvol': 0.0, 'pkt4_chr': 0.0, 'pkt4_hydp': 0.0, 'pkt4_ienergy': 0.0}
            data = cwgutils.readLinesAndSplit(os.path.join(loc, 'stability_parameters.csv'), ' ')
            sample[mName]['rmsf'] = float(data[-1][1])
            sample[mName]['emin'] = float(data[-1][2])
            sample[mName]['HBondTotal'] = float(data[-1][3])
            sample[mName]['HBondEnergy'] = float(data[-1][4])
        
	    data2 = cwgutils.readLines(os.path.join(loc, 'catalytic_pocket_parameters.csv'))
            volume = [x for x in data2 if 'Volume Encroachment' in x]
            if len(volume) == 4:
                sample[mName]['pkt1_dvol'] = float(volume[0].split(':')[1])
                sample[mName]['pkt2_dvol'] = float(volume[1].split(':')[1])
                sample[mName]['pkt3_dvol'] = float(volume[2].split(':')[1])
                sample[mName]['pkt4_dvol'] = float(volume[3].split(':')[1])
            charge = [x for x in data2 if 'Charge Density' in x]
            if len(charge) == 4:
                sample[mName]['pkt1_chr'] = float(charge[0].split(':')[1])
                sample[mName]['pkt2_chr'] = float(charge[1].split(':')[1])
                sample[mName]['pkt3_chr'] = float(charge[2].split(':')[1])
                sample[mName]['pkt4_chr'] = float(charge[3].split(':')[1])
            hphi = [x for x in data2 if 'Pocket Hydrophobicity' in x]
            if len(hphi) == 4:
                sample[mName]['pkt1_hydp'] = float(hphi[0].split(':')[1])
                sample[mName]['pkt2_hydp'] = float(hphi[1].split(':')[1])
                sample[mName]['pkt3_hydp'] = float(hphi[2].split(':')[1])
                sample[mName]['pkt4_hydp'] = float(hphi[3].split(':')[1])
            energy = [x for x in data2 if 'Internal Energy' in x]
            if len(energy) == 4:
                sample[mName]['pkt1_ienergy'] = float(energy[0].split(':')[1])
                sample[mName]['pkt2_ienergy'] = float(energy[1].split(':')[1])
                sample[mName]['pkt3_ienergy'] = float(energy[2].split(':')[1])
                sample[mName]['pkt4_ienergy'] = float(energy[3].split(':')[1])
            print "Finished collecting mutation MDS data"
        
	#Reading mutation outcomes data
	trainSet = []
	testSet = []
	mdata = cwgutils.readLines(mutData)
	for line in mdata:
	    if 'gene_name' in line:
	        basisEntry = line.split(':')[1]
	        sample[basisEntry]['outcome'] = 0.0
	        trainSet.append(basisEntry)
	for line in mdata:
	    if 'mutation' in line:
	        mName = str(basisEntry)+'_'+str(line.split(':')[1])
		print line.split(':')
            	if line.split(':')[2] == 'GOF' or line.split(':')[2] == 'SOF':
                	sample[mName]['outcome'] = 1.0
                	trainSet.append(mName)
            	elif line.split(':')[2] == 'LOF':
                	sample[mName]['outcome'] = -1.0
                	trainSet.append(mName)
            	elif line.split(':')[2] == 'COF':
                	sample[mName]['outcome'] = 0.0
                	trainSet.append(mName)
            	else:
                	testSet.append(mName)
                
	return sample, trainSet, testSet

def writeLFiles(sample, trainSet, testSet):
    '''This program writes collected data into training file and test file for the ANN'''
    with open(os.path.join(cwl, 'ann_input_train.csv'), 'a') as f:
        f.write('Blk'+',')
        for entry in sorted(trainSet):
            f.write(str(entry)+',')
        f.write('\n')
        for feature in Features:
            f.write(str(feature)+',')
            for entry in sorted(trainSet):
                f.write(str(sample[entry][feature])+',')
            f.write('\n')

    with open(os.path.join(cwl, 'ann_input_test.csv'), 'a') as f:
        f.write('Blk'+',')
        for entry in sorted(testSet):
            f.write(str(entry)+',')
        f.write('\n')
        for feature in Features:
            if feature != 'outcome':
                f.write(str(feature)+',')
                for entry in sorted(testSet):
                    f.write(str(sample[entry][feature])+',')
                f.write('\n')  
                
                   
    
        
        
        
        
        
        
        
        
        
        
        
        
