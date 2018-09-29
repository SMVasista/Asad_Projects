from __future__ import division
import os, re, sys, math, random, copy
import numpy as NP
from glob import glob
import cwgutils


def mergePDBS(elemKey, HOLDER):
	#STEP1: Identifying joining strands
	#print sorted(list(set(sorted(master)).intersection(sorted(pa[1]))))
	#print HOLDER
	if 'hybrid' not in HOLDER.keys():
		HOLDER['hybrid'] = []
		HOLDER['modified'] = []

	#R1 : Creating master template
	#print elemKey
	sEnds = list(set(sorted(elemKey[0][1])).intersection(elemKey[1][1]))[2:]

	i = 0
	j = 1

	#Applying linear transformation of j
	aaid = sEnds[0]
	print aaid

	###Identifying co-ordinate system for Ca of i
	ori = [0, 0, 0]
	for v in HOLDER[elemKey[0][0]]:
		if v[3] == aaid and v[0] == 'CA':
			sel = v
			selO = v
	ori[0] = float(sel[4])
	ori[1] = float(sel[5])
	ori[2] = float(sel[6])
	print ori
	
	###Coverting Ca of j to ori by applying translational transformation
	diff = [0, 0, 0]
	for v in HOLDER[elemKey[1][0]]:
		if v[3] == str(aaid) and v[0] == 'CA':
			sel = v
	diff[0] = float(sel[4]) - ori[0]
	diff[1] = float(sel[5]) - ori[1]
	diff[2] = float(sel[6]) - ori[2]

	print diff

	###Applying transformation of diff to j
	for elem in HOLDER[elemKey[0][0]]:
		HOLDER['hybrid'].append(elem)

	print "Finished applying linear tranformation to template"

	#Applying rotational tranformation for second aaid in sticky ends
	
	aaid = sEnds[1]
	for v in HOLDER[elemKey[0][0]]:
		if v[3] == aaid and v[0] == 'CA':
			selA = v
	for v in HOLDER[elemKey[1][0]]:
		if v[3] == aaid and v[0] == 'CA':
			selB = v

	###Identifying three vector positions using selA, selB and selO

	O = [selO[4], selO[5], selO[6]]
	A = [selA[4], selA[5], selA[6]]
	B = [selB[4], selB[5], selB[6]]

	###Identifying alpha, beta and gamma angles of both the vectors
	OA = [A[0] - O[0], A[1] - O[1], A[2] - O[2]]
	print OA
	unit_OA = float(1/(math.sqrt(OA[0]**2 + OA[1]**2 + OA[2]**2)))
	print unit_OA
	alphaA = math.acos(OA[0]*unit_OA)
	betaA = math.acos(OA[1]*unit_OA)
	gammaA = math.acos(OA[2]*unit_OA)

	OB = [B[0] - O[0] - diff[0], B[1] - O[1] - diff[1], B[2] - O[2] - diff[2]]
	unit_OB = float(1/(math.sqrt(OB[0]**2 + OB[1]**2 + OB[2]**2)))
	alphaB = math.acos(OB[0]*unit_OB)
	betaB = math.acos(OB[1]*unit_OB)
	gammaB = math.acos(OB[2]*unit_OB)

	rotAngle = [alphaA - alphaB, betaA - betaB, gammaA - gammaB]

	###Creating rotational Matrix
	Rz = NP.matrix([[math.cos(rotAngle[2]), -1*math.sin(rotAngle[2]), 0], [math.sin(rotAngle[2]), math.cos(rotAngle[2]), 0], [0, 0, 1]])
	Ry = NP.matrix([[math.cos(rotAngle[1]), 0, math.sin(rotAngle[1])], [0, 1, 0], [-1*math.sin(rotAngle[1]), 0, math.cos(rotAngle[1])]])
	Rx = NP.matrix([[1, 0, 0], [0, math.cos(rotAngle[0]), -1*math.sin(rotAngle[0])], [0, math.sin(rotAngle[0]), math.cos(rotAngle[0])]])

	R = Rz*Ry*Rx

	print R

	print "Finished rotaional tranformation for amino acid 2"

	for elem in HOLDER[elemKey[1][0]]:
		if elem[3] == sEnds[0]:
			HOLDER['modified'].append((elem[0], elem[1], elem[2], elem[3], elem[4] - diff[0], elem[5] - diff[1], elem[6] - diff[2]))
		else:
			tempMx = NP.matrix([[elem[4] - diff[0]], [elem[5] - diff[1]], [elem[6] - diff[2]]])
			rotCoOrd = NP.dot(R, tempMx)
			HOLDER['modified'].append((elem[0], elem[1], elem[2], elem[3], float(rotCoOrd[0]), float(rotCoOrd[1]), float(rotCoOrd[2])))

	#for line in HOLDER['hybrid']:
	#	print line


	#Applying second degree rotational transformation
	aaid = sEnds[1]
	for v in HOLDER[elemKey[0][0]]:
		if v[3] == aaid and v[0] == 'O':
			HselA = v
		if v[3] == aaid and v[0] == 'C':
			HselO = v
	for v in HOLDER['modified']:
		if v[3] == aaid and v[0] == 'O':
			HselB = v

	hO = [HselO[4], HselO[5], HselO[6]]
	hA = [HselA[4], HselA[5], HselA[6]]
	hB = [[HselB[4]], [HselB[5]], [HselB[6]]]
	hhB = NP.dot(R, hB)

	###Identifying alpha, beta and gamma angles of both the vectors
	hOA = [hA[0] - hO[0], hA[1] - hO[1], hA[2] - hO[2]]
	unit_hOA = float(1/(math.sqrt(hOA[0]**2 + hOA[1]**2 + hOA[2]**2)))
	halphaA = math.acos(hOA[0]*unit_hOA)
	hbetaA = math.acos(hOA[1]*unit_hOA)
	hgammaA = math.acos(hOA[2]*unit_hOA)

	hOB = [hhB[0][0] - hO[0], hhB[1][0] - hO[1], hhB[2][0] - hO[2]]
	unit_hOB = float(1/(math.sqrt(hOB[0]**2 + hOB[1]**2 + hOB[2]**2)))
	halphaB = math.acos(hOB[0]*unit_hOB)
	hbetaB = math.acos(hOB[1]*unit_hOB)
	hgammaB = math.acos(hOB[2]*unit_hOB)

	HrotAngle = [halphaA - halphaB, hbetaA - hbetaB, hgammaA - hgammaB]

	###Creating rotational Matrix
	Hz = NP.matrix([[math.cos(HrotAngle[2]), -1*math.sin(HrotAngle[2]), 0], [math.sin(HrotAngle[2]), math.cos(HrotAngle[2]), 0], [0, 0, 1]])
	Hy = NP.matrix([[math.cos(HrotAngle[1]), 0, math.sin(HrotAngle[1])], [0, 1, 0], [-1*math.sin(HrotAngle[1]), 0, math.cos(HrotAngle[1])]])
	Hx = NP.matrix([[1, 0, 0], [0, math.cos(HrotAngle[0]), -1*math.sin(HrotAngle[0])], [0, math.sin(HrotAngle[0]), math.cos(HrotAngle[0])]])

	H = Hz*Hy*Hx

	print H

	for elem in HOLDER[elemKey[1][0]]:
		if elem[3] == sEnds[0]:
			HOLDER['hybrid'].append((elem[0], elem[1], elem[2], elem[3], elem[4] - diff[0], elem[5] - diff[1], elem[6] - diff[2]))
		else:
			tempMx = NP.matrix([[elem[4] - diff[0]], [elem[5] - diff[1]], [elem[6] - diff[2]]])
			rotCoOrd = NP.dot(R, tempMx)
			hrotCoOrd = NP.dot(H, rotCoOrd)
			HOLDER['hybrid'].append((elem[0], elem[1], elem[2], elem[3], float(hrotCoOrd[0]), float(hrotCoOrd[1]), float(hrotCoOrd[2])))

	return HOLDER['hybrid']

	
def normalizeCentroids(model):
    #print model
    AAModel = []
    uID = list(set([v[3] for v in model]))
    for uid in uID:
        x = []
	y = []
	z = []
	for line in model:
	    if uid == line[3]:
		#print line
		aaName = line[1]
		x.append(line[4])
		y.append(line[5])
		z.append(line[6])
	AAModel.append([aaName, uid, NP.mean(x), NP.mean(y), NP.mean(z)])

    return AAModel
	

def extractSphere(model, aaid, radius):
	SOI = []
	for line in model:
		if str(line[1]) == aaid:
			basisCoOrd = [line[2], line[3], line[4]]
	
	for line in model:
		coOrd = [line[2], line[3], line[4]]
		if math.sqrt((coOrd[0] - basisCoOrd[0])**2 + (coOrd[1] - basisCoOrd[1])**2 + (coOrd[2] - basisCoOrd[2])**2) < radius:
			SOI.append(line)
			
	return SOI

def influenceSurface(model, aaid, cutOff):
	#'''Identifies a sphere of influence around the given amino-acid and calculates centroid based sphere co-ordinates'''
	AAModel = normalizeCentroids(model)
	#print AAModel
	r = 50*cutOff
	soi = extractSphere(AAModel, aaid, r)
	return soi
	
def identifyCatalyticPockets(data):
    '''This program identifies putative catalytic pockets from the 3D structure to calculate reactivity.'''
    seq = [x[0] for x in data]
    #Identifying dimensions
    x_min = min([x[2] for x in data])
    x_max = max([x[2] for x in data])
    y_min = min([y[3] for y in data])
    y_max = max([y[3] for y in data])
    z_min = min([z[4] for z in data])
    z_max = max([z[4] for z in data])      
    p_center = [(x_min+x_max)/2, (y_min+y_max)/2, (z_min+z_max)/2]
    udistx = (x_max - x_min)/4
    udisty = (y_max - y_min)/4
    udistz = (z_max - z_min)/4
    
    print udistx*4, udisty*4, udistz*4
        
    UNIT = {}
        
    #Identifying unit cells
    for i in range(1,4):
        for j in range(1,4):
            for k in range(1,4):
                xq = [(i-1)*udistx, i*udistx]
                yq = [(j-1)*udisty, j*udisty]
                zq = [(k-1)*udistz, k*udistz]
                    
                #Identifying all AAs lying in the box
                AA_IN_BOX = []
                for seq in data:
                    if xq[0] <= seq[2] <= xq[1] and yq[0] <= seq[3] <= yq[1] and zq[0] <= seq[4] <= zq[1]:
                        print "Identified", seq[0], "within unit", i, j ,k
                        AA_IN_BOX.append([seq[0], seq[1]])
                    
                #Calculating mean-dist within unit:
                d = 0.0
                    
                if len(AA_IN_BOX) < 2:
                    UNIT[float(0)] = [xq, yq, zq]
                else:
                    for elem1 in AA_IN_BOX:
                        for elem2 in AA_IN_BOX:
                            d += eucDist(getCoOrdinates(data, elem1[1]), getCoOrdinates(data, elem2[1]))*(1/len(AA_IN_BOX))
                    UNIT[d] = [xq, yq, zq]
    
    h = sorted(UNIT.keys(), reverse=True)
    pockets = []
    for i in h[:4]:
        pockets.append(UNIT[h[i]])
        
    return pockets
                    


########################################################################
###############Ancillary Mathematical/Geometric functions###############

def eucDist(pA, pB):
    if len(list(pA)) == 1 and len(list(pB)) == 1:
        return abs(float(pB) - float(pA))

    elif len(list(pA)) == 2 and len(list(pB)) == 2:
        return math.sqrt((pB[0] - pA[0])**2 + (pB[1] - pA[1])**2)
        
    elif len(list(pA)) == 3 and len(list(pB)) == 3:
        return math.sqrt((pB[0] - pA[0])**2 + (pB[1] - pA[1])**2 + (pB[2] - pA[2])**2)

def getCoOrdinates(coOrdData, pos):
    for element in coOrdData:
        if element[1] == pos:
            return (element[2], element[3], element[4])










