##############################################
######## AminoAcid feature ###################
########### ARKO || 12-july-2018 #############
##############################################

############################Using __main__.__dict__.items() method to call biochemical properties

Hydrophobicity = {'A':0.87,'L':2.17,'R':0.85,
'K':1.64,'N':0.09,'M':1.67,
'D':0.66,'F':2.87,'C':1.52,
'P':2.77,'Q':0,'S':0.07,'E':0.67,
'T':0.07,'G':0.1,'W':3.77,'H':0.87,
'Y':2.67,'I':3.15,'V':1.87}
		
pKa = {'A':4.76,'L':4.79,'R':4.3,'K':4.27,
'N':3.64,'M':4.25,'D':5.69,'F':4.31,
'C':3.67,'P':0,'Q':4.54,'S':3.83,
'E':5.48,'T':3.87,'G':3.77,'W':4.75,
'H':2.84,'Y':4.3,'I':4.81,'V':4.86}

Flexibility = {'A':1.315,'L':1.234,'R':1.31,'K':1.367,'N':1.38,
'M':1.269,'D':1.372,'F':1.247,'C':1.196,'P':1.342,
'Q':1.342,'S':1.381,'E':1.376,'T':1.324,'G':1.382,
'W':1.186,'H':1.279,'Y':1.199,'I':1.241,'V':1.233}

Charge = {'A':0,'L':0,'R':1,'K':1,'N':0,
'M':0,'D':-1,'F':0,'C':0,'P':0,
'Q':0,'S':0,'E':-1,'T':0,'G':0,
'W':0,'H':1,'Y':0,'I':0,'V':0}

Radii_gyration = {'A':0.77,'L':1.54,'R':2.38,'K':2.08,'N':1.45,
'M':1.8,'D':1.43,'F':1.9,'C':1.22,'P':1.25,
'Q':1.75,'S':1.08,'E':1.77,'T':1.24,'G':0.58,
'W':2.21,'H':1.78,'Y':2.13,'I':1.56,'V':1.29}

HBond_donor_prop = {'R': 2.5, 'N': 1, 'D': 0, 'Q': 0.5, 'E': 0, 'H': 1, 'K': 1.5, 'S': 1, 'T': 1, 'W': 0.5, 'Y': -1, 'A': 0, 'L': 0, 'M': 0, 'F': 0, 'C': 0, 'P': 0, 'G': 0, 'I': 0, 'V': 0}

HBond_acceptor_prop = {'R': 0, 'N': -2, 'D': -4, 'Q': -2, 'E': -4, 'H': -1, 'K': 0, 'S': -2, 'T': -2, 'W': 0, 'Y': -1, 'A': 0, 'L': 0, 'M': 0, 'F': 0, 'C': 0, 'P': 0, 'G': 0, 'I': 0, 'V': 0}

Salt_Bridge_prop = {'L': {'Q': 1, 'N': 1}, 'R' : {'Q': 1, 'N': 1}}

Disf_bond_prop = {'C': {'C': 1}}

Free_energy_transfer = {'A':0.3,'L':0.5,'R':-1.4,'K':-1.8,'N':-0.5,
'M':0.4,'D':-0.6,'F':0.5,'C':0.9,'P':-0.3,
'Q':-0.7,'S':-0.1,'E':-0.7,'T':-0.2,'G':0.3,
'W':0.3,'H':-0.1,'Y':-0.4,'I':0.7,'V':0.6}

VDW_vol = {'A':1,'L':4,'R':6.13,'K':4.77,'N':2.95,
'M':4.43,'D':2.78,'F':5.89,'C':2.43,'P':2.72,
'Q':3.95,'S':1.6,'E':3.78,'T':2.6,'G':1,
'W':8.08,'H':4.66,'Y':6.47,'I':4,'V':3}

def map_aa_names(AminoAcid, key):
    aaMap = {'G': ('Glycine','GLY'), 'P': ('Proline','PRO'), 'A': ('Alanine','ALA'), 'V': ('Valine','VAL'), 'L': ('Leucine','LEU'), 'I': ('Isoleucine','ILE'), 'M': ('Methionine','MET'), 'C': ('Cysteine','CYS'), 'F': ('Phenylalanine','PHE'), 'Y': ('Tyrosine','TYR'), 'W': ('Tryptophan','TRP'), 'H': ('Histidine','HIS'), 'K': ('Lysine','LYS'), 'R': ('Arginine','ARG'), 'Q': ('Glutamine','GLN'), 'N': ('Asparagine','ASN'), 'E': ('Glutamic Acid','GLU'), 'D': ('Aspartic Acid','ASP'), 'S': ('Serine','SER'), 'T': ('Threonine','THR')}
    if len(list(AminoAcid)) == 1 and key == '3':
        for k in aaMap.keys():
            if AminoAcid == k:
                return aaMap[k][1]
    elif len(list(AminoAcid)) == 3 and key == '1':
        for k, v in aaMap.items():
            if AminoAcid == v[1]:
                return k
    elif len(list(AminoAcid)) == 1 and key == '1':
        return AminoAcid
    elif len(list(AminoAcid)) == 3 and key == '3':
        return AminoAcid
    #else:
    #    print "Amino acid map not found"
    #    return None
                
if __name__ == '__main__':
	pass




