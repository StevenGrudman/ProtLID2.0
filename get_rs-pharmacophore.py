import numpy as np
import glob
import matplotlib.pyplot as plt
import statistics
import operator
import os

homeDir = os.getcwd()
protlidPath = os.getenv('PROTLIDv1HOME')
aminoAcids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

### Get distance between two points [x1,y1,z1] and [x2,y2,z2]
def distance(x1,y1,z1,x2,y2,z2):
	d = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
	return d

### Holds list of atoms for each amino acid
aaAtomDict = {}
### Holds list of atom types [HB,HYD,BB] for each amino acid
aaTypeDict = {}
hold = ['ALA']
hold2 = []
hold3 = []
file = f'{protlidPath}/sharedFiles/aaAtomsType.txt'
fh = open(file)
for f in fh:
	f = f.strip()
	f = f.split()
	if f[0] not in hold:
		aaAtomDict[currentAtom] = hold2
		hold2 = []
		aaTypeDict[currentAtom] = hold3
		hold3 = []
		hold.append(f[0])
	hold2.append(f[1])
	hold3.append(f[2])
	currentAtom = f[0]
fh.close()
aaAtomDict[currentAtom] = hold2
aaTypeDict[currentAtom] = hold3

### Consider frames 25 to 125 from the MD simulations. (5 ps to 25 ps)
start = 25
end = 125

### Get all interface residue numbers in contact with mesh
interfaceAtomNumbers = []
file = f'{homeDir}/RequiredFiles/mesh.ref.xyz_'
fh = open(file)
for f in fh:
	f = f.strip()
	f = f.split()
	if int(f[5]) not in interfaceAtomNumbers:
		interfaceAtomNumbers.append(int(f[5]))
fh.close()

### Get receptor Coordinates
rec_x = []
rec_y = []
rec_z = [] 
file = f'{homeDir}/RequiredFiles/rec.reformat.pdb'
fh = open(file)
for f in fh:
	f = f.strip()
	ff = f.split()
	if ff[0] == 'ATOM':
		if int(ff[1]) in interfaceAtomNumbers:
			rec_x.append(float(f[30:38]))
			rec_y.append(float(f[38:46]))
			rec_z.append(float(f[46:54]))
fh.close()

### append all converged amino acids to the list 'results'
results = []
for aa in aminoAcids:
	tightness = [[] for x in aaAtomDict[aa]]
	closestDist_array = [[] for x in aaAtomDict[aa]]
	center_array = [[] for x in aaAtomDict[aa]]
	meshPoints = []

	### Get a list of all MD simulations for a particular amino acid,
	files = glob.glob(f'{homeDir}/PDB_07/{aa}/*')
	for file in files:
		mp = file.split('/')[-1]
		mp = mp.split('.')
		meshPoints.append(f'{mp[2]}.{mp[3]}')

		### Get the coordinates of each atom of the single residue probe over the entire MD simulation 
		coords = [[] for x in aaAtomDict[aa]]
		fh = open(file)
		for f in fh:
			f = f.strip()
			ff = f.split()
			if 'ATOM' in f:
				if f[77:78] != 'H':
					coords[aaAtomDict[aa].index(ff[2])].append([float(f[30:38]),float(f[38:46]),float(f[46:54])])
		fh.close()

		### For each atom in the probe
		for c,ccc in enumerate(coords):
			### Get the average position over the simulation
			center = np.array(ccc[start:end]).mean(0)
			center_array[c].append(list(center))

			### Get the smallest distance between the average position of the atom and the receptor over the simulation
			hold = 1000
			cc = 0
			while cc < len(rec_x):
				dist = distance(center[0],center[1],center[2],rec_x[cc],rec_y[cc],rec_z[cc])
				if dist < hold:
					hold = dist
				cc += 1
			closestDist_array[c].append(hold)

			### Calculate average atom stability. In other words, save the average distance from the atom position in each frame to the average postion 
			hold_cen = []
			for x,y,z in ccc[start:end]:
				hold_cen.append(distance(x,y,z,center[0],center[1],center[2]))
			tightness[c].append(np.mean(np.array(hold_cen)))

	### Save all converged probes for a particular amino acid to the list 'results'
	cc = 0
	while cc < len(meshPoints):
		### find which atom in the probe is the most stable and closest to the receptor
		hold = []
		c = 0
		while c < len(aaAtomDict[aa]):
			hold.append(closestDist_array[c][cc]+tightness[c][cc])
			c += 1
		### See if individual probe converges based on the the criteria in the below if statement
		i = hold.index(min(hold))
		if tightness[i][cc] < 3.5 and closestDist_array[i][cc] < 5.5:
			if aaTypeDict[aa][i] != 'BB':
				results.append([aa,meshPoints[cc],aaTypeDict[aa][i],center_array[i][cc][0],center_array[i][cc][1],center_array[i][cc][2]])
			# print(aa,meshPoints[cc],aaAtomDict[aa][i],aaTypeDict[aa][i])
			# print(tightness[i][cc],closestDist_array[i][cc],aaTypeDict[aa][i])
		cc += 1


### Create dictionaries between meshpoint indices, coordinates, and closestReceptorAtom
mp_indices = []
file = glob.glob(f'{homeDir}/RequiredFiles/meshindex.*')[0]
fh = open(file)
for f in fh:
	mp_indices.append(int(f))
fh.close()

mpCoordDict = {}
mpIndexDict = {}
mp_coords = []
c = 0
file = f'{homeDir}/RequiredFiles/mesh.xyz_'
fh = open(file)
for f in fh:
	f = f.strip()
	f = f.split()
	mp_coords.append([float(f[0]),float(f[1]),float(f[2])])
	mpCoordDict[f'{float(f[0])}_{float(f[1])}_{float(f[2])}'] = mp_indices[c]
	mpIndexDict[mp_indices[c]] = f'{f[7]}.{f[8]}.{f[6]}'
	c += 1
fh.close()

### All Meshpoints within the proximity of a converged probe are assigned a fraction of the preference. For example if 3 meshpoints are within proximity, the are all assigned 1/3 of the total preference
probe_fractional_value = []
for aa,mp,HB_HYD,x,y,z in results:
	if HB_HYD == 'HB':
		distCutoff = 1.5
	elif HB_HYD == 'HYD':
		distCutoff = 3
	hold = 0
	for mp_x,mp_y,mp_z in mp_coords:
		if distance(x,y,z,mp_x,mp_y,mp_z) < distCutoff:
			hold += 1
	probe_fractional_value.append(hold)

### Fill rs-pharmacophore with preferences
hold_mp = []
### fractional sum of convereged probes at each meshpoint 
aminoAcidCounts = []
aminoAcids_HB_HYD = ['ALA_HYD','ARG_HB','ARG_HYD','ASN_HB','ASN_HYD','ASP_HB','ASP_HYD','CYS_HYD','GLN_HB','GLN_HYD','GLU_HB','GLU_HYD','HIS_HB','HIS_HYD','ILE_HYD','LEU_HYD','LYS_HB','LYS_HYD','MET_HYD','PHE_HYD','PRO_HYD','SER_HB','SER_HYD','THR_HB','THR_HYD','TRP_HB','TRP_HYD','TYR_HB','TYR_HYD','VAL_HYD']
c = 0
for aa,mp,HB_HYD,x,y,z in results:
	if HB_HYD == 'HB':
		distCutoff = 1.5
	elif HB_HYD == 'HYD':
		distCutoff = 3
	for mp_x,mp_y,mp_z in mp_coords:
		if distance(x,y,z,mp_x,mp_y,mp_z) < distCutoff:
			if mpCoordDict[f'{mp_x}_{mp_y}_{mp_z}'] not in hold_mp:
				hold_mp.append(mpCoordDict[f'{mp_x}_{mp_y}_{mp_z}'])
				aminoAcidCounts.append([0]*30)
			i = hold_mp.index(mpCoordDict[f'{mp_x}_{mp_y}_{mp_z}'])
			j = aminoAcids_HB_HYD.index(f'{aa}_{HB_HYD}')
			aminoAcidCounts[i][j] += 1/probe_fractional_value[c]
	c += 1


### total converged amino acids
totalAminoAcidCounts = np.array(aminoAcidCounts).sum(axis=0)
### proportion each individual amino acid has to sum of all convereged amino acids 
percentOfTotal = np.array([100*x/sum(totalAminoAcidCounts) for x in totalAminoAcidCounts])

### Probes with more surviving probes are given greater weights. Accounts for amino acids with low statistical power
rrr = []
for a in aminoAcidCounts:
	### Weight is a hyperbolic function 
	weight = 1-1/(0.5*percentOfTotal+1)
	rrr.append(weight*100*np.array(a)/np.array(totalAminoAcidCounts))
rrr = np.array(rrr)

### Save rs-pharmacophore. The first 3 columns are the meshpoint coordinates and the rest of the columns are the preferences in the 'aminoAcids_HB_HYD' list
newFile = open(f'{homeDir}/pharmacophore.tsv','w')
c = 0
while c < len(hold_mp):
	i = mp_indices.index(hold_mp[c])
	s = ''
	for x in rrr[c]:
		s += f'{x}\t'
	s = s[0:-1]
	newFile.write(f'{mp_coords[i][0]}\t{mp_coords[i][1]}\t{mp_coords[i][2]}\t{s}\n')
	c += 1
newFile.close()

