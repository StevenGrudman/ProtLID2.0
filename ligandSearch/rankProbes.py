import sys
import os
import operator
import subprocess
import glob
import gzip
import warnings
import shutil
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, Superimposer, PDBIO

receptor = sys.argv[1]
homeDir = os.getcwd()
protlidPath = os.getenv('PROTLIDv2HOME')
outputDir = f'{homeDir}/results'
mpDistCutoff = 4.5
rankOut = 'rankings'
os.makedirs(outputDir, exist_ok=True)
os.makedirs(f'{outputDir}/{rankOut}', exist_ok=True)

aminoAcids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

### create list of decoy ligands and cognate ligands 
ligands = []
file = f'{protlidPath}/ligandSearch/small-ig.list'
fh = open(file)
for f in fh:
	if f.strip() != receptor:
		ligands.append(f.strip())
fh.close()

file = f'{protlidPath}/ligandSearch/cognateLigands_2only/cognates.{receptor}.list'
fh = open(file)
for f in fh:
	f = f.strip()
	if f.split()[0] not in ligands:
		ligands.append(f.split()[0])
fh.close()


aminoAcids_HB_HYD = ['ALA_HYD','ARG_HB','ARG_HYD','ASN_HB','ASN_HYD','ASP_HB','ASP_HYD','CYS_HYD','GLN_HB','GLN_HYD','GLU_HB','GLU_HYD','HIS_HB','HIS_HYD','ILE_HYD','LEU_HYD','LYS_HB','LYS_HYD','MET_HYD','PHE_HYD','PRO_HYD','SER_HB','SER_HYD','THR_HB','THR_HYD','TRP_HB','TRP_HYD','TYR_HB','TYR_HYD','VAL_HYD']

### Hold cooridinates of each meshpoint in rs_MIF
mp_coords = []
### Hold preferences of each meshpoin in the rs_MIF
prefMatrix = []
file = f'{homeDir}/pharmacophore.tsv'
fh = open(file)
for f in fh:
	f = f.strip()
	f = f.split()
	hold = []
	for ff in f:
		hold.append(float(ff))
	mp_coords.append([hold[0],hold[1],hold[2]])
	prefMatrix.append(hold[3:])
fh.close()

### Get distance between two points [x1,y1,z1] and [x2,y2,z2]
def distance(x1,y1,z1,x2,y2,z2):
	d = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
	return d

### Get atom type
def HBorHYD(atom):
	atomType = 'NA'
	if atom == 'N' or atom == 'CA' or atom == 'O' or atom == 'C' or atom == 'OXT':
		atomType = 'BB'
	elif 'O' in atom or 'N' in atom:
		atomType = 'HB'
	else:
		atomType = 'HYD'
	return(atomType)

### Superimpose docked complex with the original receptor pdb the rs-MIF was created on. The rs-MIF is spatially sensitive
def superimposeDockedStructures(complx,dd):
	try:
		parser = PDBParser(QUIET=True)
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', BiopythonWarning)

			### unsip pdb and hold in 'structure'
			def parse_gzipped_pdb(gzipped_pdb_file):
				with gzip.open(gzipped_pdb_file, 'rt') as f:
					parser = PDBParser()
					structure = parser.get_structure('structure', f)
				return structure

			### Extract alpha carbon position from 'structure' 
			def get_alpha_carbons(structure,ch):
				alpha_carbons = []
				for model in structure:
					for chain in model:
						if chain.id == ch:
							for residue in chain:
								if "CA" in residue:
									alpha_carbons.append(residue["CA"])
				return alpha_carbons

			### native receptor pdb
			structure1 = parser.get_structure('structure', f'{protlidPath}/receptors/{receptor}.pdb')
			### docked complex pdb
			structure2 = parse_gzipped_pdb(f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/{dd}/{complx}.gz')

			alpha_carbons1 = get_alpha_carbons(structure1,'R')
			alpha_carbons2 = get_alpha_carbons(structure2,'R')

			super_imposer = Superimposer()
			super_imposer.set_atoms(alpha_carbons1, alpha_carbons2)
			super_imposer.apply(structure2.get_atoms())

			### Save superimposed pdb at native receptor position
			io = PDBIO()
			io.set_structure(structure2)
			with gzip.open(f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/super_{dd}/{complx}.gz', 'wt') as f_out:
				io.save(f_out)
	except:
		### In case complex is unable to be superimposed 
		return 'skip'
		pass

### Get the ProtLID score for a docked complex
def getComplexScore(complx,dd):
	file = f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/super_{dd}/{complx}.gz'
	complexCoords = []
	aminoAcids_HB_HYD_indicies = []
	take = []
	checkShift = []
	try:
		with gzip.open(file, 'rt') as fh:
			for f in fh:
				### check for bad generated PDB
				if f[0:4] == 'ATOM':
					checkShift.append(f[66:67])
				if len(set(checkShift)) > 1:
					return 0
				elif f[21:22] != 'R' and 'TER' not in f and 'END' not in f and f[0:4] == 'ATOM' and f[77:78] != 'H':
					atomType = HBorHYD(f[13:16].strip())
					complexCoords.append([float(f[30:38]),float(f[38:46]),float(f[46:54])])
					### Saves which ligand atoms to consider
					if atomType != 'BB' and f[17:20] in aminoAcids:
						take.append('yes')
						aminoAcids_HB_HYD_indicies.append(aminoAcids_HB_HYD.index(f'{f[17:20]}_{atomType}'))
					else:
						### place holders
						take.append('no')
						aminoAcids_HB_HYD_indicies.append(-1)
	except EOFError as e:
		print(f"Error: {e}. {file} cannot be opened by two different scripts at the same time")
		pass

	### for each qualifying ligand atom, add correspoding rs_MIF score to ProtLID score
	score = 0
	c = 0
	for x,y,z in complexCoords:
		cc = 0
		for xx,yy,zz in mp_coords:
			if take[c] == 'yes':
				if distance(x,y,z,xx,yy,zz) < mpDistCutoff:
					score += prefMatrix[cc][aminoAcids_HB_HYD_indicies[c]]
			cc += 1
		c += 1

	return score

skip = ['6AKQ.A','2QA9.E','2PKA.Y','1MH1.A']
dd1 = '3_emref'
dd2 = '6_caprieval'
ref = dd1.split('_')[1]
holdCombinedData = [[],[],[]]
for lig in ligands:

	### HADDOCK fails to sample the receptor 3SGQ.E with the ligands in 'skip'
	if receptor == '3SGQ.E' and lig in skip:
		continue

	### Create dictionary between rigidbody docked complex and docked complex after minimization
	rigidDict = {}
	file = f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/traceback/traceback.tsv'
	fh = open(file)
	for f in fh:
		if '00_topo' not in f and len(f) > 1:
			f = f.strip()
			f = f.split()
			rigidDict[f[4]] = f[2]
	fh.close()

	### Takes only 1 complex per cluster of docked ligand poses. Also ensures a rigidbody docked pose and minimized docked pose exists
	allQualifiedComplexes = []
	p = f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results'
	file = f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/5_clustrmsd/clustrmsd.tsv'
	fh = open(file)
	holdClust = []
	for f in fh:
		f = f.strip()
		f = f.split()
		if len(f) > 1:
			if 'rank' != f[0]:
				if int(f[3]) not in holdClust:  
					if os.path.isfile(f'{p}/1_rigidbody/{rigidDict[f[1]]}.gz') and os.path.isfile(f'{p}/3_emref/{f[1]}.gz'):
						allQualifiedComplexes.append(rigidDict[f[1]])
						allQualifiedComplexes.append(f[1])
						holdClust.append(int(f[3]))
	fh.close()


	os.makedirs(f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/super_{dd1}', exist_ok=True)
	### List of complex names and dockqscores 
	qualifiedComplexes = []
	### Holds HADDOCK model rank
	haddockScoreDict1 = {}
	### Holds HADDOCK model score
	haddockScoreDict2 = {}
	file = f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/{dd2}/capri_ss.tsv'
	fh = open(file)
	c = 1
	for f in fh:
		f = f.split('\t')
		if 'model' != f[0]:
			if f[0].split('/')[-1] in allQualifiedComplexes:
				qualifiedComplexes.append([f[0].split('/')[-1],float(f[8])])
				haddockScoreDict1[f[0].split('/')[-1]] = c
				haddockScoreDict2[f[0].split('/')[-1]] = float(f[3])
				c += 1
	fh.close()

	### Get the ProtLID scores of each complex
	poseScores = []
	for qualifiedComplex,dockq in qualifiedComplexes:
		if os.path.exists(f'{homeDir}/ligandSearch/{receptor}/{receptor}_{lig}_results/super_{dd1}/{qualifiedComplex}.gz') == False:
			hhh = superimposeDockedStructures(qualifiedComplex,dd1)
			if hhh == 'skip':
				continue
		holdCurrent = getComplexScore(qualifiedComplex,dd1)
		poseScores.append([receptor,lig,qualifiedComplex,dockq,holdCurrent,dd1])
	poseScores = sorted(poseScores, key = operator.itemgetter(4), reverse=True)

	### Save ProtLID scores for each complex
	poseFile = open(f'{outputDir}/{lig}_{ref}_scoredPoses','w')
	protLIDScoreDict1 = {}
	protLIDScoreDict2 = {}
	c = 1
	for aa,bb,cc,dd,ee,ff in poseScores:
		protLIDScoreDict1[cc] = c
		protLIDScoreDict2[cc] = ee
		poseFile.write(f'{aa}\t{bb}\t{cc}\t{dd}\t{ee}\n')
		c += 1
	poseFile.close()

	### Rank poses from ligand search by combining HADDOCK and ProtLID scores
	combinedRanks = []
	for qualifiedComplex,dockq in qualifiedComplexes:
		combinedRanks.append([qualifiedComplex,protLIDScoreDict1[qualifiedComplex] + haddockScoreDict1[qualifiedComplex],haddockScoreDict2[qualifiedComplex],protLIDScoreDict2[qualifiedComplex],dockq,dd1])
	combinedRanks = sorted(combinedRanks, key = operator.itemgetter(1), reverse=False)

	### Gets the best pose considering the top 1,10,100 poses 
	i = 0
	for num in [1,10,100]:
		hold = []
		holdDockq = -1
		cc = 0
		while cc < num:
			if combinedRanks[cc][4] > holdDockq:
				hold = [num,combinedRanks[cc][5],lig,combinedRanks[cc][0],combinedRanks[cc][4],combinedRanks[cc][1],combinedRanks[cc][2],combinedRanks[cc][3]]
				holdDockq = combinedRanks[cc][4]
			cc += 1
		holdCombinedData[i].append(hold)
		i += 1

os.makedirs(f'{outputDir}/topModels', exist_ok=True)
i = 0
### Normalize HADDOCK and ProtLID scores by highest score in each respective set
maxHadd = min([x[6] for x in holdCombinedData[i]])
maxProt = max([x[7] for x in holdCombinedData[i]])
holdCombinedData_2 = []
for aa,bb,cc,dd,ee,ff,gg,hh in holdCombinedData[i]:
	if maxHadd != 0 and maxProt != 0:
		newScore = round((gg/maxHadd + hh/maxProt)/2,3)
	elif maxHadd != 0:
		newScore = round(gg/maxHadd,3)
	elif maxProt != 0:
		newScore = round(hh/maxProt,3)
	else:
		newScore = 0
	holdCombinedData_2.append([aa,bb,cc,dd,ee,gg,hh,ff,newScore])
holdCombinedData_2 = sorted(holdCombinedData_2, key = operator.itemgetter(8), reverse=True)

### Saves the ranking of the receptors cognate ligands among the set of decoy ligands
top = 10**i
combinedFile = open(f'{outputDir}/{rankOut}/{receptor}_combined_rankings_{ref}_Top{top}.txt','w')
combinedFile.write('Top\tRefinement\tLigand\tComplex\tDockq\tHaddockScore\tProtLIDScore\tCombinedComplexRank\tCombinedScore\n')
for aa,bb,cc,dd,ee,ff,gg,hh,ii in holdCombinedData_2:
	combinedFile.write(f'{aa}\t{bb}\t{cc}\t{dd}\t{ee}\t{ff}\t{gg}\t{hh}\t{ii}\n')
	shutil.copy(f'{homeDir}/ligandSearch/{receptor}/{receptor}_{cc}_results/{dd1}/{dd}.gz', f'{outputDir}/topModels/{cc}_topModel.pdb.gz')
combinedFile.close()



