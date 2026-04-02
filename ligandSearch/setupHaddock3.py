import sys
import os
import subprocess
import shutil
import time

homeDir = os.getcwd()
homeDirName = homeDir.split('/')[-1]
protlidPath = os.getenv('PROTLIDv2HOME')
outputDir = f'{homeDir}/ligandSearch'

receptor = sys.argv[1]
### Job name cannot exceed 8 characters becasue jobId from job name needs to be extracted to make later dependency 
if len(receptor) != 6:
	sys.exit(f'EXIT: receptor name must be 6 chareacters long. Example: 1AK4.A')
pdb = receptor.split('.')[0]
nativeligand = sys.argv[2]
qc = receptor.split('.')[1]
ic = nativeligand.split('.')[1]

os.makedirs(f'{homeDir}/ligandSearch',exist_ok=True)
os.makedirs(outputDir, exist_ok=True)
os.makedirs(f'{outputDir}/err', exist_ok=True)
os.makedirs(f'{outputDir}/output', exist_ok=True)
os.makedirs(f'{outputDir}/haddock3_jobs', exist_ok=True)

### Haddock will not run if the output directory already exists.
if os.path.isdir(f'{outputDir}/{receptor}'):
	sys.exit(f'EXIT: {outputDir}/{receptor} already exists. Please delete directory to continue')
os.makedirs(f'{outputDir}/{receptor}')
os.makedirs(f'{outputDir}/{receptor}/data')
shutil.copy(f'{protlidPath}/ligands/{nativeligand}.pdb', f'{outputDir}/{receptor}/data')
shutil.copy(f'{protlidPath}/receptors/{receptor}.pdb', f'{outputDir}/{receptor}/data')

### Create new pdb complex replacing chain names with R for receptor and L for ligand
newPdbFile = open(f'{outputDir}/{receptor}/data/{pdb}.pdb','x')
file = f'{protlidPath}/pdbs/{pdb}.pdb'
fh = open(file)
for f in fh:
	if (f[0:4] == 'ATOM' or f[0:3] == 'TER') and f[21:22] == qc:
		newPdbFile.write(f'{f[0:21]}R{f[22:]}')
	if (f[0:4] == 'ATOM' or f[0:3] == 'TER') and f[21:22] == ic:
		newPdbFile.write(f'{f[0:21]}L{f[22:]}')
fh.close()
newPdbFile.close()

### Get interface between native receptor and ligand
intercaat = subprocess.check_output(['python3', f'{protlidPath}/protlid_auxPrograms/intercaat/intercaat.py', '-pdb', f'{pdb}.pdb', '-qc', qc,\
'-ic', ic, '-fp', f'{protlidPath}/pdbs/'], shell =False, text= True)
intercaat = intercaat.split('\n')
interface = []
for f in intercaat:
	f = f.strip()
	f = f.split()
	if len(f) == 14:
		if f[1] not in interface:
			interface.append(f[1])

### populate a list 'ligands' including all decoy ligands and cognate ligands
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


### create a SLURM submit script to run HADDOCK on the receptor and all the ligands in the list 'ligands'
numOfJobs = len(ligands)
submitFile = open(f'{outputDir}/submit_{receptor}_1.sh','w')
submitFile.write(f'''#!/bin/bash
#SBATCH --job-name={receptor}-1
#SBATCH --output={outputDir}/output/LOGFILE_{receptor}_1.out
#SBATCH --error={outputDir}/err/LOGFILE_{receptor}_1.err
#SBATCH --array=1-{numOfJobs}

cd $SLURM_SUBMIT_DIR

CMDFILE={outputDir}/haddock3_jobs/jobs_{receptor}_1
CMD=$(awk "NR==$SLURM_ARRAY_TASK_ID" $CMDFILE)
source ~/.bashrc
conda activate haddock3
$CMD''')
submitFile.close()


jobFile = open(f'{outputDir}/haddock3_jobs/jobs_{receptor}_1','w')
for ligand in ligands:

	jobFile.write(f'haddock3 {outputDir}/{receptor}/docking_{receptor}_{ligand}.tbl\n')
	shutil.copy(f'{protlidPath}/ligands/{ligand}.pdb', f'{outputDir}/{receptor}/data')
	### run naccess on ligand PDB to get ligand surface accessiblity 
	os.system(f'bash {protlidPath}/ligandSearch/get_SA.PDB.sh {protlidPath}/ligands/{ligand}.pdb {protlidPath}/SA_PDBs {homeDirName}')

	### create list of all surface assessible ligand residues
	sa_ligand_res = []
	file = f'{protlidPath}/SA_PDBs/{ligand}.SA'
	fh = open(file)
	for f in fh:
		if f[22:27].strip()not in sa_ligand_res:
			sa_ligand_res.append(f[22:27].strip())
	fh.close()

	### Create restraints file for site directed docking at the receptor interface
	qc = 'R'
	ic = 'L'
	newFile = open(f'{outputDir}/{receptor}/data/{receptor}_{ligand}_restraint.tbl','w')
	for i in interface:
		newFile.write('!\n')
		newFile.write(f'assign ( resid {i} and segid {qc})\n')
		newFile.write( '       (\n')
		for r in sa_ligand_res:
			newFile.write(f'         ( resid {r} and segid {ic})\n')
			if r != sa_ligand_res[-1]:
				newFile.write( '     or\n')
		newFile.write('       )  2.0 0.0 2.0\n')
	newFile.close()

	### create run file for HADDOCK.
	newFile = open(f'{outputDir}/{receptor}/docking_{receptor}_{ligand}.tbl','w')
	newFile.write(f'''# ====================================================================
# directory in which the scoring will be done
run_dir = "{outputDir}/{receptor}/{receptor}_{ligand}_results"

# execution mode
mode = "local"
ncores = 1

# molecules to be docked
molecules =  [
    "{outputDir}/{receptor}/data/{receptor}.pdb",
    "{outputDir}/{receptor}/data/{ligand}.pdb"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]

[rigidbody]
ambig_fname = "{outputDir}/{receptor}/data/{receptor}_{ligand}_restraint.tbl"
sampling = 1000

[caprieval]
reference_fname = '{outputDir}/{receptor}/data/{pdb}.pdb'
receptor_chain = 'R'
ligand_chains = ['L']

[emref]
ambig_fname = "{outputDir}/{receptor}/data/{receptor}_{ligand}_restraint.tbl"
dielec = 'rdie'
nemsteps = 300

[ilrmsdmatrix]
receptor_chain = 'R'
ligand_chains = ['L']

[clustrmsd]
n_clusters = 1000
min_population = 1
clust_cutoff = 2

[caprieval]
reference_fname = '{outputDir}/{receptor}/data/{pdb}.pdb'
receptor_chain = 'R'
ligand_chains = ['L']


# ====================================================================''')
	newFile.close()
jobFile.close()

leftBrack = '{'
rightBrack = '}'
os.system(f'sbatch {outputDir}/submit_{receptor}_1.sh')


##################

### After HADDOCK is run for all receptor-ligand complexes, submit "rankProbes.py" script to rank ligands
leftBrack = '{'
rightBrack = '}'
jobid = subprocess.check_output([f'squeue --me | grep {receptor}-1 | awk \'{leftBrack}split($1,a,\"_\"); print a[1]{rightBrack}\' | sort -r | head -1'],shell=True, text=True).strip()

submitFile = open(f'{outputDir}/submit_{receptor}_ranking.sh','w')
submitFile.write(f'''#!/bin/bash
#SBATCH --job-name={receptor}-r
#SBATCH --output={outputDir}/output/LOGFILE_{receptor}_r.out
#SBATCH --error={outputDir}/err/LOGFILE_{receptor}_r.err

cd $SLURM_SUBMIT_DIR

source ~/.bashrc
conda activate haddock3
python3 {protlidPath}/ligandSearch/rankProbes.py {receptor}
$CMD''')
submitFile.close()

os.system(f'sbatch --dependency=afterok:{jobid} {outputDir}/submit_{receptor}_ranking.sh')

