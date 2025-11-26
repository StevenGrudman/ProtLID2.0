NAME
ProtLID2.0 (Protein-Ligand Interface Design 2.0)

SYNOPSIS
ProtLID2.0 operates in 2 steps: rs-pharmacophore generation followed by a ligand search.
ProtLID2.0 has been shown to be a reliable scoring function for site-directed protein–protein docking and as a tool for identifying a receptor's cognate ligand from a sub-proteome of decoys.


REQUIREMENTS:
perl (Tested with version 5.16.3)
modeller (Tested with version 10.3)
amber (Tested with amber20)
python3 (Tested with version 3.9.12)
haddock3 (Tested with version 2024.12.0b7)
	Must create a conda env for haddock3. Name this env "haddock3" 

INSTRUCTIONS
Most commands in this readme can be copied and pasted to run ProtLID2.0. However anything between [] must be updated. For example, [pdb complex] must be replaed with a pdb complex like 1I8L.pdb.

### ADD FOLLOWING CODE TO .bashrc. Make sure to update PROTLIDv1HOME path and amber path if necessary 
#######################################
export PROTLIDv1HOME=[path to protlidv2.0]
export LD_LIBRARY_PATH=$PROTLIDv1HOME/protlid_auxPrograms/naccess2.1.1:$LD_LIBRARY_PATH
PATH=$PATH:$PROTLIDv1HOME/scripts:$PROTLIDv1HOME/code_interfaceDesign:$PROTLIDv1HOME/code_partnerSearch:$PROTLIDv1HOME/code_interfaceDefinition:$PROTLIDv1HOME/bin:$PROTLIDv1HOME/protlid_auxPrograms:$PROTLIDv1HOME/protlid_auxPrograms/naccess2.1.1
export PATH=$PATH

test -f  /usr/local/bio/amber20/amber.sh && source  /usr/local/bio/amber20/amber.sh
#######################################


### GENERATE RS-PHARMACOPHORE (!!!submits scripts to SGE!!!)
#######################################
### Create a new directory and copy a pdb complex into this directory. !All following commands should be executed from this directory!

### get interface of receptor and define mesh points
$PROTLIDv1HOME/code_interfaceDefinition/run_genIntAtomList_INTERCAATd4noIllHydRes.sh [pdb complex] [receptor chain] [receptor start] [receptor end] [ligand chain] [lignad start] [ligand end]

### set up required files for amber 
RUN_setupRequiredFiles_amber20.sh [pdb complex] [receptor chain] [receptor start] [receptor end]

### copy tleap script to current directory
cp $PROTLIDv1HOME/code_interfaceDesign/template.genInitPDB_wtLeap_oddProbe.sh submit.genInitPDB_wtLeap.sh

### submit tleap script to cluster
qsub -l -t 1-19 submit.genInitPDB_wtLeap.sh

### once tleap files are generated, generate and run MD scripts
cp $PROTLIDv1HOME/code_interfaceDesign/generateSubmitJobs_md.sh .
./generateSubmitJobs_md.sh
for i in `ls submit.1-7_pdb.???.?.sh`; do qsub $i; done

### delete files from MD simulations
python3 $PROTLIDv1HOME/cleanUp.py &

### generate rs-pharmacophore
python3 $PROTLIDv1HOME/get_rs-pharmacophore.py

######################################


### Ligand Search (!!!submits scripts to SLURM!!!)
#######################################
### Run Haddock3
python3 $PROTLIDv1HOME/ligandSearch/setupHaddock3.py [receptor] [cognate ligand]  
#######################################


OUTPUT
Results can be found in [working directory]/ligandSearch/results

REFERENCE
Email Steven for manuscript if not published
TBD
DOI: TBD

AUTHORS
Steven Grudman		steven.grudman@einsteinmed.edu	smgrudman@gmail.com
Christopher McClain	christopher.mcclain@einsteinmed.edu
Andras Fiser		andras.fiser@einsteinmed.edu



