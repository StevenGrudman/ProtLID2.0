# Comparative modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output

env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl=complete_pdb(env,'temp_centerRes_reduced.pdb') # fill rest of residue from 2 given atoms
mdl.write('temp_centerRes_full.pdb', model_format='PDB')

# selected atoms do not feel the neighborhood
#env.edat.nonbonded_sel_atoms = 2

a = automodel(env,
              alnfile  = 'temp_alignment.ali',     # alignment filename
              knowns   = ('CENTER_RES', 'REC'),   # codes of the templates
              sequence = 'PROBE_PLUS_RECEPTOR')   # code of the target

a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
a.make(exit_stage=2)                # do the actual comparative modeling, exit w/o optimization

# output: pdbfile PROBE_PLUS_RECEPTOR.ini
