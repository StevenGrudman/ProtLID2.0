# Comparative modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output

# Override the 'select_atoms' routine in the 'automodel' class:
# (To build an all-hydrogen model, derive from allhmodel rather than automodel
# here.)
class MyModel(automodel):

#    def select_atoms(self):
        # Select residues 1 and 2 (PDB numbering)
#        return selection(self.residue_range('119:B', '121:B'))

        # The same thing from chain A (required for multi-chain models):
        # return selection(self.residue_range('1:A', '2:A'))

        # Residues 4, 6, 10:
        # return selection(self.residues['4'], self.residues['6'],
        #                  self.residues['10'])

        # All residues except 1-5:
        # return selection(self) - selection(self.residue_range('1', '5'))

    def special_restraints(self,aln):
        rsr=self.restraints
        rsr.append(file='RECEPTOR_POSITION_RESTRAINT.rsr')
        rsr.append(file='REDUCED_PROBE_RESTRAINT.rsr')

env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl=complete_pdb(env,'temp_centerRes_reduced.pdb') # fill rest of residue from 2 given atoms
mdl.write('temp_centerRes_full.pdb', model_format='PDB')

# selected atoms do not feel the neighborhood
#env.edat.nonbonded_sel_atoms = 2

a = MyModel(env,
            alnfile  = 'temp_alignment.ali',     # alignment filename
            knowns   = ('CENTER_RES', 'REC'),   # codes of the templates
            sequence = 'PROBE_PLUS_RECEPTOR')   # code of the target
#              csrfile =  'POSITION_RESTRAINT.rsr')

a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
a.make(exit_stage=0)                # do the actual comparative modeling
