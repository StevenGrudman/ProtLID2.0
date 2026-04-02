# ProtLID2.0  
**Protein-Ligand Interface Design 2.0**

## Overview

ProtLID2.0 operates in **two main steps**:

1. **Residue-specific pharmacophore generation**
2. **Ligand search**

ProtLID2.0 has been shown to be a reliable scoring function for:

- **site-directed protein–protein docking**
- **identifying a receptor’s cognate ligand** from a sub-proteome of decoys

---

## Requirements

The following software is required:

- **Perl** (tested with `v5.16.3`)
- **MODELLER** (tested with `v10.3`)
- **AMBER** (tested with `Amber20`)
- **Python 3** (tested with `v3.9.12`)
- **HADDOCK3** (tested with `2024.12.0b7`)

> **Important:**  
> You must create a Conda environment for HADDOCK3 named:

```bash
haddock3
```

---

## Notes Before You Start

Most commands in this README can be copied and pasted directly.

However, **anything inside square brackets `[]` must be replaced** with your own values.

### Example

```bash
[pdb complex]
```

should be replaced with something like:

```bash
1I8L.pdb
```

---

## Environment Setup

Add the following to your `.bashrc`.

> **Important:**  
> Update the `PROTLIDv2HOME` path and the AMBER path if needed.

```bash
#######################################
export PROTLIDv2HOME=[path to protlidv2.0]
export LD_LIBRARY_PATH=$PROTLIDv2HOME/protlid_auxPrograms/naccess2.1.1:$LD_LIBRARY_PATH
PATH=$PATH:$PROTLIDv2HOME/scripts:$PROTLIDv2HOME/code_interfaceDesign:$PROTLIDv2HOME/code_partnerSearch:$PROTLIDv2HOME/code_interfaceDefinition:$PROTLIDv2HOME/bin:$PROTLIDv2HOME/protlid_auxPrograms:$PROTLIDv2HOME/protlid_auxPrograms/naccess2.1.1
export PATH=$PATH

test -f /usr/local/bio/amber20/amber.sh && source /usr/local/bio/amber20/amber.sh
#######################################
```

Then reload your shell:

```bash
source ~/.bashrc
```

---

# Step 1: Generate the RS-Pharmacophore  
> **This step submits jobs to SGE**

## 1. Create a Working Directory

Create a new directory and copy your PDB complex into it.

> **Important:**  
> All following commands in this section should be run **from this working directory**.

---

## 2. Define the Interface and Mesh Points

```bash
$PROTLIDv2HOME/code_interfaceDefinition/run_genIntAtomList_INTERCAATd4noIllHydRes.sh [pdb complex] [receptor chain] [receptor start] [receptor end] [ligand chain] [ligand start] [ligand end]
```

### Example

```bash
$PROTLIDv2HOME/code_interfaceDefinition/run_genIntAtomList_INTERCAATd4noIllHydRes.sh 1I8L.pdb A 1 106 B 3 120
```

---

## 3. Set Up Required Files for AMBER

```bash
RUN_setupRequiredFiles_amber20.sh [pdb complex] [receptor chain] [receptor start] [receptor end]
```

### Example

```bash
RUN_setupRequiredFiles_amber20.sh 1I8L.pdb A 1 106
```

---

## 4. Copy the `tleap` Submission Script

```bash
cp $PROTLIDv2HOME/code_interfaceDesign/template.genInitPDB_wtLeap_oddProbe_SLURM.sh submit.genInitPDB_wtLeap.sh
```

---

## 5. Submit `tleap` Jobs to the Cluster

```bash
qsub -t 1-19 submit.genInitPDB_wtLeap.sh
```

---

## 6. Generate and Submit MD Jobs

Once the `tleap` files are generated:

```bash
cp $PROTLIDv2HOME/code_interfaceDesign/generateSubmitJobs_md.sh .
./generateSubmitJobs_md.sh
for i in `ls submit.1-7_pdb.???.?.sh`; do qsub $i; done
```

---

## 7. Clean Up MD Simulation Files

```bash
python3 $PROTLIDv2HOME/cleanUp.py &
```

---

## 8. Generate the RS-Pharmacophore

```bash
python3 $PROTLIDv2HOME/get_rs-pharmacophore.py
```

---

# Step 2: Ligand Search  
> **This step submits jobs to SLURM**

## Run HADDOCK3 Setup

```bash
python3 $PROTLIDv2HOME/ligandSearch/setupHaddock3.py [receptor] [cognate ligand]
```

### Example

```bash
python3 $PROTLIDv2HOME/ligandSearch/setupHaddock3.py receptor.pdb ligand.pdb
```

---

## Output

Results will be written to:

```bash
[working directory]/ligandSearch/results
```

---

## Reference

If the manuscript is not yet published, please contact **Steven**.

- **Manuscript:** TBD
- **DOI:** TBD

---

## Authors

- **Steven Grudman**  
  steven.grudman@einsteinmed.edu  
  smgrudman@gmail.com

- **Christopher McClain**  
  christopher.mcclain@einsteinmed.edu

- **Andras Fiser**  
  andras.fiser@einsteinmed.edu
