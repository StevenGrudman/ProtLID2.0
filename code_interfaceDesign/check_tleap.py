#!/usr/bin/env python3

"""
Purpose:
    Check that ProtLID tleap outputs were generated correctly and are ready for MD simulations.

Expected directory structure:

    project_directory/
    ├── RequiredFiles/
    │   ├── meshindex.*.list
    │   └── mesh.ref.xyz_
    ├── COM_PRM/
    │   └── rec.<AA>.prm
    ├── RST_INIT/
    │   └── <AA>/runid_<N>/rec.<AA>.<meshpoint>.<N>.rst
    └── submit.genInitPDB_wtLeap.sh

What this program does:
    1. Find the required project folders and files.
    2. Find the mesh point list.
    3. Read the tleap submit script to get amino acid probe types and replicate count.
    4. Count mesh points.
    5. Check that required output files exist.
    6. Report missing files and say whether it is safe to continue.
"""

from pathlib import Path
from itertools import product
import sys


def print_expected_structure():
    print("~~~~~~~~~ Expected Project Structure ~~~~~~~~~")
    print("project_directory/")
    print("├── RequiredFiles/")
    print("│   ├── meshindex.*.list")
    print("│   └── mesh.ref.xyz_")
    print("├── COM_PRM/")
    print("│   └── rec.<AA>.prm")
    print("├── RST_INIT/")
    print("│   └── <AA>/runid_<N>/rec.<AA>.<meshpoint>.<N>.rst")
    print("└── submit.genInitPDB_wtLeap.sh")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


def print_resolved_paths(project_dir, required_dir, com_prm_dir, rst_init_dir, tleap_submit, reference_mesh):
    print("~~~~~~~~~ Resolved Paths ~~~~~~~~~")
    print(f"Project directory: {project_dir}")
    print(f"Required files:    {required_dir}")
    print(f"COM_PRM:           {com_prm_dir}")
    print(f"RST_INIT:          {rst_init_dir}")
    print(f"Tleap submit file: {tleap_submit}")
    print(f"Reference mesh:    {reference_mesh}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


def find_meshpoint_file(required_dir):
    mesh_files = list(required_dir.glob("meshindex.*.list"))

    if len(mesh_files) == 0:
        sys.exit(f"ERROR: No meshindex.*.list file found in {required_dir}")

    if len(mesh_files) > 1:
        print("ERROR: Multiple meshindex.*.list files found:")
        for file in mesh_files:
            print(f"    {file}")
        sys.exit("Please keep only one meshindex.*.list file in RequiredFiles.")

    return mesh_files[0]


def read_meshpoints(mesh_file):
    if not mesh_file.is_file():
        sys.exit(f"ERROR: Missing meshpoint file: {mesh_file}")

    with open(mesh_file, "r") as fh:
        meshpoints = [int(line.strip()) for line in fh if line.strip()]

    return meshpoints


def read_tleap_parameters(tleap_submit):
    if not tleap_submit.is_file():
        sys.exit(f"ERROR: Missing tleap submit file: {tleap_submit}")

    aa_probes = []
    num_replicates = None

    with open(tleap_submit, "r") as fh:
        for line in fh:
            line = line.strip()

            if line.startswith("tasks=("):
                aa_probes = line.split("(", 1)[1].rstrip(")").split()

            elif line.startswith("nRun="):
                num_replicates = int(line.split("=", 1)[1])

    if not aa_probes:
        sys.exit("ERROR: Could not find amino acid probe list: tasks=(...)")

    if num_replicates is None:
        sys.exit("ERROR: Could not find replicate count: nRun=...")

    return aa_probes, num_replicates


def print_extracted_parameters(mesh_file, meshpoints, aa_probes, num_replicates):
    print("~~~~~~~~~ Extracted Parameters ~~~~~~~~~")
    print(f"Meshpoint file:              {mesh_file}")
    print(f"Number of mesh points:       {len(meshpoints)}")
    print(f"Number of amino acid probes: {len(aa_probes)}")
    print(f"Number of replicates:        {num_replicates}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


def check_reference_mesh(reference_mesh):
    if not reference_mesh.is_file():
        print(f"ERROR: Missing reference mesh file: {reference_mesh}")
        return False

    print("INFO: Reference mesh file exists")
    return True


def check_prm_files(com_prm_dir, aa_probes):
    ok = True

    for aa in aa_probes:
        prm_file = com_prm_dir / f"rec.{aa}.prm"

        if not prm_file.is_file():
            print(f"ERROR: Missing parameter file: {prm_file}")
            ok = False

    if ok:
        print("INFO: All .prm files are present")

    return ok


def check_rst_files(rst_init_dir, aa_probes, meshpoints, num_replicates):
    ok = True

    for aa, meshpoint, rep in product(
        aa_probes,
        meshpoints,
        range(1, num_replicates + 1),
    ):
        rst_file = rst_init_dir / aa / f"runid_{rep}" / f"rec.{aa}.{meshpoint}.{rep}.rst"

        if not rst_file.is_file():
            print(f"ERROR: Missing restart file: {rst_file}")
            ok = False

    if ok:
        print("INFO: All .rst files are present")

    return ok


def main():
    project_dir = Path.cwd().resolve()
    required_dir = project_dir / "RequiredFiles"
    com_prm_dir = project_dir / "COM_PRM"
    rst_init_dir = project_dir / "RST_INIT"
    tleap_submit = project_dir / "submit.genInitPDB_wtLeap.sh"
    reference_mesh = required_dir / "mesh.ref.xyz_"

    print_expected_structure()
    print_resolved_paths(
        project_dir,
        required_dir,
        com_prm_dir,
        rst_init_dir,
        tleap_submit,
        reference_mesh,
    )

    mesh_file = find_meshpoint_file(required_dir)
    meshpoints = read_meshpoints(mesh_file)
    aa_probes, num_replicates = read_tleap_parameters(tleap_submit)

    print_extracted_parameters(
        mesh_file,
        meshpoints,
        aa_probes,
        num_replicates,
    )

    ref_ok = check_reference_mesh(reference_mesh)
    prm_ok = check_prm_files(com_prm_dir, aa_probes)
    rst_ok = check_rst_files(
        rst_init_dir,
        aa_probes,
        meshpoints,
        num_replicates,
    )

    print()

    if ref_ok and prm_ok and rst_ok:
        print("INFO: All necessary files are present. Proceed to the next step.")
    else:
        print("WARNING: Missing files detected. Resolve issues before continuing.")


if __name__ == "__main__":
    main()
