## rust_pdbloader
A bunch of Rust snippet to play with data in PDB (https://www.rcsb.org/). 

These commits were created as a pointer to the state used in CASP14.
A large update is planned.

## Requirement
Rust https://www.rust-lang.org/

## Usage (Template Based Modeling)
``` bash
git clone https://github.com/yamule/rust_pdbloader.git
cd rust_pdbloader
cargo run --release --  -residue_mapping -sidechains "./resources/pepbuilderj/resources/sampledresidues/" -backbones "./resources/pepbuilderj/resources/sampledresidues/" -in "./example_files/T1094.1IW7_C.fas" -out "./example_files/results/T1094.1IW7_C.fas.mapped.pdb"
cargo run --release --  -refinement -resource "./resources/" -angle "./resources/angle_distribution_energy.dat" -steps_checkpoint 1 -out "./example_files/results/T1094.1IW7_C.fas.refined1.pdb" -param_file "example_files/param_refine.txt"  -in "./example_files/results/T1094.1IW7_C.fas.mapped.pdb"  -flag "./example_files/results/T1094.1IW7_C.fas.mapped.pdb.flag" -build_missing_param1 "./example_files/param_build_missing.txt"  -build_missing_param2 "./example_files/param_build_missing2.txt"  -num_structurs_step1 5 -num_structurs_step2 1
```

## Energy function or terms used in this software
> solvation_parameters.dat

Lazaridis, Themis, and Martin Karplus. "Effective energy function for proteins in solution." Proteins: Structure, Function, and Bioinformatics 35.2 (1999): 133-152.


> opt_nov15_lj_param.dat
> opt_nov15_partial_charge.dat

Park, Hahnbeom, et al. "Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules." Journal of chemical theory and computation 12.12 (2016): 6201-6212.


> CHARMM

Brooks, Bernard R., et al. "CHARMM: a program for macromolecular energy, minimization, and dynamics calculations." Journal of computational chemistry 4.2 (1983): 187-217.


> EvoEF2

Huang, Xiaoqiang, Robin Pearce, and Yang Zhang. "EvoEF2: accurate and fast energy function for computational protein design." Bioinformatics 36.4 (2020): 1135-1142.
