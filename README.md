## rust_pdbloader
A bunch of Rust snippet to play with pdb. 

This commit was created as a pointer to the state used in CASP14.
A large update is planned.

## Requirement
Rust
https://www.rust-lang.org/

## Usage
git clone https://github.com/yamule/rust_pdbloader.git
cd rust_pdbloader

### Template Based Modelling
cargo run --release --  -residue_mapping -sidechains "./resources/pepbuilderj/resources/sampledresidues/" -backbones "./resources/pepbuilderj/resources/sampledresidues/" -in ./example_files/T1094.1IW7_C.fas -out ./example_files/results/T1094.1IW7_C.fas.mapped.pdb
cargo run --release --  -refinement -toph19 "./resources/toppar_c36_jul18/toppar/toph19.inp" -param19 "./resources/toppar_c36_jul18/toppar/param19.inp" -resource "./resources/" -angle "./resources/angle_distribution_energy.dat" -steps_checkpoint 1 -out ./example_files/results/T1094.1IW7_C.fas.refined1.pdb -param_file example_files/param_refine.txt  -in "./example_files/results/T1094.1IW7_C.fas.mapped.pdb"  -flag "./example_files/results/T1094.1IW7_C.fas.mapped.pdb.flag" -build_missing_param1 "./example_files/param_build_missing.txt"  -build_missing_param2 "./example_files/param_build_missing2.txt"  -num_structurs_step1 5 -num_structurs_step2 1


## Energy function or terms used in this software
solvation_parameters.dat
PROTEINS: Structure, Function, and Genetics 35:133?152 (1999)
Effective Energy Function for Proteins in Solution
Themis Lazaridis1 and Martin Karplus

opt_nov15_lj_param.dat
opt_nov15_partial_charge.dat
Park, Hahnbeom, Philip Bradley, Per Greisen, Yuan Liu, Vikram Khipple Mulligan, David E. Kim, David Baker, and Frank DiMaio. “Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules.” Journal of Chemical Theory and Computation 12, no. 12 (December 13, 2016): 6201–12. https://doi.org/10.1021/acs.jctc.6b00819.

CHARMM
Brooks, B. R. et al. CHARMM: a program for macromolecular energy, minimization, and dynamics calculations. Journal of computational chemistry 4, 187-217 (1983).

EvoEF2
Huang, X., Pearce, R. & Zhang, Y. EvoEF2: accurate and fast energy function for computational protein design. Bioinformatics 36, 1135-1142, doi:10.1093/bioinformatics/btz740 (2020).
