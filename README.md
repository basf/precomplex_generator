# Overview

The precomplex generator is a tool for obtaining suitable input structures for automated transition-state (TS) searches (precomplexes) based on single-ended reaction path optimization algorithms. Our strategy is applicable to uni- and bimolecular reactions as well as pseudo-termolecular reaction for a broad range of molecular structures. The precomplex generator rotates around selected bonds to obtain a structure in which the reaction centers are well aligned.


# Requirements:

```
- Python3.6 
- OpenBabel
- numpy 
- Cython
- wheel
- Networkx
- Pymatgen
```
# How to directly use the precomplex generator:
We provide precompiled Fortran code and Python files, which can be directly used if local library compatibility is given (RHEL7 or similar OS).

```
1. Clone the Github repository to your current directory ($PWD):
   git clone https://github.com/basf/precomplex_generator.git

2. Add the PYTHONPATH to your PYTHONPATH's:
   export PYTHONPATH=$PWD/precomplex_generator:$PYTHONPATH 

3. Add alias for convenience:
   export top_dir=$PWD 
   alias prec_gen='python3 $top_dir/precomplex_generator/precomplex_generator/main.py'

```
# How to install the precomplex generator (.egg):
```
1. Clone the Github repository to your current directory ($PWD):
   git clone https://github.com/basf/precomplex_generator.git
   
2. For installation, we provide an install.sh script. 
   Caution: This script will only work, if python and openbabel are already available!
   cd precomplex_generator
   source install.sh
```

# How to run the precomplex_generator:
```
You find examples in the EXAMPLES folder.

Only the reactant XYZ file(s) as well as the rule.txt file are required.
The rule.txt file contains the "reaction rule" in terms of atom-specific bond formations/breakages.
The COMMAND file contains the respective command you need.

For example:
prec_gen mol_name_1 mol_name_2 -r rule.txt
```

