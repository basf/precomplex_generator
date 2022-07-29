### Openbabel incompatibility

``` Python
# import openbabel as ob
from openbabel import openbabel as ob

# symbol1 = ob.OBElementTable().GetSymbol(atom.GetAtomicNum())
symbol1 = element_dict[atom.GetAtomicNum()]
```

- Requires `element_dict`, as defined in file element_dict

### Fortran linking of shared object files

- Copied all the required functions/subroutines into `final.f90` and `get_all_torsions_mod.f90` to avoid linking errors when compiling with f2py

### Relative imports

- Removed precomplex_generator specification for relative imports from working directory:

``` Python
# from precomplex_generator.tools import file_read_dict
from tools import file_read_dict
```

### Misc

``` Python
# from StructureGeneration.symmetry_tools import * -> Unresolved input within try-except block
```

- `Esterification_Cl` still making problems -> Some pymatgen error. First precomplex built, though...
- Same error with Bischler-Napieralski, as well as Ester-Hydrolysis
-
