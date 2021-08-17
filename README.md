# pymatcher
This python module closely follows the Matcher App and uses the same strategy as the Matcher Mover to output poses instead of PDBs.

The basic usage is :
```
import pymatcher as pm
poselist=pm.matcher(**kwargs)
```
Where kwargs contains the `kwargs["-s"]` with a path to the scaffold, and `kwargs["-lig_name"]` with the three letter code for the ligand being matched.  All standard match options apply and can be passed as options when initializing pyrosetta.  Please refer to the example notebooks for more detailed usage instructions.

To maintain compatibility with pyrosetta.distributed, the function returns a list of poses, which can be manipulated inside of a python script, in a wrapper function, or passed as individual poses to downstream pyrosetta.distributed protocols.
