# pymatcher
This python module closely follows the Matcher App and uses the same strategy as the Matcher Mover to output poses instead of PDBs.

The basic usage is :
```
import pymatcher as pm
poselist=pm.matcher(**kwargs)
```
Where kwargs contains the `kwargs["-s"]` with a path to the scaffold or `kwargs["-pose_in"]` with an input pose and `kwargs["-lig_name"]` with the three letter code for the ligand being matched.  All standard match options apply and can be passed as options when initializing pyrosetta.  Please refer to the example notebooks for more detailed usage instructions.

To maintain compatibility with pyrosetta.distributed, the function returns a list of poses, which can be manipulated inside of a python script, in a wrapper function, or passed as individual poses to downstream pyrosetta.distributed protocols.

Current list of kwargs:
```
kwargs["-s"] = scaffold.pdb #input pdb
kwargs["-pose_in"] = scaffold #input pose
kwargs["-scaffold_name"] = "scaffold" #string to be used as the scaffold name.
  Highly recommended with kwargs["-pose_in"].
kwargs["-lig_name"] = "lig" #three letter code for the ligand.  
  The corresponding params file should be passed as a command line argument
kwargs["-match_pos"] = [1,2,3,4] #list of positions to use for matching instead of the pos file.  
  A pos file must still be passed as a command line argument to avoid errors, 
  but this kwarg will overwrite those positions.

