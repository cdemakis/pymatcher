import pyrosetta
import time
import json
def matcher(**kwargs):
    pyrosetta.rosetta.basic.options.set_boolean_option("run:preserve_header", True)

    poselist = []
    
    matcher_task = pyrosetta.rosetta.protocols.match.MatcherTask()
    
    if kwargs["-s"] is None:
        print("No input scaffolds.  Please use kwargs[-s] to provide input structures.")
        return poselist

    scaffold = pyrosetta.pose_from_pdb(kwargs["-s"])
    scaffold.update_residue_neighbors()

    ligpose = pyrosetta.rosetta.core.pose.Pose()
    chem_mng = (
        pyrosetta.rosetta.utility.SingletonBase_core_chemical_ChemicalManager_t().get_instance()
    )
    ligres = pyrosetta.rosetta.core.conformation.ResidueFactory().create_residue(
    ligpose.conformation()
    .modifiable_residue_type_set_for_conf()
    .name_map(kwargs["-lig_name"])
)
    ligpose.append_residue_by_jump(ligres, 1)

    if pyrosetta.rosetta.basic.options.truefalseoption("match:ligand_rotamer_index"):
        pyrosetta.rosetta.protocols.match.set_ligpose_rotamer(ligpose)

    oats = pyrosetta.rosetta.utility.vector1_core_id_AtomID(3)
    cent, nbr1, nbr2 = ligres.type().select_orient_atoms()
    oats[1] = pyrosetta.rosetta.core.id.AtomID(cent, 1)
    oats[2] = pyrosetta.rosetta.core.id.AtomID(nbr1, 1)
    oats[3] = pyrosetta.rosetta.core.id.AtomID(nbr2, 1)

    matcher_task.set_downstream_pose(ligpose, oats)
    matcher_task.set_upstream_pose(scaffold)

    matcher_task.initialize_from_command_line()
   
    matcher_start_time=time.time()
    matcher = pyrosetta.rosetta.protocols.match.Matcher()
    matcher.initialize_from_task(matcher_task)

    matcher_task.output_writer_name("PoseMatchOutputWriter")
    match_proc = pyrosetta.rosetta.protocols.match.output.ProcessorFactory.create_processor(
        matcher, matcher_task
    )

    matcher.find_hits()
    matcher_end_time=time.time()
    proc_start_time=time.time()
    matcher.process_matches(match_proc)
    
    def get_theozyme_resis(pose):
        ss = pyrosetta.rosetta.std.stringstream()
        pyrosetta.rosetta.core.io.pdb.dump_pdb(pose, ss)
        resis = {}
        for line in ss.str().splitlines():
            if "REMARK 666" in line:
                resis[line.strip().split()[-3]] = line.strip().split()[-4]
    
        return resis

    if match_proc.output_writer().match_groups_ushits().u()>=1:
        for m in range(1, match_proc.output_writer().match_groups_ushits().u() + 1):
            matchedpose = scaffold.clone()
            match_proc.output_writer().insert_match_into_pose(matchedpose, m)
            theozyme=get_theozyme_resis(matchedpose)
            matcher_data=json.dumps({"match_group":m, "theozyme":theozyme, "ligand":kwargs["-lig_name"], "scaffold":kwargs["-s"], "cst_file":pyrosetta.rosetta.basic.options.get_file_option("match:geometric_constraint_file")})
            pyrosetta.rosetta.core.pose.setPoseExtraScore(matchedpose,"matcher_data",matcher_data)
            poselist.append(matchedpose)
    else:
        print("Matcher found no matches")
        
    proc_end_time=time.time()

    print(f"Matcher ran for {proc_end_time - matcher_start_time} seconds, where finding hits took {matcher_end_time - matcher_start_time} seconds and processing the matches took {proc_end_time - proc_end_time} seconds.")
    return poselist

