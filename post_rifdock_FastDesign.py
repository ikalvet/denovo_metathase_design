#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:19:22 2019

@author: ikalvet
"""

import pyrosetta as pyr
import pyrosetta.rosetta
import os
import time
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.protocols.grafting.simple_movers import DeleteRegionMover
from pyrosetta.rosetta.core.select import get_residues_from_subset
import sys
import numpy as np
from shutil import copy2
import argparse

SCRIPT_DIR = os.path.dirname(__file__)

polar_residues = ['HIS', 'LYS', 'SER', 'THR', 'TYR', 'GLN', 'GLU',
                  'ASN', 'ASP', 'ARG', 'CYS']

aromatic_residues = ['PHE', 'TRP']

cuts = [6.0, 8.0, 10.0, 12.0]

cst_B, cst_A = np.polyfit(np.array([0.01, 0.5]), np.log(np.array([0.01, 10.0])), 1)
cst_A = np.e**cst_A


def pymol_string(list):
    return "+".join([str(x) for x in list])


def get_rif_residues(pdbfile):
    """
    Reads the PDB file as a text to find
    which residues were introduced by rifdock
    """
    pdblines = open(pdbfile, 'r').readlines()

    # Extracting RIF residues from PDBfile
    rifres = None
    for line in pdblines:
        if 'rif_residues' in line:
            rifreslinesplit = line.split()
            if rifreslinesplit.__len__() == 1:
                continue
            else:
                rifres = rifreslinesplit[1].split(',')
                rifres = [int(x) for x in rifres]
            break
    return rifres


def get_ligand_heavyatoms(pose):
    """
    Returns a list of non-hydrogen atomnames
    """
    if pose.residue(pose.size()).is_ligand() is False:
        print("Last residue is not ligand!")
        return None
    heavyatoms = []
    for n in range(1, pose.residue(pose.size()).natoms() + 1):
        element = pose.residue(pose.size()).atom_type(n).element()
        if element != 'H':
            heavyatoms.append(pose.residue(pose.size()).atom_name(n).lstrip().rstrip())
    return heavyatoms


def get_packer_layers(pose, target_atoms, do_not_design, cuts):
    """
    Finds residues that are within certain distances from target atoms, defined though <cuts>.
    Returns a list of lists where each embedded list contains residue numbers belongin to that layer.
    Last list contains all other residue numbers that were left over.
    get_packer_layers(pose, target_atoms, do_not_design, cuts) -> list
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose)
        target_atoms (list) :: list of integers.
                               Atom numbers of a residue or a ligand that are used to calculate residue distances.
        do_not_design (list) :: list of integers. Residue numbers that are not allowed to be designed.
                                Used to override layering to not design matched residues.
        cuts (list) :: list of floats.
    """
    assert len(cuts) > 3, f"Not enough layer cut distances defined {cuts}"
    print("Determining design/repack/do-not-touch layers "
          "based on cuts: {}".format(cuts))

    residues = []
    ligand_atoms = []
    for a in target_atoms:
        ligand_atoms.append(pose.residue(pose.size()).xyz(a))

    for i in cuts:
        residues.append([])
    residues.append([])  # Additional list for all other residues that will be kept static

    for resno in range(1, pose.size()):
        if resno not in do_not_design:
            resname = pose.residue(resno).name()
            CA = pose.residue(resno).xyz('CA')
            CA_distances = []

            if 'GLY' not in resname:
                CB = pose.residue(resno).xyz('CB')
                CB_distances = []

            for a in ligand_atoms:
                CA_distances.append((a - CA).norm())
                if 'GLY' not in resname:
                    CB_distances.append((a - CB).norm())
            CA_mindist = min(CA_distances)

            if 'GLY' not in resname:
                CB_mindist = CB_distances[CA_distances.index(CA_mindist)]

            # Figuring out into which cut that residue belongs to,
            # based on the smallest CA distance and whether the CA is further away from the ligand than CB or not.
            # PRO and GLY are disallowed form the first two cuts that would allow design,
            # they can only be repacked.
            if CA_mindist <= cuts[0] and resname not in ['PRO', 'GLY']:
                residues[0].append(resno)
            elif CA_mindist <= cuts[1] and resname not in ['PRO', 'GLY'] and CB_mindist < CA_mindist:
                residues[1].append(resno)
            elif CA_mindist <= cuts[2]:
                residues[2].append(resno)
            elif CA_mindist <= cuts[3] and CB_mindist < CA_mindist:
                residues[3].append(resno)
            else:
                residues[-1].append(resno)
    return residues


def get_layer_selections(pose, rifres, constraints, heavyatoms):
    # Getting the redesign / repack / do-not-touch layers for FastDesign
    # (It is not possible to create these layers using ResidueSelectors the way EnzDes does it)
    residues = get_packer_layers(pose, heavyatoms, rifres, cuts)

    # Creating a ResidueSelector for repack layer
    SEL_repack_residues = ResidueIndexSelector()
    for res in residues[2] + residues[3] + list(constraints.keys()):
        SEL_repack_residues.append_index(res)

    # Creating a ResidueSelector for do-not-touch layer
    SEL_do_not_repack = ResidueIndexSelector()
    for res in residues[4]:
        SEL_do_not_repack.append_index(res)

    # Creating a list of residues that will be mutated.
    # Including RIF residues that ended up not being used in the constraints
    # This means that almost all apolar RIF residues will be redesigned
    mutate_residues = residues[0] + residues[1] +\
                        [x for x in rifres if x not in constraints]

    SEL_mutate_residues = ResidueIndexSelector()
    for res in mutate_residues:
        SEL_mutate_residues.append_index(res)
    return SEL_mutate_residues, SEL_repack_residues, SEL_do_not_repack, residues


def separate_protein_and_ligand(pose):
    """
    Separates the ligand and the rest of the pose to 'infinity'
    (to 666 angstrom).
    Assumes ligand is the last residue in the pose.
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose)
    """
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"

    tmp_pose = pose.clone()
    lig_seqpos = tmp_pose.size()
    lig_jump_no = tmp_pose.fold_tree().get_jump_that_builds_residue(lig_seqpos)
    rbt = pyr.rosetta.protocols.rigid.RigidBodyTransMover(tmp_pose, lig_jump_no)
    rbt.step_size(666)
    rbt.apply(tmp_pose)
    return tmp_pose


def mutate_residues_to_ala(pose, mutate_residues):
    """
    Mutates given residues to alanines.
    Arguments:
        pose (obj, pyrosetta.rosetta.core.pose.Pose)
        mutate_residues (obj, pyrosetta.rosetta.core.select.residue_selector.ResidueSelector) :: list of integers
    """
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"
    assert isinstance(mutate_residues, pyrosetta.rosetta.core.select.residue_selector.ResidueSelector), "mutate_residues: Invalid input type"
    # assert all([isinstance(x, int) for x in mutate_residues]), "mutate_residues: expected a list of integers"

    pose2 = pose.clone()
    mutres = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
    mutres.set_selector(mutate_residues)
    mutres.set_res_name('ALA')
    mutres.apply(pose2)
    return pose2


def run_fastdesign(pose, scorefunction, mutate_residues=None,
                   repack_residues=None, do_not_repack=None, cst=None,
                   cartesian=False, mutate_ALA=False, rifres=None):
    """
    Performs FastDesign on an input pose.
    Constraints can be applied through the 'cst' argument
    Repack-only residues are defined through the 'repack_residues' argument using a ResidueSelector
    Do-not-repack residues are defined through the 'do_not_repac'k argument using a ResidueSelector
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose)
        scorefunction (obj, pyrosetta.rosetta.core.scoring.ScoreFunction)
        repack_residues (obj, pyrosetta.rosetta.core.select.residue_selector.ResidueSelector)
        do_not_repack (obj, pyrosetta.rosetta.core.select.residue_selector.ResidueSelector)
        cst (obj, pyrosetta.rosetta.protocols.constraint_generator.AddConstraints)
        cartesian (bool)
        mutate_ALA (bool)
    """
    nt = type(None)
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"
    assert isinstance(scorefunction, pyrosetta.rosetta.core.scoring.ScoreFunction), "scorefunction: Invalid input type"
    assert isinstance(mutate_residues, (pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, nt)), "mutate_residues: Invalid input type"
    assert isinstance(repack_residues, (pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, nt)), "repack_residues: Invalid input type"
    assert isinstance(do_not_repack, (pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, nt)), "do_not_repack: Invalid input type"
    assert isinstance(cst, (pyrosetta.rosetta.protocols.constraint_generator.AddConstraints, nt)), "cst: Invalid input type"
    assert isinstance(cartesian, bool), "cartesian: Invalid input type"
    assert isinstance(mutate_ALA, bool), "mutate_ALA: Invalid input type"

    pose2 = pose.clone()

    ligand = ResidueIndexSelector()
    ligand.append_index(pose2.size())

    # If requested, all designable residues are first mutated to alanines
    if mutate_ALA:
        assert mutate_residues is not None, "Need to define mutable residues!"
        pose2 = mutate_residues_to_ala(pose2, mutate_residues)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Setting up residue-level extra rotamers. Not sure if better than standard?
    erg_RLT = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGenericRLT()
    erg_RLT.ex1(False)
    erg_RLT.ex2(False)

    if mutate_residues is not None:
        DIN = pyrosetta.rosetta.core.pack.task.operation.DisallowIfNonnativeRLT()
        DIN.disallow_aas("MCGP")
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(DIN, mutate_residues))
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(erg_RLT, mutate_residues))

    # if rifres is not None:
    #     SEL_rifres = ResidueIndexSelector()
    #     for res in rifres:
    #         SEL_rifres.append_index(res)
    #     restr_RLT = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
    #     tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(restr_RLT, SEL_rifres))
    #     tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(erg_RLT, SEL_rifres))

    if repack_residues is not None:
        # disable design on the only repack part
        restr_RLT = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(restr_RLT, repack_residues))
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(erg_RLT, repack_residues))

    if do_not_repack is not None:
        # don't repack the rest of the protein
        prvnt_RLT = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prvnt_RLT, do_not_repack, False))

    packer_task = tf.create_task_and_apply_taskoperations(pose2)

    # Set up a MoveMap
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(True)

    # Turning off bb and chi torsions for the ligand:
    mm.set_bb(pose2.size(), False)
    mm.set_chi(pose2.size(), False)

    # Applying constraints, if defined
    if cst is not None:
        cst.apply(pose2)

    # Setting up FastDesign
    fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(scorefxn_in=scorefunction, standard_repeats=1)
    fast_design.cartesian(cartesian)
    fast_design.set_task_factory(tf)
    fast_design.set_movemap(mm)
    if cartesian:
        fast_design.min_type("lbfgs_armijo_nonmonotone")  # For Cartesian scorefunctions
    elif cartesian is False:
        fast_design.min_type("dfpmin_armijo_nonmonotone")  # For non-Cartesian scorefunctions
    fast_design.apply(pose2)
    return pose2


def no_ligand_repack(pose, scrfxn, repack_residues_list=None, outfile=None):
    """
    Performs repacking of an input pose.
    If repack_residues_list is provided then it generates a ResidueSelector
    out of that list, allowing only these residues to be repacked.
    All other residues are nto touched.
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose)
        scrfxn (obj, pyrosetta.rosetta.core.scoring.ScoreFunction)
        repack_residues_list (list) :: list of integers
        outfile (str) :: name of the dumped PDB file
    """
    nt = type(None)
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"
    assert isinstance(scrfxn, pyrosetta.rosetta.core.scoring.ScoreFunction), "scrfxn: Invalid input type"
    assert isinstance(repack_residues_list, (list, nt)), "repack_residues_list: Invalid input type"
    assert isinstance(outfile, (str, nt)), "outfile: Invalid input type"

    # Building a ResidueSelector out of a list of residue numbers
    # that should be repacked
    if repack_residues_list is not None:
        repack_residues = ResidueIndexSelector()
        for r in repack_residues_list:
            repack_residues.append_index(r)
        do_not_repack = NotResidueSelector(repack_residues)
    else:
        # If no list of resno's is provided, exit!
        sys.exit("Not implemented yet!")

    tmp_pose = pose.clone()

    # Now lets see how to set up task operations
    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These three are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Setting up residue-level extra rotamers. Not sure if better than standard?
    erg_RLT = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGenericRLT()
    erg_RLT.ex1(True)
    erg_RLT.ex2(True)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(erg_RLT, repack_residues))

    # We are not designing, so this will also be used
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())

    # disable design on the repack-only part
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repack_residues))

    # don't repack the rest of the protein
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), do_not_repack, False))

    # Have to convert the task factory into a PackerTask
    task = tf.create_task_and_apply_taskoperations(tmp_pose)

    pack_mover = PackRotamersMover(scrfxn, task)
    pack_mover.apply(tmp_pose)

    if outfile is not None:
        tmp_pose.dump_pdb(outfile)

    return tmp_pose


def get_constraintpairs(pose, residues, target_resno, target_atomnames):
    """
    Generates constraintpair dictionary out of a pose, list of residue numbers
    and a list of target atomnames.
    Each item in the dictionary represents a residue-target pair.
    Each pair is described by the interacting residue atom name,
        residue atom number, target_atom name, interatomic distance.
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose)
        residues (list) :: list of integers
        target_resno (int)
        target_atomnames (list) :: list of strings
    """
    TRP_TYR_polar_atoms = {'TRP': ['NE1', 'N', 'O'],  # TODO: also include N and C-term
                           'TYR': ['OH', 'N', 'O'],  # TODO: also include N and C-term
                           'PHE': ['N', 'O']}  # TODO: also include N and C-term

    constraintpairs = {}

    for res in residues:
        dists = {}
        for n in range(1, pose.residue(res).natoms()+1):
            dists[n] = {}
            if not pose.residue(res).atom_is_hydrogen(n):
                for trgt in target_atomnames:
                    dists[n][trgt] = (pose.residue(target_resno).xyz(trgt) -
                                      pose.residue(res).xyz(n)).norm()

        closest_res_atom = None
        closest_tgt_atom = None
        shortest_dist = 1000.0

        for n in dists.keys():
            for trgt in dists[n].keys():
                if dists[n][trgt] < shortest_dist:
                    shortest_dist = dists[n][trgt]
                    closest_res_atom = n
                    closest_tgt_atom = trgt

        closest_res_atom_name = pose.residue(res).atom_name(closest_res_atom).rstrip().lstrip()

        # Checking if TYR or TRP is interacting through its polar heteroatom
        # and with a heteroatom of a ligand
        cst_pi_stack = False
        if pose.residue(res).name3() in ['TRP', 'TYR', 'PHE']:
            if closest_res_atom_name in TRP_TYR_polar_atoms[pose.residue(res).name3()]:
                if pose.residue(res).atom_is_backbone(closest_res_atom):
                    cst_pi_stack = False
                elif pose.residue(target_resno).atom_type(target_atomnames.index(closest_tgt_atom)+1).element() != 'C':
                    cst_pi_stack = False
                else:
                    cst_pi_stack = True
            else:
                cst_pi_stack = True

        # If TRP, PHE or TYR is found to be pi-stacking then ligand aromatic
        # carbons are used for constraints. Two constraints are generated
        # to ensure that the relative orientation of rings is maintained.
        if cst_pi_stack:
            print("Residue {}: "
                  "recalculating constraints for pi stacking".format(res))
            interactions = []
            for n in dists.keys():
                for trgt in dists[n].keys():
                    trgt_no = target_atomnames.index(trgt) + 1
                    trgt_type_name = pose.residue(target_resno).atom_type(trgt_no).atom_type_name()
                    res_atom_name = pose.residue(res).atom_name(n).rstrip().lstrip()
                    interactions.append([n, trgt, dists[n][trgt],
                                         trgt_type_name, res_atom_name])
            sorted_interactions = sorted(interactions, key=lambda x: x[2])

            count = 0
            pi_stack_csts = []
            for i in sorted_interactions:
                if i[3] == 'aroC':
                    if count == 1:
                        if i[1] == pi_stack_csts[0]['target_atom'] or i[-1] == pi_stack_csts[0]['res_atom_name']:
                            continue
                    pi_stack_csts.append({'res_atom_no': i[0],
                                          'res_atom_name': i[-1],
                                          'target_atom': i[1],
                                          'dist': i[2]})

                    count += 1
                if count == 2:
                    break

            if any([cst['dist'] > 4.0 for cst in pi_stack_csts]):
                print("Residue {}: pi stacking not constrained, "
                      "residue too far.".format(res))
                continue

            constraintpairs[res] = pi_stack_csts
        else:
            tgt_atom_no = pose.residue(target_resno).atom_index(closest_tgt_atom)
            if any([pose.residue(x).atom_type(y).element() == 'C' for x, y in zip([res, target_resno], [closest_res_atom, tgt_atom_no])]):
                print("Residue {}{}: found close atoms {}-{}. "
                      "One of them is carbon - not constraining this "
                      "interaction.".format(pose.residue(res).name1(), res,
                                            closest_res_atom_name,
                                            closest_tgt_atom))
                continue

            constraintpairs[res] = [{'res_atom_no': closest_res_atom,
                                     'res_atom_name': closest_res_atom_name,
                                     'target_atom': closest_tgt_atom,
                                     'dist': shortest_dist}]

    return constraintpairs


def build_cst_xml_obj(pose, constraints, score_wts="beta", weight=10.0,
                      function='HARMONIC', width=0.5):
    """
    This function generates an XML object from string.
    The purpose of this XML object is to create a DistanceConstraintGenerator
    for PyRosetta.
    It takes a pose and constraint dictionary as input.
    The constraints are created between the ligand (last residue in the pose)
    and the residues that are keys in the constraints dictionary.
    Arguments:
        pose (pyrosetta pose object)
        constraints (dict)
        weight (float)
        function (str)
        width (float)
    """

    str1 = """<SCOREFXNS><ScoreFunction name="beta" weights="{}">
    			<Reweight scoretype="atom_pair_constraint" weight="{}" />
    		    </ScoreFunction>
                </SCOREFXNS>
                <RESIDUE_SELECTORS>""".format(score_wts, weight)

    ligand_sel = "<Index name=\"ligand\" resnums=\"{}\" />".format(pose.size())

    for n, resno in enumerate(constraints.keys()):
        res_sel_line = "<Index name=\"r{}\" resnums=\"{}\" />".format(n+1, resno)
        res_neighbor_line = ''
        str1 += res_sel_line + res_neighbor_line

    str3 = """</RESIDUE_SELECTORS>
                <MOVERS>
                <AddConstraints name="add_cst">"""
    for n, resno in enumerate(constraints.keys()):
        for cst in constraints[resno]:
            str3 += """<DistanceConstraintGenerator
                        name="test{0}"
                        residue_selector1="ligand"
                        residue_selector2="r{0}"
                        atom_name1="{1}"
                        atom_name2="{2}"
                        function="{4} {3:.1f} {5:.1f}" />""".format(n+1, cst['target_atom'],
                                                                    cst['res_atom_name'],
                                                                    cst['dist'],
                                                                    function, width)

    str3 += """</AddConstraints>
            </MOVERS>
            <PROTOCOLS>
            <Add mover="add_cst" />
            </PROTOCOLS>"""

    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(str1 + ligand_sel + str3)

    return obj


def print_scorefile_line(pose, name, additional_scores=None):
    """
    Generates a formatted line for scorefile entry.
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose) :: scored pose
        name (str) :: name of the pose object
        additional_scores (dict) :: additional scores not in the pose object
    """
    scoreline = ""
    scoreline += "{:>14.3f}".format(pose.scores['total_score'])
    scoreline += "{:>14.2f}".format(pose.scores['total_score']/pose.size())
    for k in pose.scores.keys():
        if k != 'total_score':
            labelwidth = 14
            if len(k) >= labelwidth:
                labelwidth = len(k) + 1
            scoreline += "{0:>{1}.2f}".format(pose.scores[k], labelwidth)
    if additional_scores is not None:
        for k in additional_scores.keys():
            labelwidth = 14
            if len(k) >= labelwidth:
                labelwidth = len(k) + 10
            if not isinstance(additional_scores[k], (float, int)):
                while len(additional_scores[k]) >= labelwidth:
                    labelwidth += 1
            if isinstance(additional_scores[k], float):
                scoreline += "{0:>{1}.2f}".format(additional_scores[k], labelwidth)
            else:
                if isinstance(additional_scores[k], str):
                    while len(additional_scores[k]) >= labelwidth:
                        labelwidth += 1
                scoreline += "{0:>{1}}".format(additional_scores[k], labelwidth)
    scoreline += "{0:>{1}}\n".format(name, len(name)+4)
    return scoreline


def print_scorefile_header(pose, name, additional_scores=None):
    """
    Generates a formatted header line for the scorefile.
    Arguments:
        pose (object, pyrosetta.rosetta.core.pose.Pose) :: scored pose
        name (str) :: name of the pose object
        additional_scores (dict) :: additional scores not in the pose object
    """
    nt = type(None)
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"
    assert isinstance(name, str), "name: Invalid input type"
    assert isinstance(additional_scores, (dict, nt)), "additional_scores: Invalid input type"
    
    scoreheader = ""
    scoreheader += "{:>14}".format('total_score')
    scoreheader += "{:>14}".format('score_per_res')
    for k in pose.scores.keys():
        if k != 'total_score':
            labelwidth = 14
            if len(k) >= labelwidth:
                labelwidth = len(k) + 1
            scoreheader += "{0:>{1}}".format(k, labelwidth)
    if additional_scores is not None:
        for k in additional_scores.keys():
            labelwidth = 14
            if len(k) >= labelwidth:
                labelwidth = len(k) + 10
            if not isinstance(additional_scores[k], (float, int)):
                if labelwidth < len(additional_scores[k]):
                    labelwidth = len(additional_scores[k]) + 10
            scoreheader += "{0:>{1}}".format(k, labelwidth)
    scoreheader += "{0:>{1}}\n".format('description', len(name)+10)
    return scoreheader


def cst_score_from_value(cst_value, measured_value):
    """
    Returns the constraint score for a given value
    using an exponential function until error of 0.5, and logarithmic after that
    """
    delta = abs(cst_value - measured_value)
    if delta <= 0.5:
        score = cst_A * np.e**(cst_B * delta)
    else:
        score = 10.0 + np.log(delta/0.5)*10
    return score


def calculate_cst_score(pose, constraints):
    cst_score = 0.0
    for res in constraints:
        for cst in constraints[res]:
            resa = cst['res_atom_name']
            tgt = cst['target_atom']
            cst_value = cst['dist']
            measured_value = (pose.residue(pose.size()).xyz(tgt) -
                              pose.residue(res).xyz(resa)).norm()
            cst_score += cst_score_from_value(cst_value, measured_value)
    return cst_score


def getSASA(pose, resno=None):
    """
    Takes in a pose and calculates its SASA.

    If <resno> is provided then the SASA of that residue is returned.
    then only the SASA of the ligand inside the pocket is returned.

    Procedure by Brian Coventry
    """
    assert isinstance(pose, pyrosetta.rosetta.core.pose.Pose), "pose: Invalid input type"

    atoms = pyrosetta.rosetta.core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())

    for i in range(1, pose.size()+1):
        atoms.resize(i, pose.residue(i).natoms(), True)

    surf_vol = pyrosetta.rosetta.core.scoring.packing.get_surf_vol(pose, atoms, 1.4)

    if resno is not None:
        res_surf = 0
        for i in range(1, pose.residue(resno).natoms()+1):
            res_surf += surf_vol.surf(resno, i)
        return res_surf
    else:
        return surf_vol


def rmsd_no_super(pose1, pose2, residues):
    """
    This function is supposed to mimic the corresponding function found in Rosetta.
    """
    sum2 = 0.0
    natoms = 0
    for res in residues:
        num_atoms = pose1.residue(res).natoms()
        for atomno in range(1, num_atoms+1):
            diff = pose1.residue(res).xyz(atomno) - pose2.residue(res).xyz(atomno)
            sum2 += diff.length_squared()
            natoms += 1
    return np.sqrt(sum2/natoms)


def trash_pdbfile(pdbfile):
    if not os.path.exists('bad_docks'):
        os.mkdir('bad_docks')

    copy2(pdbfile, f'bad_docks/{pdbfile}')



"""
Parsing user input
"""

def main(args):

    pdbfile = args.pdb
    repack = args.repack
    repack_dump_pdb = args.repack_dump_pdb
    outdir = args.outdir
    nstruct = args.nstruct
    suffix = ""
    if args.suffix is not None:
        suffix += f"_{args.suffix}"
    ala_design = args.no_ala_design
    min_polar = args.min_polar
    dump_debug_pdb = args.dump_intermediate_pdb

    if args.weights is not None:
        # args.weights = f"-score:weights {args.weights}"
        score_wts = f"-score:weights {args.weights}"

    test_run = args.test
    cst_weights = args.ramp_cst_weights
    extra_res_fa = ""
    if args.params is not None:
        extra_res_fa = "-extra_res_fa "
        extra_res_fa += " ".join(args.params)


    # Reporting user input:
    print("\nInput parameters/settings:")
    print("  Input PDB: " + pdbfile)
    if args.params is not None:
        print("  Used params files: " + " ".join(args.params))
    print("  Output directory: " + os.getcwd() + '/' + outdir)
    print("  Output suffix: " + suffix)
    print("  Number of design iterations: {}".format(nstruct))
    print("  No-ligand-repack performed: {}".format(repack))
    print("  No-ligand-repack output PDB generated: {}".format(repack_dump_pdb))
    print(f"  Designable positions will be mutated to ALA before design: {ala_design}")
    print(f"  Minimum polar RIF residues required to start design: {min_polar}")
    print("  Constraint weights used: {}".format(', '.join([str(c) for c in cst_weights])))
    print(f"  Intermediate PDBs will be dumped (debug mode): {dump_debug_pdb}\n")
    print(f"  Scoring weights file used: {score_wts}\n")

    """
    Getting Rosetta started
    """
    print("")
    DAB = args.dalphaball
    assert os.path.exists(DAB), "Please provide correct path to DAlphaBall.gcc"
    pyr.init(f"-dalphaball {DAB} {extra_res_fa} -beta {score_wts} -corrections::beta_nov16")

    scorefxn = pyr.get_fa_scorefxn()

    sc = pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter()


    pose = pyr.pose_from_file(pdbfile)


    """
    Calculating the ligand SASA
    """
    # Creating a pose object that just contains the ligand:
    ligand_pose = pose.split_by_chain()[pose.num_chains()]
    ligand_sasa = getSASA(ligand_pose).tot_surf

    # Finding what is the RESNO of the ligand
    # Can't always assume it's the last residue of the pose
    ligand_seqpos = None
    for res in pose.residues:
        if res.is_ligand():
            ligand_seqpos = res.seqpos()
            break
    if ligand_seqpos is None:
        sys.exit("Can't find ligand in the pose. "
                 f"Are you sure you have the correct PDB: {pdbfile}?")

    # Calculating the SASA of the ligand in the docked pose
    # If too exposed, do not design!
    docked_ligand_surf = getSASA(pose, resno=ligand_seqpos)

    rel_sasa = docked_ligand_surf/ligand_sasa

    if rel_sasa > 0.20:
        trash_pdbfile(pdbfile)
        sys.exit("{}: ligand is too exposed. "
                 "SASA = {:.2f}".format(pdbfile, rel_sasa))


    """
    Sorting out residues:
        - which rif residues are there
        - which polar rif residues are there
        - what are the names of ligand heavyatoms
        - what are the redesign / repack / do-not-touch layers for FastDesign
        - creating ResidueSelectors for FastDesign
        - finding polar contact constraints for FastDesign
    """
    rifres = get_rif_residues(pdbfile)

    if rifres is None:
        trash_pdbfile(pdbfile)
        sys.exit("{}: No rif residues found in pdb file.".format(pdbfile))


    # Only keeping polar rif residues and PHE.
    polar_rifres = []
    keep_rifres = []
    for res in rifres:
        if pose.residue(res).name3() in polar_residues + aromatic_residues:
            keep_rifres.append(res)
        if pose.residue(res).name3() in polar_residues:
            polar_rifres.append(res)

    if len(polar_rifres) < min_polar:
        print("{}: Not enough polar rif residues found "
              "({} < {})".format(pdbfile, len(polar_rifres), min_polar))
        sys.exit(1)


    """
    Setting up constraints between the ligand and the protein
    """
    # Getting the names of ligand heavyatoms (non-H)
    heavyatoms = get_ligand_heavyatoms(pose)

    if heavyatoms is None:
        sys.exit("No heavyatoms found for ligand. You messed up big time!?")

    # Getting constraint definitions (as a dictionary) between the ligand and RIF residues
    # Constraints are set only between heavyatoms!
    print("Calculating distance constraints based on polar (and aromatic) RIF residues.")
    constraintpairs = get_constraintpairs(pose, keep_rifres, ligand_seqpos, heavyatoms)

    print("\n{} constraints used during the design:".format(sum([len(constraintpairs[k]) for k in constraintpairs.keys()])))
    for resno in constraintpairs.keys():
        for cst in constraintpairs[resno]:
            print("  ligand-{}{}: atoms {}-{}, d = {:.1f}".format(pose.residue(resno).name1(),
                                                                  resno,
                                                                  cst['target_atom'],
                                                                  cst['res_atom_name'],
                                                                  cst['dist']))

    # Building an XML object that only contains the constraints
    # It has to be done this way because PyRosetta does not have
    # DistanceConstraintGenerator mover properly implemented
    obj = build_cst_xml_obj(pose, constraintpairs, score_wts=args.weights, weight=cst_weights[0])

    # Extracting the constraints mover from the XML object
    cst = obj.get_mover('add_cst')


    # Extracting the constraint-weighed scorefunction form the XML object
    scorefxn_cst = obj.get_score_function('beta')


    """
    Setting up layers
    """

    SEL_mutate_residues, SEL_repack_residues, SEL_do_not_repack, residues = get_layer_selections(pose, keep_rifres, constraintpairs, heavyatoms)

    mutate_residues = get_residues_from_subset(SEL_mutate_residues.apply(pose))
    print(f"{len(mutate_residues)} residues will be redesigned:")
    print("  {}".format('+'.join([str(x) for x in mutate_residues])))
    print("{} residues will be repacked:".format(len(residues[2] + residues[3] + list(constraintpairs.keys()))))
    print("  {}".format(pymol_string(get_residues_from_subset(SEL_repack_residues.apply(pose)))))
    print(f"{len(residues[4])} residues will not be touched:")
    print("  {}".format(pymol_string(get_residues_from_subset(SEL_do_not_repack.apply(pose)))))

    # Finding the closest contacts between polar RIF residues and the ligand
    # These contacts are in the form of atompairs that are used to generate
    # constraints for FastDesign
    print("Polar RIF residues found: {}".format('+'.join([str(x) for x in polar_rifres])))


    # Getting no-ligand-repack residues.
    # Basically a list of all residues that are close-ish to the ligand.
    nlr_repack_residues = residues[0] + residues[1] + residues[2] + residues[3] + keep_rifres


    if dump_debug_pdb:
        pose_ala = mutate_residues_to_ala(pose, SEL_mutate_residues)
        pose_ala.dump_pdb(outdir + pdbfile.replace('.pdb', '_fdes{}_ala.pdb'.format(suffix)))

    if test_run:
        sys.exit("Script executed in 'test_run' mode. FastDesign will not be performed.")

    """
    Running FastDesign
    """
    print("Starting FastDesign")
    redesigns = {}
    for n in range(nstruct):
        start = time.time()
        name = os.path.basename(pdbfile).replace('.pdb', '_fdes{}_{}.pdb'.format(suffix, n))
        pose2 = pose.clone()

        # Mutating all designable residues to alanines, if requested
        if ala_design:
            pose2 = mutate_residues_to_ala(pose2, SEL_mutate_residues)

        # Performing FastDesign
        redesigns[n] = run_fastdesign(pose2, scorefxn_cst,
                                      mutate_residues=SEL_mutate_residues,
                                      repack_residues=SEL_repack_residues,
                                      do_not_repack=SEL_do_not_repack, cst=cst,
                                      rifres=polar_rifres)


        # If more constraint weights are requested then FastDesign is performed
        # that many more times with each constraint weight.
        # Only the final pose is saved
        if len(cst_weights) > 1:
            for i, w in enumerate(cst_weights[1:]):
                print("Designing with cst weight = {:.1f}".format(w))
                i += 1

                print(f"cst score after designing wth weights {cst_weights[i-1]}: "
                      f"{calculate_cst_score(redesigns[n], constraintpairs)}")

                if dump_debug_pdb:
                    redesigns[n].dump_pdb(outdir + name.replace('.pdb', '_w{:.0f}.pdb'.format(cst_weights[i-1])))
                tmp_pose = redesigns[n].clone()

                # Recalculating the design/repack/do-not-touch layers
                SEL_mutate_residues, SEL_repack_residues, SEL_do_not_repack, residues = get_layer_selections(tmp_pose, keep_rifres, constraintpairs, heavyatoms)

                obj_w = build_cst_xml_obj(pose, constraintpairs,
                                          score_wts=args.weights, weight=w)
                cst_w = obj_w.get_mover('add_cst')
                scorefxn_cst_w = obj_w.get_score_function('beta')

                redesigns[n] = run_fastdesign(tmp_pose, scorefxn_cst_w,
                                              mutate_residues=SEL_mutate_residues,
                                              repack_residues=SEL_repack_residues,
                                              do_not_repack=SEL_do_not_repack,
                                              rifres=polar_rifres, cst=cst_w)

        # Dumping the final design as PDB
        redesigns[n].dump_pdb(outdir + name)

        ######################################################
        ## Calculating custom scores/metrics for the design ##
        ######################################################
        additional_scores = {}

        # Calculating the score for how well the constraints were respected
        additional_scores['cst_score'] = calculate_cst_score(redesigns[n], constraintpairs)

        # Shooting the ligand to space
        pose3 = separate_protein_and_ligand(redesigns[n])

        # Re-scoring the design
        scorefxn(redesigns[n])

        # Scoring the separated pose
        scorefxn(pose3)

        # Calculating the SASA of the ligand inside the designed protein
        ligand_surf = getSASA(redesigns[n], resno=ligand_seqpos)

        additional_scores['interf_E'] = redesigns[n].scores['total_score'] - pose3.scores['total_score']
        additional_scores['ligand_E'] = redesigns[n].energies().residue_total_energy(ligand_seqpos)
        additional_scores['L_SASA'] = ligand_surf/ligand_sasa
        additional_scores['bb_rmsd'] = pyr.rosetta.core.scoring.bb_rmsd(redesigns[n], pose)

        # Calculating shape complementarity
        lig_sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(redesigns[n].size())
        protein_sel = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
        sc = sc.fresh_instance()
        sc.use_rosetta_radii(True)
        sc.selector1(protein_sel)
        sc.selector2(lig_sel)
        shapecomp = sc.compute(redesigns[n])
        additional_scores['sc'] = shapecomp.sc

        # Reporting the average score of relevant RIF residues
        rifres_scores = []
        for res in keep_rifres:
            rifres_scores.append(redesigns[n].energies().residue_total_energy(res))
        additional_scores['rifres_E_avg'] = np.average(rifres_scores)

        # If requested, performing no-ligand-repack of the newly designed pose
        # This means removing the ligand from the pose and repacking the residues
        # that were allowed to be redesigned and repacked in FastDesign
        # The repacked pose is dumped and rms scores are appended to the scorefile
        if repack:
            # Removing constraints
            if pose3.constraint_set().has_constraints():
                pose3.constraint_set().clear()

            # Running no-ligand-repack
            nlr_pose = no_ligand_repack(pose3, scorefxn,
                                        repack_residues_list=nlr_repack_residues)

            if True in (repack_dump_pdb, dump_debug_pdb):
                nlr_pose.dump_pdb(outdir + name.replace('.pdb', '_nlr.pdb'))

            # Calculating additional scores for the no-ligand-repack task
            additional_scores['nlr_dE'] = nlr_pose.scores['total_score'] - pose3.scores['total_score']

            # Calculating the nlr_rms of the whole protein
            nlr_rms = rmsd_no_super(pose3, nlr_pose, nlr_repack_residues)
            nlr_rms = pyr.rosetta.core.scoring.all_atom_rmsd(pose3, nlr_pose)
            additional_scores['nlr_totrms'] = nlr_rms

            # Calculating the nlr-rms of each relevant RIF residue (in the 'keep_rifres' list)
            additional_scores['nlr_rms_rifres'] = ""
            rifres_rms = {}
            for m in keep_rifres:
                rifres_rms[m] = rmsd_no_super(pose3, nlr_pose, [m])
                if additional_scores['nlr_rms_rifres'] != "":
                    additional_scores['nlr_rms_rifres'] += "+"
                additional_scores['nlr_rms_rifres'] += "{:.2f}".format(rifres_rms[m])

            # Reporting the average nlr_rms of relevant RIF residues
            additional_scores['nlr_rms_rif_avg'] = np.average([val for key, val in rifres_rms.items()])

        # Reporting how many polar RIF residues are in the design
        if len(polar_rifres) != 0:
            additional_scores['N_polar'] = len(polar_rifres)
        else:
            additional_scores['N_polar'] = 0

        # Reporting relevant RIF residue numbers
        additional_scores['rifres'] = "+".join([str(x) for x in keep_rifres])

        end = time.time()
        print("{}: Iter {} took {:.2f} seconds. "
              "Score = {}".format(pdbfile.rstrip('.pdb'), n, end - start,
                                  redesigns[n].scores['total_score']))

        for k in additional_scores:
            if not isinstance(additional_scores[k], str):
                print(f"{k:<16}: {additional_scores[k]:>7.3f}")
            else:
                print(f"{k:<16}: {additional_scores[k]}")

        # Starting a scorefile if it does not exist yet
        # Writing column headers to the scorefile
        scorefilename = outdir + os.path.basename(pdbfile).replace('.pdb', f'{suffix}.sc')
        print(f"Writing to scorefile: {scorefilename}")

        if not os.path.exists(scorefilename):
            with open(scorefilename, 'w') as file:
                file.write(print_scorefile_header(redesigns[n], pdbfile, additional_scores))

        # Writing the scores from the latest FastDesign iteration to the scorefile
        # TODO: this bit is still a bit wonky, some lines are not written properly
        # when multiple jobs write to the same file.
        with open(scorefilename, 'a') as file:
            file.write(print_scorefile_line(redesigns[n], name, additional_scores))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--pdb", type=str, required=True, help="Input PDB file") #
    parser.add_argument("--nstruct", type=int, default=1, help="How many design iterations?") #
    parser.add_argument("--suffix", type=str, help="Suffix to be added to the end the output filename") #
    parser.add_argument("--outdir", type=str, default="./", help="Directory where the results should be saved to") #
    parser.add_argument("--weights", type=str, default="beta_nov16", help="path to an optional scorefunction weights file") #
    parser.add_argument("--params", type=str, nargs="+", help="params file(s) of ligands and non-canonical residues in the input PDB") #
    parser.add_argument("--dalphaball", type=str, default=f"{SCRIPT_DIR}/DAlphaBall.gcc" help="Path to compiled DAlphaBall.gcc file. By default it uses the file provided in the repository. If it throws erros then download and compile your own DAlphaBall from https://github.com/outpace-bio/DAlphaBall.git")

    parser.add_argument("--min_polar", type=int, default=1, help="minimum number of polar residues contacting the ligand required (default = 1)")  #
    parser.add_argument("--repack", action="store_true", default=False, help="no-ligand-repack will be performed on the FastDesign output & no-ligand-repack RMS values will be appended to the scorefile") #
    parser.add_argument("--repack_dump_pdb", action="store_true", default=False, help="no-ligand-repack output will also be dumped as PDB") #
    parser.add_argument("--ramp_cst_weights", type=float, nargs="+", default=[10.0], help="Constraint weight ramping schedule")
    parser.add_argument("--no_ala_design", action="store_true", default=False, help="disables the step where all designable residues are mutated to ALA before FastDesign") #

    parser.add_argument("--dump_intermediate_pdb", action="store_true", default=False, help=" dumps PDB files at each step. Used for debugging.")  #
    parser.add_argument("--test", action="store_true", default=False, help="does a test run of setting up layers. Will not perform FastDesign")

    args = parser.parse_args()

    main(args)



