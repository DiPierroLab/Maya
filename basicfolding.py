#!/usr/bin/env python
# coding: utf-8
import optparse
import math
import random
import pyrosetta
import pyrosetta.rosetta as rsa
from pyrosetta.rosetta import core, protocols
import os, glob
import Bio
from Bio import SeqIO
from pyrosetta import *
init()

#### Fill in the directory and date
#### ****Change iterations! 
directory = '/Users/mlane/PyRosetta4/test_mini/'
date = '6.4.21'

fastafiles = glob.glob(str(directory) + '*.fasta')

for fastafile in fastafiles:
    for record in SeqIO.parse(fastafile, 'fasta'):
        protein_id = record.id
        protein = protein_id.split('_')[0]
        sequence = str(record.seq)
        
    in_file = '%s_in_%s.pdb' %(protein, date)
    cen_file = '%s/cen_files/%s_cen_%s.pdb' %(directory, protein, date)
    out_file = '%s/output_files/%s_out_%s.pdb' %(directory, protein, date)
    
    pose = pose_from_sequence(sequence, 'fa_standard')
    pose.pdb_info().name('%s' %protein)
    dump_pdb(pose, in_file)
    input_pose = pose_from_pdb(in_file)
    
    for i in range(1, pose.total_residue() + 1):
        pose.set_phi(i, -180)
        pose.set_psi(i, 180)
        pose.set_omega(i, 180)
    
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    to_centroid.apply(pose)
    
    pmm = PyMOLMover()
    pmm.keep_history(True)
    pmm.apply(pose)
    movemap = MoveMap()
    movemap.set_bb(True)

    long_frag_file = '%s/%s_aat000_09_05.200_v1_3.txt' %(directory, protein)
    long_frag_length = 9
    long_inserts = 1
    short_frag_file = '%s/%s_aat000_03_05.200_v1_3.txt' %(directory, protein)
    short_frag_length = 3
    short_inserts = 3    
    fragset_long = rsa.core.fragment.ConstantLengthFragSet(long_frag_length)
    fragset_long.read_fragment_file(long_frag_file)
    fragset_short = rsa.core.fragment.ConstantLengthFragSet(short_frag_length)
    fragset_short.read_fragment_file(short_frag_file)
    
    long_frag_mover = rsa.protocols.simple_moves.ClassicFragmentMover(fragset_long, movemap)
    short_frag_mover = rsa.protocols.simple_moves.ClassicFragmentMover(fragset_short, movemap)
    pmm.send_movemap(pose, movemap)
    insert_long_frag = rsa.protocols.moves.RepeatMover(long_frag_mover, long_inserts)
    insert_short_frag = rsa.protocols.moves.RepeatMover(short_frag_mover, short_inserts)     
    folding_mover = rsa.protocols.moves.SequenceMover()
    
    def randTrial(pose):    
        folding_mover.add_mover(insert_long_frag)
        folding_mover.add_mover(insert_short_frag)
        long_frag_mover.apply(pose)
        short_frag_mover.apply(pose)
        pmm.apply(pose)
        return pose

    scorefxn_low = create_score_function('score3')
    scorefxn_high = get_fa_scorefxn()
    def score(pose):
        return scorefxn_low(pose)

    def decision(before_pose, after_pose):
        E = score(after_pose) - score(before_pose)
        if E < 0:
            return after_pose
        elif random.uniform(0, 1) >= math.exp(-E/1):
            return before_pose
        else:
            return after_pose

    def basic_folding(pose):
        lowest_pose = Pose()    
        for i in range(3):
            if i == 0:
                lowest_pose.assign(pose)
            before_pose = Pose()
            before_pose.assign(pose)
            after_pose = Pose()
            after_pose.assign(randTrial(pose))
           
            pose.assign(decision(before_pose, after_pose))
            if score(pose) < score(lowest_pose):
                lowest_pose.assign(pose)    
        return lowest_pose
    
    cen = basic_folding(pose)

    dump_pdb(cen, cen_file)
    final_pose = pose_from_pdb(cen_file)
    to_fullatom.apply(final_pose)
    
    fast_relax_mover = rsa.protocols.relax.FastRelax(scorefxn_high)
    fast_relax_mover.apply(final_pose)
    
    dump_pdb(final_pose, out_file)
    output_pose = pose_from_pdb(out_file)

    def align_and_get_rmsds(input_pose, output_pose):
        pyrosetta.rosetta.core.pose.full_model_info.make_sure_full_model_info_is_setup(input_pose)
        rmsds = []
        pyrosetta.rosetta.core.pose.full_model_info.make_sure_full_model_info_is_setup(output_pose)
        rmsds += [pyrosetta.rosetta.protocols.stepwise.modeler.align.superimpose_with_stepwise_aligner(input_pose, output_pose)]
        return rmsds
    
    rmsd = align_and_get_rmsds(input_pose, output_pose)
    with open('%s/rmsd_files/rmsd.txt' %directory, 'a') as writer:
        writer.write('\nThe rmsd of %s is %s' %(protein, rmsd))

print(fastafiles)

