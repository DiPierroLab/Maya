#!/usr/bin/env python
# coding: utf-8
import optparse
import math
import random
from datetime import datetime
import pyrosetta.rosetta as rsa
import os, glob
from Bio import SeqIO
from Bio import PDBList
from pyrosetta import *
init()

#### Fill in your directory 
#### ****Change iterations!

directory = '/path/to/test/files/'
scratchdir = '/scratch/username/'
date = datetime.now().strftime('%Y_%m_%d-%H:%M_%p')

fastafiles = glob.glob(str(directory) + '*.fasta') # read all .fasta files

for fastafile in fastafiles:
    for record in SeqIO.parse(fastafile, 'fasta'): # I included this part because the record.id would print 
        protein_id = record.id                     # more than necessary (i.e. ">6Q21_1|Chains" instead of 
        protein = protein_id.split('_')[0]         # "6Q21".)
        sequence = str(record.seq)
        
    in_file = '%s%s_in_%s.pdb' %(scratchdir, protein, date) # *Not necessary*
    cen_file = '%s%s_cen_%s.pdb' %(scratchdir, protein, date) #Not necessary*
    out_file = '%soutput_files/%s_out_%s.pdb' %(directory, protein, date) # Output structures will go into
                                                                          # output_files folder
    # Creating initial pose and structure from sequence 
    pose = pose_from_sequence(sequence, 'fa_standard')
    pose.pdb_info().name('%s' %protein)
    dump_pdb(pose, in_file)
    input_pose = pose_from_pdb(in_file)
    
    for i in range(1, pose.total_residue() + 1):
        pose.set_phi(i, -180)
        pose.set_psi(i, 180)
        pose.set_omega(i, 180)
    
    # Switching to centroid model (allows cleaner perturbations) 
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    to_centroid.apply(pose)
    
    # Sending structure to PyMOL for visualization 
    pmm = PyMOLMover()
    pmm.keep_history(True)
    pmm.apply(pose)
    movemap = MoveMap()
    movemap.set_bb(True)

    # Loading fragment information (*NOTE: These .txt files must be downloaded from robetta.org
    #                               I used the old website http://robetta.bakerlab.org/fragmentsubmit.jsp)
    long_frag_file = '%s%s_aat000_09_05.200_v1_3.txt' %(directory, protein)
    long_frag_length = 9
    long_inserts = 1
    short_frag_file = '%s%s_aat000_03_05.200_v1_3.txt' %(directory, protein)
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
    
    # Subroutine 1|Random move: Randomly perturb fragments to make new structure 
    def randTrial(pose):    
        folding_mover.add_mover(insert_long_frag)
        folding_mover.add_mover(insert_short_frag)
        long_frag_mover.apply(pose)
        short_frag_mover.apply(pose)
        pmm.apply(pose)
        return pose

    # Subroutine 2|Scoring: Returns the numerical energy score of pose
    scorefxn_low = create_score_function('score3')
    scorefxn_high = get_fa_scorefxn()
    def score(pose):
        return scorefxn_low(pose)

    # Subroutine 3|Decision: Accepting/rejecting the new (perturbed) structure based on Metropolis criterion
    def decision(before_pose, after_pose):
        E = score(after_pose) - score(before_pose) # E is the change in total energy as the structure changes
        if E < 0:                                  # If E is less then zero, keep the new structure
            return after_pose
        elif random.uniform(0, 1) >= math.exp(-E/1):
            return before_pose
        else:
            return after_pose

    # Execution 
    def basic_folding(pose):
        lowest_pose = Pose()    
        for i in range(100): # number of times the 3 subroutines are performed 
            if i == 0:
                lowest_pose.assign(pose)
            before_pose = Pose()
            before_pose.assign(pose)
            after_pose = Pose()
            after_pose.assign(randTrial(pose))
           
            pose.assign(decision(before_pose, after_pose)) # Only keep lowest energy conformation that is
            if score(pose) < score(lowest_pose):           # acheived at any point during the simulation 
                lowest_pose.assign(pose)    
        return lowest_pose
    
    cen = basic_folding(pose)

    # Switching back to full atom model 
    dump_pdb(cen, cen_file)
    final_pose = pose_from_pdb(cen_file)
    to_fullatom.apply(final_pose)
    
    # Refining the structure 
    fast_relax_mover = rsa.protocols.relax.FastRelax(scorefxn_high)
    fast_relax_mover.apply(final_pose)
    
    # Creating output pose and .pdb file
    dump_pdb(final_pose, out_file)
    output_pose = pose_from_pdb(out_file)
    
### This part is a work in progress. I want to use retrieve_pdb_file() to get the PDBs for each protein 
##. from the PDB website. Even though I specify the file format, the files it downloads are .ent files,
##. which aren't compatible arguments with the following pdb_from_pdb() object. Any suggestions/advice
##. is greatly appreciated! 

#     pdbl = PDBList()
#     native_pdb = pdbl.retrieve_pdb_file('%s' %protein, pdir='/stratch/username/', file_format='pdb')
#     native_pose = pose_from_pdb(native_pdb)

#     def align_and_get_rmsds(native_pose, output_pose): # final step: aligning structures and calculating RMSD
#         pyrosetta.rosetta.core.pose.full_model_info.make_sure_full_model_info_is_setup(native_pose)
#         rmsds = []
#         pyrosetta.rosetta.core.pose.full_model_info.make_sure_full_model_info_is_setup(output_pose)
#         rmsds += [pyrosetta.rosetta.protocols.stepwise.modeler.align.superimpose_with_stepwise_aligner(native_pose, output_pose)]
#         return rmsds
    
#     rmsd = align_and_get_rmsds(native_pose, output_pose)
#     with open('%soutput_files/rmsd_files/rmsd_%s.txt' %(directory, date), 'a') as writer: # writing the RMSDs to a .txt file
#         writer.write('\nThe rmsd of %s is %s' %(protein, rmsd))

print(fastafiles) # Check .fasta files

