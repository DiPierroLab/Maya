# Read Me

This version of the code was written to be used on the Discovery cluster. It's based on the tutorial from PyRosetta's github found on https://rosettacommons.github.io/PyRosetta.notebooks/.  
**Disclaimer: I am not a coder, so any advice/suggestions would be greatly appreciated!** 

Input : .fasta files and fragment .txt files  
Output: .pdb files and rmsd.txt files

*Make sure protein IDs are written in all CAPS (otherwise the cluster won't be able to match the fragment files to the right fasta files) 
*Make sure the directory has all the fasta files and fragment files (fragment files can be generated using robetta) 
*The output structures should be in an output directory that has a subdirectory for the RMSDs



The in_file and cen_file files that are created are not necessary for the code. I made them because I wanted to see the initial input structure, which is essentially a squiggly line, and the centroid model, which looks the same as the final output structure. 



About the code:

In Subroutine 3, the decision is made using the Metropolis criterion. From the PyRosetta github: 
  "The Metropolis criterion has a probability of accepting a move as P=exp(−ΔG/kT). When ΔE≥0, the Metropolis criterion probability of accepting the move is P=exp(−ΔG/kT). When ΔE<0, the Metropolis criterion probability of accepting the move is P=1. Use kT=1 Rosetta Energy Unit (REU)."
