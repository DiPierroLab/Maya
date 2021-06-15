# Read Me
This version of the basicfolding.py script was written to be used on the Discovery cluster (information about the cluster can be found on Northeastern's [Research Computing website](https://rc-docs.northeastern.edu/en/latest/index.html#)). The script is based on the [tutorial from PyRosetta's github](https://rosettacommons.github.io/PyRosetta.notebooks/)
*Disclaimer: I am not a coder, so any advice/suggestions would be greatly appreciated!* 

### The Basics
**input**: fasta (.fasta) files and fragment (.txt) files
* fasta file: protein-ID.fasta
* 3mer file: protein-ID_aat000_03_05.200_v1_3.txt
* 9mer file: protein-ID_aat000_09_05.200_v1_3.txt

**output**: pdb (.pdb) files and rmsd (.txt) files
* pdb file: protein-ID_out_2021_%m_%d-%H:%M_%p.pdb
* rmsd file: rmsd__2021_%m_%d-%H:%M_%p.txt

**protein-ID**: 4-digit protein id from PDB website

**best structure/output structure/out structure**: structure created by basicfolding.py

**native structure**: structure as see on PDB website

**Pose**
  * represents a single molecule and is constructed from a PDB file
  * contains various types of info that describe the structure including core components Energies, PDBInfo and Conformation

### Random notes
* Make sure the protein IDs are written in all CAPS (otherwise the cluster won't be able to match the fragment files to the right fasta files) 
* Make sure the directory has all the fasta files and fragment files (fragment files can be generated using [robetta](http://robetta.bakerlab.org/fragmentsubmit.jsp)) 
* The output structures should be in an output directory that has a subdirectory for the RMSDs
* The in_file and cen_file files that are created are not the final structures that should be analyzed. I made them because I wanted to see the initial input structure, which is essentially a squiggly line, and the centroid model, which looks the same as the final output structure. 

### The Details
About the rmsd_files subdirectory in the output: It doesn't work yet!
* Originally, I aligned the best structure to the native structure manually. Meaning, I aligned the structures myself using PyMOL. First, I transferred the output_files directory from the cluster to my local laptop. Then, I downloaded all the pdbs from the PDB website. These known protein structures are called native structures. Next, I opened each output pdb and aligned it to the native structure. It took time because I had to align proteins one by one. However, the last section of the script is written to automate this process. I use BioPython's retrieve_pdb_file() to get the native structure from the PDB website and PyRosetta's object align_and_get_rmsds() to align the structures. The rmsds should be written for each protein in one .txt file inside the subdirectory. A new .txt file is written each time you run the script, NOT for each individual protein.
* The problem: The issue I ran into stems from the retrieve_pdb_file() object. For some reason, the retrieved files have the extension .ent, which is not an acceptable argument for pose_from_pdb() and pose_from_pdb() is necessary for align_and_get_rmsds(), i.e. for aligning the structures.

In Subroutine 3, the decision is made using the Metropolis criterion
* From the PyRosetta github: 
  "The Metropolis criterion has a probability of accepting a move as P=exp(−ΔG/kT). When ΔE≥0, the Metropolis criterion probability of accepting the move is P=exp(−ΔG/kT). When ΔE<0, the Metropolis criterion probability of accepting the move is P=1. Use kT=1 Rosetta Energy Unit (REU)."
