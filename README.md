# Read Me
The goal of my script basicfolding.py is to start with an amino acid sequence and end with a folded protein structure. I used PyRosetta software to write the code and PyMOL to visualize the folded proteins. The script is based on the [tutorial from PyRosetta's github](https://rosettacommons.github.io/PyRosetta.notebooks/). While the majority of this repo is about the script itself, this version of basicfolding.py was written to be used on the Discovery cluster (information about the cluster can be found on Northeastern's [Research Computing website](https://rc-docs.northeastern.edu/en/latest/index.html#)). Because of that, I included a section at the end of this README about the sbatch command and shell script I wrote to submit basicfolding.py as a job to the cluster. 

You'll notice the last chunk of code is commented out. This is because I'm trying to use BioPython's retrieve_pdb_file() to get the native structure of a protein from the PDB website. But retrieve_pdb_file() gets me a .ent file, which is not acceptable argument since the following objects only take pdb files. There's more about this issue in **The Details** section below.

*Disclaimer: I am not a coder, so any advice/suggestions would be greatly appreciated!* 

### Table of Contents 
* **The Basics** - definitons of important terms/words I use throughout the repo for explanations
* **Random notes** - little things I made note of during the process of developing the script and running it on the cluster
* **The Details** - deeper explanations about parts of the script itself or issues I ran into while writing it
* **The Cluster** - using sbatch and writing a shell script for Discovery 

## The Basics
**input**: fasta (.fasta) files and fragment (.txt) files
* fasta file: protein-ID.fasta
* 3mer file: protein-ID_aat000_03_05.200_v1_3.txt
* 9mer file: protein-ID_aat000_09_05.200_v1_3.txt

**output**: pdb (.pdb) files and rmsd (.txt) files
* pdb file: protein-ID_out_2021_%m_%d-%H:%M_%p.pdb
* rmsd file: rmsd__2021_%m_%d-%H:%M_%p.txt

**protein-ID**: 4-digit protein id from PDB website

**best structure/output structure/out structure**: structure created by basicfolding.py

**native structure**: known structure as see on PDB website

**Pose**: a PyRosetta object
  * represents a single molecule and is constructed from a PDB file
  * contains various types of info that describe the structure including core components Energies, PDBInfo and Conformation

## Random notes
* Make sure the protein IDs are written in all CAPS (otherwise the cluster won't be able to match the fragment files to the right fasta files) 
* Make sure the directory has all the fasta files and fragment files (fragment files can be generated using [robetta](http://robetta.bakerlab.org/fragmentsubmit.jsp)) 
* The output structures should be in an output directory that has a subdirectory for the RMSDs
* The in_file and cen_file files that are created are not the final structures that should be analyzed. I made them because I wanted to see the initial input structure, which is essentially a squiggly line, and the centroid model, which looks the same as the final output structure. 

## The Details
About the rmsd_files subdirectory in the output: It doesn't work yet!
* Originally, I aligned the best structure to the native structure manually. Meaning, I aligned the structures myself using PyMOL. First, I transferred the output_files directory from the cluster to my local laptop. Then, I downloaded all the pdbs from the PDB website. These known protein structures are called native structures. Next, I opened each output pdb and aligned it to the native structure. It took time because I had to align proteins one by one. However, the last section of the script is written to automate this process. I use BioPython's retrieve_pdb_file() to get the native structure from the PDB website and PyRosetta's object align_and_get_rmsds() to align the structures. The rmsds should be written for each protein in one .txt file inside the subdirectory. A new .txt file is written each time you run the script, NOT for each individual protein.
* The problem: The issue I ran into stems from the retrieve_pdb_file() object. For some reason, the retrieved files have the extension .ent, which is not an acceptable argument for pose_from_pdb() and pose_from_pdb() is necessary for align_and_get_rmsds(), i.e. for aligning the structures.

Using fragments instead of angles in Subroutine 1
* In the tutorial, protein angles (phi, psi, omega) are randomly perturbed and the energy score is calculated based on these varying angles. This creates fine structures but after doing some refining research, I realized I could improve my structures using homology modeling. Homology modeling (aka "comparative modeling") is the process of creating a protein structure based on its amino acid sequence and an experimental structure of a related homologous protein (meaing the two proteins have a common ancestor). Since I started with amino acid sequences of proteins that already have experiments structures, I didn't need to worry about the percentage of sequence identity. Instead, I started with a fasta file and put the sequence into [robetta](http://robetta.bakerlab.org/fragmentsubmit.jsp) to get the 3mer and 9mer fragment files. Then, I used the fragment files in the script to randomly perturb fragments of the protein, rather than angles. This resulted in much better structures with lower RMSDs.

In Subroutine 3, the decision is made using the Metropolis criterion
* From the PyRosetta github: 
  "The Metropolis criterion has a probability of accepting a move as P=exp(−ΔG/kT). When ΔE≥0, the Metropolis criterion probability of accepting the move is P=exp(−ΔG/kT). When ΔE<0, the Metropolis criterion probability of accepting the move is P=1. Use kT=1 Rosetta Energy Unit (REU)."
  
## The Cluster
This is my shell script:
```
#!/bin/bash
#SBATCH --job-name=basicfolding
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=folding.%N.%j.out
#SBATCH --error=folding.%N.%j.err
source /home/<username>/anaconda3/bin/activate <virtualenvironment>
python basicfolding.py
```

Once I wrote this in a text (.txt) file, I saved it as folding.sh. Then, I copied it to the cluster and used the command
```
sbatch folding.sh
```
to run the job in the background. The job took about ~hours for 15 proteins. See the RC website information about [sbatch](https://rc-docs.northeastern.edu/en/latest/using-discovery/sbatch.html) and [how to check on your job](https://rc-docs.northeastern.edu/en/latest/using-discovery/usingslurm.html) while it's queued or running.
