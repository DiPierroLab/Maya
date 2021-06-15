This folder contains three proteins and their respective 3mer and 9mer fragment files. I used this to test the basicfolding.py script after making small adjustments. Because it's smaller, it takes less time to run which is helpful when editing. 


The output of running the script goes into the output_files directory. The pdb files show the best structure. There's also a subdirectory called rmsd_files, which contains the root mean square deviation (RMSD) of the best structure aligned to the native structure. 

*Note about the rmsd_files subdirectory: It doesn't work yet!*

*Originally, I aligned the best structure to the native structure manually. Meaning, I aligned the structures myself using PyMOL. First, I transferred the output_files directory from the cluster to my local laptop. Then, I downloaded all the pdbs from the PDB website. These known protein structures are called native structures. Next, I opened each output pdb and aligned it to the native structure. It took time because I had to align proteins one by one. However, the last section of the script is written to automate this process. I use BioPython's retrieve_pdb_file() to get the native structure from the PDB website and PyRosetta's object align_and_get_rmsds() to align the structures. The rmsds should be written for each protein in one .txt file inside the subdirectory. A new .txt file is written each time you run the script, NOT for each individual protein.*

*The problem: The issue I ran into stems from the retrieve_pdb_file() object. For some reason, the retrieved files have the extension .ent, which is not an acceptable argument for pose_from_pdb() and pose_from_pdb() is necessary for align_and_get_rmsds(), i.e. for aligning the structures.*
