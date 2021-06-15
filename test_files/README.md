This folder contains the fasta files for 8 randomly selected proteins and the corresponding fragment files. I downloaded the fasta files from the PDB website https://www.rcsb.org/ and used Robetta http://robetta.bakerlab.org/fragmentsubmit.jsp to get the fragment files. 


The output of running the script goes into the output_files directory. The pdb files show the best structure. There's also a subdirectory called rmsd_files, which contains the root mean square deviation (RMSD) of the best structure aligned to the native structure.


*Note: Though the proteins were randomly selected, they all only have Chain A, as stated in the fasta files. (Example: Protein 6Q21 has four chains, Chains A, B, C and D.)*

*Note: In my experiment, I used 15 proteins, including 6Q21 and the script works the same. Unfortunately, I couldn't upload the other 7 fragment files because they were too big.*
