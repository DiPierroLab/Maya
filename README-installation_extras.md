# License and Installation
Rosetta Commons is the hub for Rosetta modeling software. Under the "License and Download" tab on their [Software page](https://www.rosettacommons.org/software/license-and-download), you can find all the information and links for each software package (Rosetta Software Suite, PyRosetta and Foldit Standalone). I used [PyRosetta](http://www.pyrosetta.org/home) and all the following links are specific to my process obtaining a license and downloading the PyRosetta 4 version. But if you want another software package, start on the [Software page](https://www.rosettacommons.org/software/license-and-download) and use links that correspond with the software you want.


Before installion, you have to get a [license](https://els2.comotion.uw.edu/product/pyrosetta) in order to download the PyRosetta package. After obtaining your license, installation instructions, information about different builds and download links are available on the [PyRosetta Downloads page](http://www.pyrosetta.org/dow). I used the Mac OS X 10.7 "Lion/MountainLion" (64-bit), Python-3.7.MinSizeRel version and followed the GNU/Linux and Mac OS X instructions under the section 
"Installation with an environment manager."

The downloaded zip file is called PyRosetta4.MinSizeRel.python37.mac.release-285.tar.bz2, which I renamed PyRosetta4 for simplicity. Unfortunately the file is too large to upload directly in this repo, but below is an outline of what's inside:

- apps
  - delta_score_per_mutation.py
  - ...
  - surface_docking
- demo
  - __init__.py
  - D010_Pose_structure.py
  - ...
  - untested
- PyMOL-Rosetta-relay-client.python2.py
- PyMOL-Rosetta-relay-client.python3.py
- PyMOL-RosettaServer.py
- PyMOL-RosettaServer.python3.py
- pyrosetta (shortcut)
- rosetta (shortcut)
- self-test.py
- setup
  - ez_setup.py
  - ...
  - setup.py
- test
  - __init__.py
  - __pycache__
  - ...
  - Workshop9_my_shapes.py
- version.json

# Extras
#### PyMOL
I used PyMOL to visual the structures. Like PyRosetta, you have to [register](https://pymol.org/edu/?q=educational/) to get the app (it should be free if you're a student and using PyMOL for educational/academic research purposes). The [PyMOL website](https://pymol.org/2/) also provides the link for obtaining a license and downloads for Windows/Mac/Linux and instructions for how to do so. 

#### Other notes
I used Jupyter Notebook to write basicfolding.py. [Installation instructions](https://jupyter.org/install) are on the [Project Jupyter website](https://jupyter.org/). I used conda for installation but you can also [try the different Jupyter servers in a web browser](https://jupyter.org/try) if you don't want to install anything. 
