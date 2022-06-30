 # LICENSING

The program is licensed under the Apache License, Version 2.0. See LICENSE for the full license text.

# AUTHORS of the program

1. Han Asi <mglasgan@gmail.com> - pipeline program, modification of inp.yaml, modification of Monte-Carlo sampling

2. Bhaskar Dasgupta <bhaskardg08@gmail.com> - inp.yaml, Monte-Carlo sampling

3. Tetsuro Nagai <tnagai@fukuoka-u.ac.jp> - tools dealing with the XFEL difraction images of GMM

4. Osamu Miyashita <osamu.miyashita@riken.jp>

5. Florence Tama <florence.tama@nagoya-u.jp>



   This XFEL fitting program has been developed at Nagoya University and RIKEN.

# HOW TO INSTALL

Run install.sh in folder 'project'

# HOW TO RUN

- Before running the code, user need to make a output directory named s1. (refer to 'input.yaml'> 'OUTPUT'>'IMAGEPATH','STATPATH','GMDATPATH')

- Command to run: 

  python pipeline.py inp.yaml Beg_3.gmm End_814.gmm
  
- Beg_3.gmm and End_814.gmm are input gmm file of initial and target conformation. They are created from pdb files using following commands: 
  
  gmconvert A2G -ipdb Beg.pdb -ogmm Beg.gmm -ng 3
  
  gmconvert A2G -ipdb End.pdb -ogmm End.gmm -ng 819
  
  (gmconvert software: https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html)

- Initial parameters, such as RESTRAINTS, det_width,det_dist, wave_length, SEED, lowest_k and circle_width are included in inp.yaml.

- Output setting is in inp.yaml.

# NOTES

Tested with python2.7, anaconda environment only. May not work correctly with python3.









