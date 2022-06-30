 # LICENSING

The program is licensed under the Apache License, Version 2.0. See LICENSE for the full license text.

# AUTHORS of the program

1. Bhaskar Dasgupta <bhaskardg08@gmail.com> - inp.yaml, Monte-Carlo sampling

2. Tetsuro Nagai <tnagai@fukuoka-u.ac.jp> - tools dealing with the XFEL difraction images of GMM

3. Han Asi <mglasgan@gmail.com> - pipeline program, modification of inp.yaml, modification of Monte-Carlo sampling

4. Yuki Mochizuki - initial study on XFEL diffraction pattern calculation from GMM

5. Osamu Miyashita <osamu.miyashita@riken.jp>

6. Florence Tama <florence.tama@nagoya-u.jp>



   This XFEL fitting program has been developed at Nagoya University and RIKEN.

# HOW TO INSTALL

Run install.sh in folder 'project'

# HOW TO RUN

- Before running the code, user need to make a output directory named s1. (refer to 'input.yaml'> 'OUTPUT'>'IMAGEPATH','STATPATH','GMDATPATH')

- Command to run: 

  python pipeline.py inp.yaml Beg_EF2.gmm End_EF2_814.gmm
  
- Beg_EF2.gmm and End_EF2_814.gmm are input gmm file of initial and target conformation. They are created from pdb files using following commands: 
  
  afmEmulator -f emulate_Beg_EF2.yaml
  
  (afmEmulator is in project/bin, emulate_Beg_EF2.yaml is located at /pipeline/input_files)
  
  gmconvert A2G -ipdb End_EF2.pdb -ogmm End_EF2_814.gmm -ng 819
  
  (gmconvert software: https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html)

- Initial parameters, such as RESTRAINTS, det_width,det_dist, wave_length, SEED, lowest_k and circle_width are included in inp.yaml.

- Output setting is in inp.yaml.

# NOTES

Tested with python2.7, anaconda environment only. May not work correctly with python3.









