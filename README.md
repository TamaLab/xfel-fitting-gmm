 # LICENSING

The program is licensed under the Apache License, Version 2.0. See LICENSE for the full license text.

# AUTHORS of the program

1. Han Asi <mglasgan@gmail.com> - pipeline program, modification of inp.yaml, modification of Monte-Carlo sampling

2. Bhaskar Dasgupta <bhaskardg08@gmail.com> - inp.yaml, Monte-Carlo sampling

3. Tetsuro Nagai <tnagai@fukuoka-u.ac.jp> - tools dealing with the XFEL difraction images of GMM

4. Osamu Miyashita <osamu.miyashita@riken.jp>

5. Florence Tama

   This XFEL fitting program has been developed at Nagoya University and RIKEN.

# HOW TO INSTALL

Run install.sh in folder 'project'

# HOW TO RUN

- Before running the code, user need to make a output directory named s1. (refer to 'input.yaml'> 'OUTPUT'>'IMAGEPATH','STATPATH','GMDATPATH')

- Command to run: 

  python pipeline.py inp.yaml Beg_3.gmm End_814.gmm

- Initial parameters, such as RESTRAINTS, det_width,det_dist, wave_length, SEED, lowest_k and circle_width are included in inp.yaml.

- Output setting is in inp.yaml.

# NOTES

Tested with python2.7, anaconda environment only. May not work correctly with python3.









