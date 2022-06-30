 # XFEL-Fitting-GMM
 
 A Hybrid Approach to Study Large Conformational Transitions of Biomolecules from Single Particle XFEL Diffraction Data
 
 Asi, H., Dasgupta, B., Nagai, T., Miyashita, O. & Tama, F.

*submitted*

# AUTHORS of the program

1. Bhaskar Dasgupta <bhaskardg08@gmail.com> - inp.yaml, Monte-Carlo sampling

2. Tetsuro Nagai <tnagai@fukuoka-u.ac.jp> - tools dealing with the XFEL difraction images of GMM

3. Han Asi <mglasgan@gmail.com> - pipeline program, modification of inp.yaml, modification of Monte-Carlo sampling

4. Yuki Mochizuki - initial study on XFEL diffraction pattern calculation from GMM

5. Osamu Miyashita <osamu.miyashita@riken.jp>

6. Florence Tama <florence.tama@nagoya-u.jp>

   This XFEL fitting program has been developed at Nagoya University and RIKEN.

# RELATED PUBLICATIONS

- Dasgupta, B., Miyashita, O., Uchihashi, T. & Tama, F. Reconstruction of Three-Dimensional Conformations of Bacterial ClpB from High-Speed Atomic-Force-Microscopy Images. Front Mol Biosci 8, 704274 (2021).
- Dasgupta, B., Miyashita, O. & Tama, F. Reconstruction of low-resolution molecular structures from simulated atomic force microscopy images. Biochim Biophys Acta Gen Subj 1864, 129420 (2020).
- Nagai, T., Mochizuki, Y., Joti, Y., Tama, F. & Miyashita, O. Gaussian mixture model for coarse-grained modeling from XFEL. Optics Express 26, 26734 (2018).

# HOW TO INSTALL

In folder 'project':
```
conda env create -f environment.yml

bash install.sh 
```
# HOW TO RUN

- Before running the code, user need to make a output directory named s1. (refer to 'input.yaml'> 'OUTPUT'>'IMAGEPATH','STATPATH','GMDATPATH')

- Beg_EF2_7.gmm and End_EF2_814.gmm are input gmm file of initial and target conformation. They are created from pdb files using following commands: 
```  
  afmEmulator -f emulate_Beg_EF2.yaml
```  
  (afmEmulator is in project/bin, emulate_Beg_EF2.yaml is located at /pipeline/input_files/EF2)
```  
  gmconvert A2G -ipdb End_EF2.pdb -ogmm End_EF2_814.gmm -ng 819
```  
  (gmconvert software: https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html)- Command to run XFEL GMM fitting: 
```
  python pipeline.py inp.yaml Beg_EF2_7.gmm End_EF2_814.gmm
```  
- Initial parameters, such as RESTRAINTS, det_width,det_dist, wave_length, SEED, lowest_k and circle_width are included in inp.yaml.

- Output setting is in inp.yaml.

# NOTES

Tested with python2.7, anaconda environment only. May not work correctly with python3.

# LICENSING

The program is licensed under the Apache License, Version 2.0. See LICENSE for the full license text.







