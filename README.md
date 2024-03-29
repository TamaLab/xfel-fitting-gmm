 # XFEL-Fitting-GMM

Asi H, Dasgupta B, Nagai T, Miyashita O and Tama F (2022) A hybrid approach to study large conformational transitions of biomolecules from single particle XFEL diffraction data. Front. Mol. Biosci. 9:913860. doi: 10.3389/fmolb.2022.913860


## AUTHORS of the program

1. Bhaskar Dasgupta <bhaskardg08@gmail.com> - inp.yaml, Monte-Carlo sampling

2. Tetsuro Nagai <tnagai@fukuoka-u.ac.jp> - tools dealing with the XFEL difraction images of GMM

3. Han Asi <mglasgan@gmail.com> - pipeline program, modification of inp.yaml, modification of Monte-Carlo sampling

4. Yuki Mochizuki - initial study on XFEL diffraction pattern calculation from GMM

5. Osamu Miyashita <osamu.miyashita@riken.jp>

6. Florence Tama <florence.tama@nagoya-u.jp>

This XFEL fitting program has been developed at Nagoya University and RIKEN.

## RELATED PUBLICATIONS

- Dasgupta, B., Miyashita, O., Uchihashi, T. & Tama, F. Reconstruction of Three-Dimensional Conformations of Bacterial ClpB from High-Speed Atomic-Force-Microscopy Images. Front Mol Biosci **8**, 704274 (2021).
- Dasgupta, B., Miyashita, O. & Tama, F. Reconstruction of low-resolution molecular structures from simulated atomic force microscopy images. Biochim Biophys Acta Gen Subj **1864**, 129420 (2020).
- Nagai, T., Mochizuki, Y., Joti, Y., Tama, F. & Miyashita, O. Gaussian mixture model for coarse-grained modeling from XFEL. Optics Express **26**, 26734 (2018).

## HOW TO INSTALL

In the folder 'project':
```
conda env create -f environment.yml
bash install.sh 
```
## HOW TO RUN

- Before running the code, user need to make an output directory named 's1'. (refer to 'inp.yaml'>'OUTPUT'>'IMAGEPATH', 'STATPATH', 'GMDATPATH')

- [Beg_AK_3.gmm](pipeline/Beg_AK_3.gmm) and [End_AK_214.gmm](pipeline/End_AK_214.gmm) are input gmm files of initial and target conformation. They are created from pdb files using following commands: 
```  
  afmEmulator -f emulate_Beg_AK.yaml
```  
  ([afmEmulator](project/bin/afmEmulator) is in /project/bin, [emulate_Beg_AK.yaml](pipeline/input_files/AK/emulate_Beg_AK.yaml) is located at /pipeline/input_files/AK)
```  
  gmconvert A2G -ipdb End_AK.pdb -ogmm End_AK_214.gmm -ng 214
```  
  (gmconvert software: https://pdbj.org/gmfit/doc_gmconvert/README_gmconvert.html)

- Command to run XFEL GMM fitting: 
```
  python pipeline.py inp.yaml Beg_AK_3.gmm End_AK_214.gmm
```  
- Initial parameters, such as RESTRAINTS, det_width,det_dist, wave_length, SEED, lowest_k, and circle_width are included in [inp.yaml](pipeline/inp.yaml).

- Output setting is in [inp.yaml](pipeline/inp.yaml).

## NOTES

Tested with python 2.7, anaconda environment only. May not work correctly with python 3.

## LICENSING

The program is licensed under the Apache License, Version 2.0. See LICENSE for the full license text.







