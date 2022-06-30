import AFM
from __future__ import print_function
from AFM_Header import *

import ctypes

clib_dir="/home/asi/Desktop/agk/AFM"
import OtherTools as TOOLS

k_inc=1.0/1.0 #[1/A]

# making ellipsoid from gmm
I = "/home/asi/Desktop/test2/rbp/1urp_a.gmm"
R = "/home/asi/Desktop/test2/rbp/2dri_a.gmm"
gm = AFM.GMData(fromFile=True, filename = I)
gm2 = AFM.GMData(fromFile=True, filename = R)
gm.init_ctypes_of_gmm()
gm2.init_ctypes_of_gmm()

#rlist,thetalist,coslist,sinlist
ntheta=360
#rlist=np.array([0.01,0.015,0.02,0.025,0.03,0.035,0.04])
rlist = np.arange(0.02,0.04,0.005)
dtheta=2*np.pi/ntheta
#theta_list=np.linspace(0,2*np.pi,2*np.pi/dtheta,endpoint=False) #y>0 is sufficient IF Centrosymmetry
theta_list=np.arange(0.0,2*np.pi,dtheta)
#print("the in-plane angluar values over which CC is calculated")
#print(theta_list/const.degree)
coslist,sinlist=np.cos(theta_list),np.sin(theta_list)

# the f1 values (observed data) stored
f1list=np.zeros((len(rlist),len(theta_list)))
for i_r in range(len(rlist)):
  for each in range(len(theta_list)):
    kx,ky,kz=TOOLS.get_coord_on_ewald_sphere_from_zero_curvature_coord(rlist[i_r]*coslist[each],rlist[i_r]*sinlist[each],k_inc)
    f1list[i_r][each]=gm.call_get_abs_struct_factor(kx,ky,kz)
f1list=np.log10(f1list)*2.0


# initial:f2
f2list=np.zeros((len(rlist),len(theta_list))) #container of diffraction data
for i_r in range(len(rlist)):
  for each in range(len(theta_list)):
    kx,ky,kz=TOOLS.get_coord_on_ewald_sphere_from_zero_curvature_coord(rlist[i_r]*coslist[each],rlist[i_r]*sinlist[each],k_inc)
    f2list[i_r][each]=gm2.call_get_abs_struct_factor(kx,ky,kz)

f2list=np.log10(f2list)*2.0


# making ellipsoid from gmm
R = "/home/asi/Desktop/test2/rbp/2dri_a.gmm"
gm2 = AFM.GMData(fromFile=True, filename = R)
gm2.init_ctypes_of_gmm()
mat2 = gm2.call_get_abs_struct_factor_all_mesh_on_Ewald_sphere(200,det_width,det_dist,wave_length)
referenceMat= np.log10(mat2)
