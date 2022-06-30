#########################
# this script is a part of GMM_DIFRACT_TOOLS. 
# Written by T. Nagai 
# tested with python3 only; may not work well with python2
# Here some tools are defined
#########################
import numpy as np
import sys
import ctypes
import scipy.constants as const
import scipy
from  scipy import interpolate
from  scipy import stats

clib_dir="/home/asi/Desktop/GMM_DIFRACT_TOOLS-master/lib/"


import numba
@numba.jit
def get_overlap_with_jit(f1list,f2list):
        H=0.0
        for each in range(len(f1list)):
            H+=get_overlap_with_jit_subrotine(f1list[each],f2list[each])
        return H/len(f1list)

@numba.jit
def get_overlap_with_jit_subrotine(data1,data2):
    ave1=np.mean(data1)
    ave2=np.mean(data2)

    std1=np.std(data1)
    std2=np.std(data2)

    cross_term=np.dot(data1,data2)
    cross_term/=len(data1)
    cross_term-=ave1*ave2
    
    return cross_term/(std1*std2)


def get_overlap_with_scipy(f1list,f2list):
    """ this function is equivalent to get_overlap_with_jit and might be faster only when f1 has very very large entries. 
    However, for most of cases, the jit implementation should be faster """
    H=0
    for each in range(len(f1list)):
        H+=scipy.stats.pearsonr(f1list[each], f2list[each])[0]
    return H/len(f1list)

# this is not used anymore. 
#def read_spsim_data_and_return_interpolation(file="./001_log_pt.txt"):
#    """tentatively data must be 200x200 pixcel covering [-0.1,0.1]**2 """
#    data1=np.loadtxt(file) #200x200
#    data1aligned=np.rot90(data1.T).T
#    grid_x200=np.linspace(-0.1,0.1,200)
#    grid_y200=np.linspace(-0.1,0.1,200)
#    f1 = interpolate.interp2d(grid_x200, grid_y200, data1aligned.T, kind='cubic') #interpolate takes transposed coordinate....
#    return f1

def read_gnuplot_optimized_spsim_data_and_return_interpolation(file="./001_log_pt.txt"):
    """tentatively data is must be 200x200 pixcel covering [-0.1,0.1]**2 """
    data1=np.loadtxt(file) #200x200
    data1aligned=data1.T
    grid_x200=np.linspace(-0.1,0.1,200)
    grid_y200=np.linspace(-0.1,0.1,200)
    f1 = interpolate.interp2d(grid_x200, grid_y200, data1aligned.T, kind='cubic') #interpolate takes transposed coordinate....
    return f1

def convert_Euler_angles_to_matrix(gamma,beta,alpha):
    R1=R_Rodrigues(0,0,1,gamma)
    R2=R_Rodrigues(0,1,0,beta)
    R3=R_Rodrigues(0,0,1,alpha)
    return np.dot(R1,R2,R3)

#Rodrigues's rotation matrix
def R_Rodrigues(x,y,z,omega):
    """omega must be in the unit of radian, and (x,y,z) should be unit vector, although this will be checked and normalized"""
    x2,y2,z2=x*x,y*y,z*z
    r2=x2+y2+z2
    if np.abs(r2-1.0) >1e-10:
        sys.stderr.write('The vector will be normalized; Check the unity of unit vector  %f\n'%r2)
        r=np.sqrt(r2)
        x,y,z=x/r,y/r,z/r
    C,S=np.cos(omega), np.sin(omega)
    matrix = np.array([[C+x2*(1-C),     x*y*(1-C)-z*S, x*z*(1-C)+y*S],
                       [y*x*(1-C)+z*S,  C+y2*(1-C),    y*z*(1-C)-x*S ],
                       [z*x*(1-C)-y*S,  z*y*(1-C)+x*S, C+z2*(1-C)]])
    return matrix

def get_coord_on_ewald_sphere_from_zero_curvature_coord(kx_prime,ky_prime,k_inc,isfromZplus=True):
    chi=np.sqrt(k_inc**2+kx_prime**2+ky_prime**2)
    kx=k_inc*kx_prime/chi
    ky=k_inc*ky_prime/chi
    kz=k_inc-k_inc**2/chi
    if not isfromZplus:
        kz*=-1.0
    return kx,ky,kz

#def detector_coord_to_coord_on_ewald_sphere(dx,dy,det_dist,k_inc,isfromZplus=True):
def get_coord_on_ewald_shere_from_detector_coord(dx,dy,det_dist,k_inc,isfromZplus=True):
    kx=k_inc*dx/np.sqrt(det_dist*det_dist+dx*dx+dy*dy)
    ky=k_inc*dy/np.sqrt(det_dist*det_dist+dx*dx+dy*dy)
    kz=k_inc-np.sqrt(k_inc*k_inc-kx*kx-ky*ky)
    if not isfromZplus:
        kz*=-1.0
    return kx,ky,kz
