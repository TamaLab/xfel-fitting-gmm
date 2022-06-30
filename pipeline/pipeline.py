import sys
import AFM
import numpy as np
import yaml

inp_yaml_name = sys.argv[1]
with open(inp_yaml_name, 'r') as stream:
  opt = yaml.load(stream)

det_width = float(opt['det_width'])
det_dist  = float(opt['det_dist']) 
wave_length = float(opt['wave_length']) 

I_filename = sys.argv[2]
#I = "Beg_3.gmm"
gm = AFM.GMData(fromFile=True, filename = I_filename)
ell = gm.create_ellipsoid_data()
ellipsoid = {1: ell}  

RESTR_TYPES = {}; #
RESTR_NAME = [{'NAME': 'DOMAIN_EXCL'}]  #
for i in range(len(RESTR_NAME)): RESTR_TYPES[i] = RESTR_NAME[i] #
Restr = RESTR_TYPES #
keyName = 'DOMAIN_EXCL' 
cutoff =  opt["RESTRAINTS"]["DOMAIN_EXCL"]["STRENGTH"]     
rs = AFM.EllipsoidExcludedVolumeRestraint(ell,cutoff=cutoff)
Restr[0]['RESTR'] = rs #2 dicts: dictionary of a dictionary
Restr[0]['EVAL'] = True
Restr[0]['INPUTKEY'] = "ellipsoid" #
Restr[0]['SCORE'] = None #

R_filename = sys.argv[3]
# R = "End_814.gmm"
gm2 = AFM.GMData(fromFile=True, filename = R_filename)
gm2.init_ctypes_of_gmm()
mat2 = gm2.call_get_abs_struct_factor_all_mesh_on_Ewald_sphere(200,det_width,det_dist,wave_length)
referenceMat= np.log10(mat2)

dw = det_width  # change k: 2019.4.1
ds = det_dist
wave = wave_length
cd = float(opt["circle_width"])
lowest_k = float(opt["lowest_k"])
k_inc = 1.0/wave
k= (k_inc*dw/np.sqrt(ds*ds + 2*dw*dw))
ntheta=360
rlist = np.arange(lowest_k,k+cd,cd)
dtheta=2*np.pi/ntheta
theta_list=np.arange(0.0,2*np.pi,dtheta)
coslist,sinlist=np.cos(theta_list),np.sin(theta_list)
f2list=np.zeros((len(rlist),len(theta_list)))
for i_r in range(len(rlist)):
  for each in range(len(theta_list)):
    kx,ky,kz=AFM.TOOLS.get_coord_on_ewald_sphere_from_zero_curvature_coord(rlist[i_r]*coslist[each],rlist[i_r]*sinlist[each],k_inc)
    f2list[i_r][each]=gm2.call_get_abs_struct_factor(kx,ky,kz)
f2list=np.log10(f2list)*2.0

movers_opt = opt["OPTIMIZE"]["MOVERS"] # debug
mover_level = 0
ndom = len(ell)
if movers_opt['RESIDUE'] == 'None':
  movers_opt['RESIDUE'] = []
if movers_opt['DOMAIN'] == 'All':
  movers_opt['DOMAIN'] = range(1,ndom+1)
  mover_level += 2
elif movers_opt['DOMAIN'] == 'None':
  movers_opt['DOMAIN'] = []  
else:
  mlist = []
  movers_opt['DOMAIN'] = str(movers_opt['DOMAIN'])
  for m in movers_opt['DOMAIN'].split(","):
    try:
      im = int(m)
    except ValueError:
      raise Exception("Mover List Domain includes non-integer")
    if im > ndom or im < 1:
      raise Exception("Mover List Domain improper")
    mlist.append(im)
  movers_opt['DOMAIN'] = mlist
  mover_level += 2
if movers_opt['NPROPOSE'] == 'Auto':
  if mover_level == 1: 
    movers_opt['NPROPOSE'] = len(movers_opt['RESIDUE'])
  elif mover_level == 2:
    movers_opt['NPROPOSE'] = len(movers_opt['DOMAIN'])
  else:
    l1 = len(movers_opt['RESIDUE'])
    l2 = len(movers_opt['DOMAIN'])
    movers_opt['NPROPOSE'] = int((l1 + l2)/2)
if movers_opt['NPROPOSE'] == 'All':
  if mover_level == 1: 
    movers_opt['NPROPOSE'] = len(movers_opt['RESIDUE'])
  elif mover_level == 2:
    movers_opt['NPROPOSE'] = len(movers_opt['DOMAIN'])
  else:
    l1 = len(movers_opt['RESIDUE'])
    l2 = len(movers_opt['DOMAIN'])
    movers_opt['NPROPOSE'] = l1 + l2
if movers_opt['NPROPOSE'] == 'None':
  raise Exception("Nothing to move, modify NPROPOSE")
mlist = list(movers_opt["DOMAIN"]);  del movers_opt["DOMAIN"];  movers_opt[1] = mlist;
mlist = list(movers_opt["RESIDUE"]); del movers_opt["RESIDUE"]; movers_opt[0] = mlist;
opt["OPTIMIZE"]["MOVERS"] = movers_opt

gopt = AFM.OptionAFM()
nerr = gopt.reset_options(opt) # RESET options by the input "opt"

Data =  {"ellipsoid": ellipsoid,"refernMat": referenceMat, "restr": Restr,"f2list": f2list}
mcs = AFM.MCSampler(gopt.opt,Data,timing=False)
mcs.start()

