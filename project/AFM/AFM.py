from __future__ import print_function
from AFM_Header import *

import ctypes

import OtherTools as TOOLS # 12/4/2019

clib_dir="/home/asi/Desktop/agk/AFM"


def joke():
    return markdown(u'Nothing new, If you are not part of the solution,'
                    u'you are part of the precipitate.')

# General public functions


def printf(str, *args):
    """
      A C-like printf
    """
    print(str % args, end='', flush = True)

def invert_decomposition(rot,radii):
    """
      Recalculate undecomposed matrix from radii and rotation matrix
      rot: Rotation matrix (3x3)
      radii: Radii corresponding to matrix
      returns: Matrix that results input rotation matrix and radii
    """
    val = radii * radii
    cv = (rot.dot(np.diag(val))).dot(rot.T)
    return cv

def checkRotation(self,rot):
    """
      Check Rotation matrix (proper rotation and fix chirality)
      rot: Input rotation matrix (3x3)
      returns: Output refined rotation matrix
    """
    if np.linalg.det(rot) < 0.0: rot = -rot
    if sum(np.cross(rot[:,0],rot[:,1]) * rot[:,2]) < 0.0: rot[:,2] = -rot[:,2]
    return rot

def sort_decomposition(self,rot,radii,decreasing=False):
    """
      Sort decomposed rotation matrix and radii, according to
      increasing/decreasing order of radii
      rot: Input rotation matrix (3x3)
      radii: Input radii (vector of length 3)
      decreasing: Order of radii.
      returns: tuple of rotation matrix and radii
    """
    if decreasing: 
      windex = np.argsort(-radii)
    else:
      windex = np.argsort(radii)
    radii = radii[windex]
    rot = rot[:,windex]
    return rot,radii
     
class PDB:

  """
    This class is for reading pdb file and transforming its
    coordinates if required
  """

  def __init__(self,filename,asis=True,transf=None):
    """
      PDB class constructor

      filename: Name of the PDB file
      asis: Whether PDB coordinates are transformed or not
      transf: The transformation information

      The asis and transf arguments dictates how PDB coordinates
      are transformed. 
      If asis = True, no transformation is applied to the PDB 
      coordinates. If asis = False, then a transformation is
      applied, the transformation can be specified by the argument
      transf (see transform_coord). When asis = False, and transf 
      = None, then an approx. transformation is applied. See below.

      Normally when a PDB coordinate file is given which will be used
      later to compare with a AFM-like image, we don't know how 3D
      coordinates should be transformed to match to the AFM-image. E.g.
      consider a case when a AFM-like image has a shape and it is in the
      middle of the image with some finite border. So, in the 2D projected
      plane the imaged object has nonzero x,y-values for the topleft and 
      bottomright pixels. Those xy-values for the topleft pixel define a
      translation along xy-direction that should be applied to the PDB
      coordinates to get a initial approximate 3D coordinates. A reasonable
      guess for z-translation is at least minimum z-value should be zero, 
      then the 3D object has a positive z(max), and max(z) - min(z) can be 
      easily mapped to a color. If asis = False, transf = None, this 
      z-translation is applied (see reset_coord). 

    """
    
    if not isfile(filename):
      raise Exception("PDBFILE %s not found" % filename)
    self.filename = filename
    self.object = self.read_pdb(filename)

    coord = self.get_coordinates(self.object)    
    self.objMin,self.objMax,self.objLen = self.get_size(coord)

    if not asis:

      if transf is None:
        self.object,coord = self.reset_coord(self.object,coord)

      else:

        print("transforming\n")
        rot = transf.rot
        trans = transf.trans-self.objMin
        coord = self.transform_coord(coord,trans=trans,rot=rot)
        self.object,coord = self.set_coordinates(self.object,coord)

    self.objMin,self.objMax,self.objLen = self.get_size(coord)
    
  def read_pdb(self,filename):
    """
      This function reads the pdbfile using Bio.PDB.PDBParser
      filename: Name of the PDB file
      return: a structure (3D atomic model)
    """
    return Bio.PDB.PDBParser().get_structure(file=filename,id='')    
  
  def write_pdb(self,structure,filename):
    """
      This function writes the pdbfile using Bio.PDB.PDBIO()
      structure: The structure to be written
      filename: Name of the output file
    """
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(filename)    

  def reset_coord(self,structure,coord):
    """
      This function resets the coordinate such that it's min(x), min(y),
      and min(z) are at the origin (0,0,0).
      structure: The structure as parsed by Bio.PDB
      coord: Coordinate matrix
             (x,y,z values, shape N x 3, N = number of atoms)
      returns: reset coordinates
    """
    print("resetting\n")
    cmin, cmax, clen = self.get_size(coord)    
    coord = self.transform_coord(coord,trans=-cmin,rot=np.diag([1.0,1.0,1.0]))    
    return self.set_coordinates(structure,coord)

  def transform_coord(self,vertices,trans,rot):
    """
      This transforms the coordinate matrix given translation and rotation
      vertices: Set of 3D points
                (x,y,z values, shape N x 3, N = number of atoms)
      trans: Translation vector (3D)
      rot: Rotation matrix (3x3)
      returns: Transformed coordinates
    """
    vertices_rotated = self.rotate_coord(vertices,rot)
    return self.translate_coord(vertices_rotated,trans)

  def rotate_coord(self,vertices,rot):
    """
      This rotates the coordinate matrix given rotation matrix
      vertices: Set of 3D points
                (x,y,z values, shape N x 3, N = number of atoms)
      rot: Rotation matrix (3x3)
      returns: Rotated coordinates
    """
    # In order to rotate, first translate the coordinates to the origin,
    # then rotate, and then back-translate the coordinates
    center = np.mean(vertices,axis=0)
    verticesOrigin = self.translate_coord(vertices,-center)
    verticesRotatedOrigin = (rot.dot(verticesOrigin.T)).T
    return self.translate_coord(verticesRotatedOrigin,center)

  def translate_coord(self,vertices,trans):
    """
      This translates the coordinate matrix given translation vector
      vertices: Set of 3D points
                (x,y,z values, shape N x 3, N = number of atoms)
      trans: Translation vector (3D)
      returns: Translated coordinates
    """
    return vertices + trans

  def get_coordinates(self,structure):
    """
      This creates the coordinates matrix
      structure: The structure as parsed by Bio.PDB
      returns: The coordinate
               (x,y,z values, shape N x 3, N = number of atoms)
    """
    return np.array([a.get_coord() for a in structure.get_atoms()])

  def set_coordinates(self,structure,coord):
    """
      This updates coordinates in the input structure using input
      coordinates.
      structure: The structure as parsed by Bio.PDB
      coord: The coordinate
             (x,y,z values, shape N x 3, N = number of atoms)
      returns: The structure with new coordinate and the coordinate itself
    """
    k = 0
    for a in structure.get_atoms():
      a.set_coord(coord[k,:])
      k += 1    
    return (structure,coord)

  def get_size(self,coord):    
    """
      This calculates size of the structure (3D atomic model) from range of
      the coordinates.
      coord: The coordinate
             (x,y,z values, shape N x 3, N = number of atoms)
      returns: tuple of min(xyz), max(xyz) and delta(xyz)
    """
    coordMin = np.min(coord,axis=0)
    coordMax = np.max(coord,axis=0)
    coordDel = coordMax - coordMin
    return (coordMin,coordMax,coordDel)

class Representation:
  """

    This class creates ellipsoid representation of the 3D atomic model. 
    It is assumed that the 3D atomic models is composed of N number of ellipsoids,
    where each ellipsoid is associated to a domain (i.e. the 3D atomic model is
    composed of N domains). This set of ellipsoids is similar to a mixture of 
    Gaussian probability distribution (but not identical, because a Gaussian
    Mixture Model is an an optimized set of Gaussians done by EM-algorithm. 
    Thus the set of ellipsoids here is unoptimized set of Gaussians.
    To make this set of ellipsoids really like a mixture of Gaussians, at least 
    we have to define a mixture weight of each Gaussian (or ellipsoid). Without
    any optimization such mixture weight is defined here uniformly.
    It should be noted that here we assume that we know domain-decomposition of the
    given PDB structure and we affix an ellipsoid to each domain, and avoid any
    optimization.
    The domain can be defined as a sequence of residues (with more than one residues).
    It is also possible to define one ellipsoid per residue (i.e. very detailed 
    volumetric representation). Such settings can be defined in input yaml file.

  """

  def __init__(self,pdb,opt):
    """
      This is the constructor of the class Representation.
      pdb: The input PDB object of class PDB or None.
      opt: A dictionary of options.
      Particularly, here 'ELLIPSOID' key in opt dictionary is useful.
      ELLIPSOID: NLEVEL: can be either 1 or 2. 
                 If NLEVEL is 1, ELLIPSOID->LEVEL defines the granularity of the
                 ellipsoids. If ELLIPSOID->LEVEL = 'residue', one residue is 
                 mapped to one ellipsoid. If ELLIPSOID->LEVEL = 'domain', then
                 one domain is mapped to one ellipsoid.
                 If NLEVEL is 2, it is a multiscale model. This means in lowest
                 level each residue is mapped to a ellipsoid, and in the next upper
                 level a set of residue-level ellipsoids  (or domain) is mapped to a
                 domain-level ellipsoids. Here domain-level ellipsoids are built
                 from residue-level ellipsoids. [see nested_domain_ellipsoids]
                 [This multiscale model representation currently is not used in
                  sampling]
      ELLIPSOID: SCALESD: a scalling(s) which is applied to the standard-deviations of
                 the ellipsoids.
      ELLIPSOID: WEIGHT: Whether each ellipsoid members are weighed by mass, or by volume.
                 The volume-based weighing is not implemented for residue-level ellipsoids.
      Also, the key 'DOMAIN' in opt dictionary is used to define DOMAINS.
      DOMAIN->SEGMENT is a list of residue-sets.
    """

    self.opt = opt
    level_dict = {'residue': 0,'domain': 1}

    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    scalesd = self.opt["ELLIPSOID"]["SCALESD"]
    weight = self.opt["ELLIPSOID"]["WEIGHT"]
    level = self.opt["ELLIPSOID"]["LEVEL"]

    if pdb == None:
      input_fmt = self.opt['INPUTFILE']['FORMAT']
      input_objfile = self.opt['INPUTFILE']['NAME']
      if input_fmt == 'gmm':
        gmm = GMData(fromFile=True,filename=input_objfile)
        npt = self.opt['ELLIPSOID']['NUM_DATA_POINT']
        self.DomEllipsoids = gmm.create_ellipsoid_data(scale=scalesd[0],npoint=npt)
      else:
        raise Exception("Error: unknown inputfile format")
      return

    level_int = [None] * nlevel
    for l in range(nlevel):
      level_value = level[l]
      if level_dict.has_key(level_value):
        level_int[l] = level_dict[level_value]
      else:
        raise Exception("Option ELLIPSOID/LEVEL includes invalid name %s" % level_value)

    if nlevel > 1:
      # Multiscale
      for l in range(1,nlevel):
        lthis = level_int[l]; lprev = level_int[l-1];
        if lthis <= lprev:
          raise Exception("Option ELLIPSOID/LEVEL includes invalid hierarchy %s <= %s"\
                           % level_value[l],level_value[l-1])
      if nlevel == 2:
        # First set residue-level ellipsoids
        self.residue_ellipsoids(pdb,scalesd[0],weight[0])
        # Then set domain-level ellipsoids and nested set of residue-level ellipsoids
        self.nested_domain_ellipsoids(pdb,scalesd[1],weight[1])
    elif nlevel == 1:
      # Singlescale
      if level[0] == 'residue':
        # Residue-level ellipsoids
        self.residue_ellipsoids(pdb,scalesd[0],weight[0])
      if level[0] == 'domain':
        # Domain-level ellipsoids
        self.domain_ellipsoids(pdb,scalesd[0],weight[0])
  
  def set(self,ellipsoids):
    """
      OBSOLETE
      This function copies input ellipsoids to the class variables.
      ellipsoids: The input ellipsoids that will be copied.      
    """

    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]

    if nlevel == 2:

      if level[0] == 'residue' and level[1] == 'domain':
        self.ResEllipsoids = copy.copy(ellipsoids[0])
        self.DomEllipsoids = copy.copy(ellipsoids[1])
      if level[0] == 'domain' and level[1] == 'residue':
        self.DomEllipsoids = copy.copy(ellipsoids[0])
        self.ResEllipsoids = copy.copy(ellipsoids[1])

    elif nlevel == 1:

      if level[0] == 'residue':
        self.ResEllipsoids = copy.copy(ellipsoids)
      if level[0] == 'domain':
        self.DomEllipsoids = copy.copy(ellipsoids)

  def residue_ellipsoids(self,pdb,sd,w):
    """
      This creates residue-level ellipsoids
      pdb: The input PDB object of type class PDB
      sd: The factor of standard-deviation for ellipsoids
      w: 'mass' or 'volume' for mass or volume-based weighing of constituent
          atoms. 'volume'-based weighing is not implemented in this case.
    """

    opt = self.opt
    # Following used Bio.PDB get_residues, pdb.object is Bio.PDB object
    self.residues = [r for r in pdb.object.get_residues()]

    # ResEllipsoids = Array of Ellipsoids for each residue
    ResEllipsoids = [None] * len(self.residues) 
    if w == 'volume': 
      print("Radius of atoms not known, falling back to mass-weighting")

    for ir,r in enumerate(self.residues):

      # This uses Ellipsoid-class
      ResEllipsoids[ir] = Ellipsoid()

      # Setup coordinate-matrix from constituent atoms
      coord = []; weight = []      
      for a in r.get_atoms():
        coord.append(a.get_coord())        
        weight.append(a.mass)
      coord = np.array(coord)

      # Now set ellipsoid data
      ResEllipsoids[ir].set_data(data = coord,mass = weight,scale = sd)
      ResEllipsoids[ir].calculate_aux()

    nres = len(ResEllipsoids)
    print("Number of residue ellipsoids: %d" % nres)    

    # Copy to class-variable
    self.ResEllipsoids = copy.copy(ResEllipsoids)

  def create_domainDict(self,segments):
    """
      This is a helper function to create domain-level ellipsoids
      segments: The segments as defined in DOMAIN->SEGMENT
      returns: A dictionary with domain-level information
    """

    # Initialization    
    ndom = len(segments)

    # domainDic is a hash-of-hash-of-hash. Topmost key is integer index 
    # of domain. Middle key is index of residue within the segment.
    # Lowest key contains residue-name, status, and residue-index.
    # status = True, if corresponding residue exists in PDB, else, False.    
    domainDic = {}    

    for k in range(ndom):      
      domainDic[k] = [None] * len(segments[k])
      for ir,rname in enumerate(segments[k]):
        status = False; index = None
        for ires,res in enumerate(self.residues):
          if res.get_id()[1] == rname:
            status = True; index = ires
        domainDic[k][ir] = {'name': rname, 'status': status, 'index': index}

    # Check for conflicts
    for i in range(ndom):
      di = domainDic[i]
      for j in range(i+1,ndom):
         dj = domainDic[j]
         for dii in di:
           for djj in dj:
             if dii['name'] == djj['name']:
               print("Warning: Conflict in segment definition for %d,%d : %d %d"\
                      % (i,j,dii['name'],djj['name']))
               djj['status'] = False

    return domainDic

  def nested_domain_ellipsoids(self,pdb,sd,w):
    """
      This is domain-ellipsoids builder from residue-ellipsoids
      pdb: The input PDB object of type class PDB
      sd: The factor of standard-deviation for ellipsoids
      w: 'mass' or 'volume' for mass or volume-based weighing of constituent
          atoms.
    """

    opt = self.opt
    segments = opt["DOMAIN"]["SEGMENT"]
    ndom = len(segments)

    # Create domainDic information
    domainDic = self.create_domainDict(segments)        

    self.DomEllipsoids = [None] * ndom
    # Residue-level scalesd
    rscale = opt["ELLIPSOID"]["SCALESD"][0]

    for k in range(ndom):  

      d = domainDic[k]
      coord = []; weight = []; dmass = [];

      # This uses Ellipsoid-class
      self.DomEllipsoids[k] = Ellipsoid() 
    
      # Set-up data for ellipsoid definition
      ResMeshMin = []; ResMeshMax = [];
      for delm in d:
        if delm['status']:
          i = delm['index']
          coord.append(self.ResEllipsoids[i].mean)
          rmin = self.ResEllipsoids[i].dataMin     
          rmax = self.ResEllipsoids[i].dataMax
          ResMeshMin.append(rmin); ResMeshMax.append(rmax);
          if w == 'volume':
            rad = (rmax - rmin) * 0.5 * rscale
            weight.append(np.prod(rad))
          elif w == 'mass':
            mw = self.ResEllipsoids[i].mass
            weight.append(mw)

      # Define ellipsoid based on data
      self.DomEllipsoids[k].set_data(data = coord,mass = weight,scale = sd)
      self.DomEllipsoids[k].calculate_aux()

      # Set-up range of data for each ellipsoid
      ResMeshMin = np.array(ResMeshMin); ResMeshMax = np.array(ResMeshMax);
      dmin = np.min(ResMeshMin,axis=0)
      dmax = np.max(ResMeshMax,axis=0)
      self.DomEllipsoids[k].dataMin = dmin
      self.DomEllipsoids[k].dataMax = dmax     
     
    print("Number of domain ellipsoids: %d" % ndom)

  def domain_ellipsoids(self,pdb,sd,w):
    """
      This is domain-ellipsoids builder from PDB
      pdb: The input PDB object of type class PDB
      sd: The factor of standard-deviation for ellipsoids
      w: 'mass' or 'volume' for mass or volume-based weighing of constituent
          atoms.
    """

    opt = self.opt
    segments = opt["DOMAIN"]["SEGMENT"]
    ndom = len(segments)

    # Following used Bio.PDB get_residues, pdb.object is Bio.PDB object
    self.residues = [r for r in pdb.object.get_residues()]

    # Create domainDic information
    domainDic = self.create_domainDict(segments)  
      
    self.DomEllipsoids = [None] * ndom
    for k in range(ndom):

      d = domainDic[k]
      coord = []; weight = [];

      # This uses Ellipsoid-class    
      self.DomEllipsoids[k] = Ellipsoid()
      
      # Set-up data for ellipsoid definition
      for delm in d:
        if delm['status']:
          i = delm['index']                
          for a in self.residues[i].get_atoms():
            coord.append(a.get_coord())   
            weight.append(a.mass)
      coord = np.array(coord)

      # Define ellipsoid based on data
      self.DomEllipsoids[k].set_data(data = coord,mass = weight,scale = sd)
      self.DomEllipsoids[k].calculate_aux()

      # Set-up range of data for each ellipsoid
      dmin = np.min(coord,axis=0); dmax = np.max(coord,axis=0);
      self.DomEllipsoids[k].dataMin = dmin
      self.DomEllipsoids[k].dataMax = dmax

    print("Number of domain ellipsoids: %d" % ndom)

  def _get_size(self,mp):
    """
      Get size of the mesh-representation of ellipsoids
      mp: mesh points
    """        
    objMin = np.min(mp,axis=0)
    objMax = np.max(mp,axis=0)
    objLen = objMax - objMin
    return (objMin,objMax,objLen)

  def _get_vert_elllipsoid(self):
    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]       
    if nlevel == 2:
      nres = len(self.ResEllipsoids)                  
      vert = []
      for i in range(nres):
        vert.append(self.ResEllipsoids[i].data[:,:])
      vert = np.concatenate(vert,axis=0)  
    if nlevel == 1:
      if level[0] == 'domain':        
        ndom = len(self.DomEllipsoids)
        vert = []
        for i in range(ndom):          
          vert.append(self.DomEllipsoids[i].data[:,:])          
        vert = np.concatenate(vert,axis=0)        
      elif level[0] == 'residue':
        nres = len(self.ResEllipsoids)        
        vert = []
        for i in range(nres):
          vert.append(self.ResEllipsoids[i].data[:,:])
        vert = np.concatenate(vert,axis=0)
    return vert

  def _get_center_ellipsoid(self):
    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]
    if nlevel == 2:
      nres = len(self.ResEllipsoids)
      vert = np.empty([nres,3])
      for i in range(nres):
        vert[i,:] = self.ResEllipsoids[i].mean[:]
    if nlevel == 1:
      if level[0] == 'domain':
        ndom = len(self.DomEllipsoids)
        vert = np.empty([ndom,3])
        for i in range(ndom):
          vert[i,:] = self.DomEllipsoids[i].mean[:]
      if level[0] == 'residue':
        nres = len(self.ResEllipsoids)
        vert = np.empty([nres,3])
        for i in range(nres):
          vert[i,:] = self.ResEllipsoids[i].mean[:]
    return vert

  def _transform_data_representation(self,trans,rot):
    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]
    if nlevel == 2:
      nres = len(self.ResEllipsoids)      
      for i in range(nres):
        vert = self.ResEllipsoids[i].data[:,:]
        cen = np.mean(vert,axis=0)
        vert = (rot.dot((vert - cen).T)).T + cen + trans
        self.ResEllipsoids[i].data =  vert[:,:]
        self.ResEllipsoids[i].dataMin = dmin
        self.ResEllipsoids[i].dataMax = dmax
      ndom = len(self.DomEllipsoids)
      for i in range(ndom):
        vert = self.DomEllipsoids[i].data[:,:]
        cen = np.mean(vert,axis=0)
        vert = (rot.dot((vert - cen).T)).T + cen + trans
        self.DomEllipsoids[i].data = vert[:,:]
        dmin = np.min(vert,axis=0); dmax = np.max(vert,axis=0);
        self.DomEllipsoids[i].dataMin = dmin
        self.DomEllipsoids[i].dataMax = dmax
    elif nlevel == 1:
      if level[0] == 'residue':
        nres = len(self.ResEllipsoids)
        for i in range(nres):
          vert = self.ResEllipsoids[i].data[:,:]
          cen = np.mean(vert,axis=0)
          vert = (rot.dot((vert - cen).T)).T + cen + trans
          self.ResEllipsoids[i].data =  vert[:,:]
          dmin = np.min(vert,axis=0); dmax = np.max(vert,axis=0);
          self.ResEllipsoids[i].dataMin = dmin
          self.ResEllipsoids[i].dataMax = dmax
      elif level[0] == 'domain':
        ndom = len(self.DomEllipsoids)
        for i in range(ndom):
          vert = self.DomEllipsoids[i].data[:,:]
          cen = np.mean(vert,axis=0)
          vert = (rot.dot((vert - cen).T)).T + cen + trans
          self.DomEllipsoids[i].data = vert[:,:]
          dmin = np.min(vert,axis=0); dmax = np.max(vert,axis=0);
          self.DomEllipsoids[i].dataMin = dmin
          self.DomEllipsoids[i].dataMax = dmax

  def _transform_mean_representation(self,trans,rot):
    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]
    print(trans,rot)
    if nlevel == 2:
      nres = len(self.ResEllipsoids)
      vert_res = np.zeros((nres,3))
      for i in range(nres):
        vert_res[i,:] = self.ResEllipsoids[i].mean[:]
      cen = np.mean(vert_res,axis=0)
      vert_res = (rot.dot((vert_res - cen).T)).T + cen + trans
      for i in range(nres):
        self.ResEllipsoids[i].mean[:] = vert_res[i,:]

      ndom = len(self.DomEllipsoids)
      vert_dom = np.zeros((ndom,3))
      for i in range(ndom):
        vert_dom[i,:] = self.DomEllipsoids[i].mean[:]
      cen = np.mean(vert_dom,axis=0)
      vert_dom = (rot.dot((vert_dom - cen).T)).T + cen + trans      
      for i in range(ndom):
        self.DomEllipsoids[i].mean[:] = vert_dom[i,:]

    elif nlevel == 1:
      if level[0] == 'residue':
        nres = len(self.ResEllipsoids)
        vert_res = np.zeros((nres,3))
        for i in range(nres):
          vert_res[i,:] = self.ResEllipsoids[i].mean[:]
        cen = np.mean(vert_res,axis=0)
        vert_res = (rot.dot((vert_res - cen).T)).T + cen + trans
        for i in range(nres):
          self.ResEllipsoids[i].mean[:] = vert_res[i,:]

      elif level[0] == 'domain':
        ndom = len(self.DomEllipsoids)
        vert_dom = np.zeros((ndom,3))
        for i in range(ndom):
          vert_dom[i,:] = self.DomEllipsoids[i].mean[:]        
        cen = np.mean(vert_dom,axis=0)
        vert_dom = (rot.dot((vert_dom - cen).T)).T + cen + trans      
        for i in range(ndom):
          self.DomEllipsoids[i].mean[:] = vert_dom[i,:]

  def apply_transformation(self,asis,transf):
    """
      Transform ellipsoids centers and data by the argument 'transf'
      transf: Transformation (shift + rotation)      
    """
    nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
    level = self.opt["ELLIPSOID"]["LEVEL"]

    self.SetupMesh(filename=None,do_projection=False)
    mp = self.mesh.meshPoints
    vert  = self._get_vert_elllipsoid()
    cvert = self._get_center_ellipsoid()

    self.repMin,self.repMax,self.repLen = self._get_size(mp)
    self.datMin,self.datMax,self.datLen = self._get_size(vert)
    self.cenMin,self.cenMax,self.cenLen = self._get_size(cvert)

    print("Mesh points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.repMin[0],self.repMin[1],self.repMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.repMax[0],self.repMax[1],self.repMax[2]))
    print("Data points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.datMin[0],self.datMin[1],self.datMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.datMax[0],self.datMax[1],self.datMax[2]))
    print("Center points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.cenMin[0],self.cenMin[1],self.cenMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.cenMax[0],self.cenMax[1],self.cenMax[2]))

    if not asis:
      if transf is None:
        rot = np.diag([1.0,1.0,1.0])
        trans1 = -self.repMin
        trans2 = -self.repMin
      else:
        rot = transf.rot        
        trans1 = transf.trans - self.repMin
        trans2 = transf.trans - self.repMin
      self._transform_data_representation(trans1,rot)
      self._transform_mean_representation(trans2,rot)

    self.SetupMesh(filename=None,do_projection=False)
    mp = self.mesh.meshPoints
    vert  = self._get_vert_elllipsoid()
    cvert = self._get_center_ellipsoid()
    self.repMin,self.repMax,self.repLen = self._get_size(mp)
    self.datMin,self.datMax,self.datLen = self._get_size(vert)
    self.cenMin,self.cenMax,self.cenLen = self._get_size(cvert)
    print("Mesh points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.repMin[0],self.repMin[1],self.repMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.repMax[0],self.repMax[1],self.repMax[2]))
    print("Data points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.datMin[0],self.datMin[1],self.datMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.datMax[0],self.datMax[1],self.datMax[2]))
    print("Center points bounding box detail:")
    print("  Min: %.3f %.3f %.3f" % (self.cenMin[0],self.cenMin[1],self.cenMin[2]))
    print("  Max: %.3f %.3f %.3f" % (self.cenMax[0],self.cenMax[1],self.cenMax[2]))
        
  def SetupMesh(self,filename=None,do_projection=True,dpi=None):
    """
      This sets up Mesh representation from ellipsoid data. This is a
      wrapper method to use Mesh class instance to build mesh, then
      AFM-like projection, save the projection, and return projected
      matrix.
      filename: The output file to save the AFM-like projection
      A mesh-representation is built on ellipsoid data, where each
      ellipsoid is wrapped by a surface cut at a specific distance.
      The meshpoints are those surface points.
      Here, MESH key from input option settings are used.
      MESH->LEVEL defines which level of representation to be used to
      build mesh. This should match (include) ELLIPSOID->LEVEL.
    """

    opt = self.opt
    meshLevel = opt["MESH"]["LEVEL"]
    self.mesh = Mesh(opt)

    # Build the mesh from domain-level or residue-level ellipsoids.
    if meshLevel == 'domain':
      self.mesh.buildMesh(self.DomEllipsoids)      
    elif meshLevel == 'residue':
      self.mesh.buildMesh(self.ResEllipsoids)
    else:
      raise Exception("No appropriate mesh-level found")

    if not do_projection:
      return

    # Build the AFM-like projection
    if opt["MODE"] == "Debug": 
      self.mesh.buildProjection(num_thread=opt["NTHR"],timing=True)
    else: 
      self.mesh.buildProjection(num_thread=opt["NTHR"],timing=False)

    # Save the projection
    self.mesh.saveProjection(filename=filename,dpi=dpi)    

    # Get projection and digitize the projected matrix depending
    # on IMGCOMPARE->MODE.
    # Also here, zbins defined in mesh are used, which are already set
    # in opt. That means, before calling this function make sure that
    # zbins are defined properly.    
    self.colorMat = self.mesh.get_projected_matrix(invert=False)
    if opt["IMGCOMPARE"]["MODE"] == "up_is_dark":
      self.projMat = self.mesh.get_projected_matrix(invert=True)      
      self.colorMat = 255 - (np.digitize(self.colorMat,self.mesh.zbins) - 1)
      print("made with up is dark")
    elif opt["IMGCOMPARE"]["MODE"] == "up_is_light":
      self.projMat = self.mesh.get_projected_matrix(invert=False)
      self.colorMat = np.digitize(self.colorMat,self.mesh.zbins) - 1
      print("made with up is light")
    else:
      raise Exception("No appropriate IMGCOMPARE/MODE found")

    # Set dpi parameter
    self.opt['dpi'] = self.mesh.imageParam['dpi']

  def dumpMeshVRML(self,filename):
    """
      This function write mesh information/ellipsoid information as VRML
      filename: Output VRML file 
    """

    self.mesh.buildVRML()
    self.mesh.writeVRML(filename)

  def dumpMeshGMM(self,filename):
    """
      This function write mesh information/ellipsoid information as gmm-file
      GMM-file is a modified flat PDB-like format.
      filename: Output GMM file 
    """

    self.mesh.buildGMData()
    self.mesh.writeGMData(filename)    

class Transformation:
  """
    This class is a container for rotation or translation
  """
  def __init__(self,rot,trans,rotationType='q'):
    """
      The constructor that takes rotation and translation
      information and optional type of rotation. 'q' indicates
      quaternion.
      rot: rotation matrix (rotationType='m') 
           or vector (rotationType='q') 
      trans: translation vector
      rotationType: 'q' or 'm'. 'q' means quaternion (4D vector),
                     'm' means 3x3 rotation matrix
    """
    self.rot = rot
    self.trans = trans
    self.rotationType = rotationType

  def get_rotation(self):
    """
      Return saved rotation
      return: rotation
    """
    return self.rot

  def get_translation(self):
    """
      Return saved translation
      return: translation
    """
    return self.trans

  def set_rotation(self,rot):
    """
      Sets input rotation
      rot: Input rotation
    """
    self.rot = rot

  def set_translation(self,trans):
    """
      Sets input translation
      trans: Input transtion
    """
    self.trans = trans

  def change_rotation(self,fromtype='q',totype='m'):
    """
      This converts between two types of representation of
      rotation. 'q' to 'm' or 'm' to 'q'.
      fromtype: 'q' or 'm' Saved type of rotation
      totype: 'q' or 'm' converted type of rotation
      return: converted rotation
    """

    if fromtype == self.rotationType:
      if totype == fromtype:
        self.get_rotation()
      elif totype == 'm':
        return self.QuaternionToMatrix(self.rot)
      elif totype == 'q':
        return self.MatrixToQuaternion(self.rot)
      else:
        return None
    else:
      return None

  def QuaternionToMatrix(self,rot):
    """
      Conversion from quaternion representation of rotation to matrix
      rot: Input quaternion
      return: Output rotation matrix
    """
    R = np.zeros((3,3))
    a = rot[0]; b = rot[1]; c = rot[2]; d = rot[3];
    R[0,0] = 1.0 - 2.0*c*c - 2.0*d*d; R[0,1] = 2.0*b*c - 2.0*a*d;       R[0,2] = 2.0*b*d + 2.0*a*c;
    R[1,0] = 2.0*b*c + 2.0*a*d;       R[1,1] = 1.0 - 2*b*b - 2*d*d;     R[1,2] = 2.0*c*d - 2.0*a*b;
    R[2,0] = 2.0*b*d - 2.0*a*c;       R[2,1] = 2.0*c*d + 2.0*a*b;       R[2,2] = 1.0 - 2.0*b*b - 2.0*c*c;
    return R

  def MatrixToQuaternion(self,rot):
    """
      Conversion from matrix representation of rotation to quaternion
      rot: Input rotation matrix
      returns: Output quaternion
    """
    # This follows IMP-2.8's src, modules/algebra/src/Rotation3D.cpp, line:73-89
    m11 = rot[0,0]; m12 = rot[0,1]; m13 = rot[0,2];
    m21 = rot[1,0]; m22 = rot[1,1]; m23 = rot[1,2];
    m31 = rot[2,0]; m32 = rot[2,1]; m33 = rot[2,2];
    a = np.abs(1.0 + m11 + m22 + m33) / 4.0
    b = np.abs(1.0 + m11 - m22 - m33) / 4.0
    c = np.abs(1.0 - m11 + m22 - m33) / 4.0
    d = np.abs(1.0 - m11 - m22 + m33) / 4.0
    sumq = a + b + c + d
    a = np.sqrt(a / sumq);
    b = np.sqrt(b / sumq);
    c = np.sqrt(c / sumq);
    d = np.sqrt(d / sumq);
    if (m32 - m23 < 0.0): b = -b;
    if (m13 - m31 < 0.0): c = -c;
    if (m21 - m12 < 0.0): d = -d;
    return (a,b,c,d) 

  def QuaternionToAxisAngle(self,rot):
    """
      This is an useful converter function to quaternion to axis-angle
      representation.
      rot: Input quaternion
      returns: Output axis-angle representation of rotation (also 4D vector)
    """
    a = rot[0]; b = rot[1]; c = rot[2]; d = rot[3];
    if np.abs(a) > 0.9999: 
      return (1.0,0.0,0.0,0.0)
    else:
      angle = np.arccos(a) * 2.0
      sineValue = np.sin(angle / 2.0)
      x = b / sineValue; y = c / sineValue; z = d / sineValue;
      vec = np.array([x,y,z])
      normvec = np.linalg.norm(vec)
      if normvec > 1e-4:
        vec = vec / normvec
      return (vec[0],vec[1],vec[2],angle)     
      
class Ellipsoid:
  """
    This class defines each ellipsoid. An ellipsoid (3D) is a 3D Gaussian
    probability distribution. An ellipsoid is defined by a mean (the center
    of Gaussian), a covariance matrix (3x3). However, to avoid repeatative
    computation also a number of other auxilliary parameters are defined, 
    e.g. rotation matrix, radii, scaled radii, quaternion representation of
    rotation matrix, a multivariate probability distribution calculator 
    class.
    The class parameter variable chi2CoverRatio is set to 0.828203. This is
    used to set a cutoff to the probability distribution around the center
    according to a inverse chi-square distribution.
  """
  tol = 1e-8
  chi2CoverRatio = 0.828203
  def __init__(self):
    """
      The class constructor. This just defines a few internal variable to
      be used in later calls.
    """
    self.mean = None; self.covar = None;
    self.size = None; self.precision = None;
    self.invertible = False; self.mass = 1.0; 
    self.weight = 1.0
    self.decomposition = False
    self.quat = None; self.rotation = None; self.radii = None;
    self.data = None; self.dataMin = None; self.dataMax = None;

  def set_data(self,data,mass,weight = None,scale = 1.0):
    """
      The method to make ellipsoid from coordinate data. The
      coordinate data composed of location of N constituent data
      points in 3D. Each data-point is associated to a weight (mass).
      data: The data-matrix of size Nx3
      mass: Mass (weight) of each data-point
      weight: Weight of each data-point (OBSOLETE parameter)
      scale: Scale of standard-deviation applied
    """
    if weight is None: weight = mass      
    self.weight = weight

    # Given data and weight (reliability weights) calculate unbiased estimates
    self.size = len(data)    
    self.mass = sum(np.array(weight))

    # Calculate covariance matrix and mean
    D = np.array(data)
    W = np.array(weight) / self.mass
    mean = np.sum((D.T * W).T,axis=0)
    dd = ((D - mean).T * np.sqrt(W)).T
    covar = np.dot(dd.T, dd) / (1-np.sum(W*W))
    self.data = data
    self.mean = mean
    self.covar = covar

    # Check invertibility of covariance
    self.invertible = self.isInvertible()

    # Calculate inverse of covariance matrix
    self.precision = self.CalcPrecision()

    # Calculate standard deviations associated to covariance matrix
    self.sd = np.sqrt(np.diag(self.covar))

    # Sets up multivariate normal distribution class
    self.pnorm = multivariate_normal(mean=self.mean,cov=scale*self.covar,
                                     allow_singular=True)

    # Calculate Moment of Inertia matrix (now OBSOLETE)
    self.moi = np.trace(covar)*np.identity(3) - covar  

    # Set up data ranges  
    self.dataMin = np.min(self.data,axis=0)
    self.dataMax = np.max(self.data,axis=0)

    # Set up cutoff scale
    # See gmconvert source, Ellipsoid.c, function: InvCumurative_ChiSquare_DOF3
    chi2ppf = chi2(df=3); ppf_val = np.sqrt(chi2ppf.ppf(Ellipsoid.chi2CoverRatio))
    self.scale = scale * ppf_val

  def set_fromGDF(self,mean,covmat,W,scale=1.0,resample_point=1000,cutoff_pdf=None):
    """
      The method to set ellipsoid from Gaussian probability distribution (GDF)
      data. see class GMData.
      mean: Mean of GDF
      covmat: Covariance matrix of GDF
      W: Mixture weight of GDF. Used to resample data-points.
      scale: Scale of standard-deviation applied
      resample_point: Number of points to resample based on GDF
      cutoff_pdf: Resample point cutoff
    """
    # Sets data directly
    self.set_mean(mean); self.set_covar(covmat);

    # Set up total mass and weight of each data-point
    self.size = resample_point    
    self.weight = np.ones(self.size)    
    #self.mass = sum(np.array(self.weight)) # Bug: mass
    self.mass = W

    # Check invertibility of covariance matrix
    self.invertible = self.isInvertible()

    # Calculate inverse of covariance matrix
    self.precision = self.CalcPrecision()

    # Sets up multivariate normal distribution class
    self.pnorm = multivariate_normal(mean=self.mean,cov=scale*self.covar,
                                     seed=100,allow_singular=True)    

    # Calculate Moment of Inertia matrix (now OBSOLETE)
    self.moi = np.trace(covmat)*np.identity(3) - covmat   

    # Sample data according to GDF 
    self.sample_data(W,cutoff_pdf=cutoff_pdf)

    # Set up data ranges 
    self.dataMin = np.min(self.data,axis=0)
    self.dataMax = np.max(self.data,axis=0)

    # Set up cutoff scale
    chi2ppf = chi2(df=3); ppf_val = np.sqrt(chi2ppf.ppf(Ellipsoid.chi2CoverRatio))
    self.scale = scale * ppf_val

  def sample_data(self,W,cutoff_pdf=None):
    """
      Resampling data based on GDF
      W: Mixture weight of GDF
      cutoff_pdf: Resample point cutoff
    """
    if cutoff_pdf is None:
      self.data = self.pnorm.rvs(self.size)
    else:      
      self.data = np.zeros((self.size,3))
      for i in range(self.size):
        continue_sample = True
        while continue_sample:
          point = self.pnorm.rvs(1)
          pdf = self.pnorm.pdf(point) * W
          if pdf > cutoff_pdf:
            continue_sample = False
        self.data[i,:] = point[:]   
    
  def set_mean(self,mean):
    """
      Sets mean of ellipsoid
      mean: Input mean     
    """
    self.mean = mean

  def set_covar(self,covar):
    """
      Sets covariance and SD of ellipsoid
      covar: Input covariance
    """
    self.covar = covar
    self.sd = np.sqrt(np.diag(self.covar))

  def set_mass(self,mass):
    """
      Set total mass of ellipsoid
      mass: Mass of ellipsoid
    """
    self.mass = mass

  def set_size(self,size):
    """
      Set size of data in ellipsoid
      size: data-size
    """
    self.size = size

  def set_weight(self,weight):
    """
      Set ellipsoid data-point weight
      weight: Weight of each data-point in ellipsoid
    """
    self.weight = weight

  def set_pnorm(self):
    """
      Set multivariate normal distribution class
    """
    if self.covar is None or self.mean is None: return False
    self.pnorm = multivariate_normal(mean=self.mean,cov=self.covar)

  def isInvertible(self):
    """
      Checks invertibility of covariance matrix
      returns: True if invertible, else False
    """
    if self.covar is None or self.size is None: return False
    if self.size >= 3: 
      self.dcovar = np.linalg.det(self.covar)
      if np.abs(self.dcovar) <= Ellipsoid.tol:
        return False
      else:
        return True
    else:
      return False

  def CalcPrecision(self):
    """
      Calculate inverse of covariance matrix
      returns: inverse of covariance matrix
    """
    if not self.covar is None:      
      if self.invertible:
        precision = np.linalg.inv(self.covar)
      else:
        precision = np.linalg.pinv(self.covar,rcond = Ellipsoid.tol)
    else:
      precision = None
    return precision

  def decompose(self,mat,decreasing=False):
    """
      Performs eigenvalue decomposition of input matrix and sort result
      mat: Input matrix
      decreasing: Order of eigenvalues, and correspondingly eigenvectors
      returns: tuple of eigenvalues and eigenvectors
    """     
    w, v = np.linalg.eig(mat)
    if decreasing:
      windex = np.argsort(-w)
    else:
      windex = np.argsort(w)
    val = np.zeros(3); vec = np.zeros((3,3));
    for i,iw in enumerate(windex):
      val[i] =  w[iw]
      vec[:,i] = v[:,iw]
    return (val,vec)

  def do_decomposition(self,mat,decreasing=False):
    """
      Wrapper to the decompose method and performs checks. If check fails
      fix the decomposition and invert decomposition to recreate fixed matrix
      mat: Input matrix
      decreasing: Order of eigenvalues, and correspondingly eigenvectors
      returns: tuple of eigenvalues,eigenvectors and fixed matrix
    """
    s, rotation = self.decompose(mat,decreasing=decreasing)    
    rotation = self.fixChirality(rotation)
    radii = np.sqrt(s)
    rotation = rotation
    self.decomposition = True
    mat = invert_decomposition(rotation,radii)
    return (radii,rotation,mat)

  def fixChirality(self,r):
    """
      Fix Chirality related error in decomposition
      r: Input rotation matrix (orthogonal matrix)
      returns: fixed matrix
    """
    if np.linalg.det(r) < 0: r = -1.0 * r
    if sum(np.cross(r[:,0],r[:,1]) * r[:,2]) < 0.0: r[:,2] = -r[:,2]    
    return r

  def get_rotation_matrix(self):
    """
      Calculates rotation matrix
      returns: rotation matrix
    """
    try:
      rotation = self.rotation
    except:
      self.calculate_aux()
      rotation = self.rotation
    return rotation

  def get_radii(self,scaled=True):
    """
      Get radii from covariance matrix
      scaled: True, returns scaled radii, else, return unscaled radii
      returns: radii
    """
    if scaled:
      try:
        radii = self.scaled_radii
      except:
        self.calculate_aux()
        radii = self.scaled_radii
    else:           
      try:
        radii = self.radii
      except:
        self.calculate_aux()
        radii = self.radii
    return radii

  def calculate_aux(self):
    """
      Calculate auxilliary class variables from covariance matrix
      e.g. rotation, radii, and quaterion corresponding to rotation
    """
    self.radii, self.rotation, self.covar = self.do_decomposition(self.covar,decreasing=True)
    self.scaled_radii = self.radii * self.scale
    if self.is_rotation_valid(self.rotation):
      tr = Transformation(rot=self.rotation,trans=None,rotationType='m')
      self.quat = list(tr.change_rotation(fromtype='m',totype='q'))
    else: self.quat = [1.,0.,0.,0.]

  def set_radii(self,radii,scaled=True):
    """
      Sets radii
      radii: Input radii
      scaled: True, then sets scaled radii, else sets unscaled radii
    """
    if scaled:
      self.scaled_radii = radii
    else:
      self.radii = radii

  def is_rotation_valid(self,rm):
    """
      Checks validity of rotation matrix
      returns: True, if valid, else False
    """
    if not np.all(np.abs(np.dot(rm,rm.T) - np.identity(3)) < 1e-8):
      return False
    if not np.all(np.abs(np.dot(rm.T,rm) - np.identity(3)) < 1e-8):
      return False
    return True

  def get_rotation4D(self):
    """
      Calculate 4D rotation vector (quaternion)
      returs: Quaternion corresponding to covariance matrix
    """
    try:
      quat = self.quat
    except:
      self.calculate_aux()
      quat = self.quat
    return quat

  def set_rotation4D(self,q):
    """
      Sets quaternion
      q: Input quaternion
    """
    self.quat = q

class Mesh:
  """
    Mesh builds surface around each ellipsoid and using the surface model
    make a AFM-like projection. Also it has methods to despatch VRML file
    writer, and GMM-file writer. During sampling this class methods can be
    used to update existing mesh and recast projection.
  """
  DegToRadian = np.pi / 180.0  
  def __init__(self,opt=None,MeshPoints=None):   
    """
      The constructor taking Option dictionary. This dictionary should be
      build by reading input yaml file, and using OptionAFM class, given
      in AFM_Header.py. Also it is possible to set MeshPoints directly as
      optional argument.
      opt: Input option dictionary
      MeshPoints: An array containing surface points of the ellipsoid
                  (N x M x 3 shaped, where M = number of points per ellipsoid,
                                           N = number of ellipsoid)
    """
    self.set_opt_settings(opt)    
    if MeshPoints is None:
      self.meshPoints = None
    else:
      self.set_MeshPoints(MeshPoints)

  def set_opt_settings(self,opt):
    """
      This saves only required parameters from option dictionary as class 
      variable.
      opt: Input option dictionary      
    """
    self.nGaussian = 0; self.imageParam = {}; self.has_image_param = False;
    self.imageMag = opt["MESH"]["MAGNIFICATION"]
    self.ScaleRadii = opt["MESH"]["SCALERADII"]
    if "XGRID" in opt["MESH"] or "YGRID" in opt["MESH"] or "ZGRID" in opt["MESH"]:
      self.xbins = opt["MESH"]["XGRID"]
      self.ybins = opt["MESH"]["YGRID"]
      self.zbins = opt["MESH"]["ZGRID"]    
      self.nx = len(self.xbins); self.ny = len(self.ybins); self.nz = len(self.zbins);
      self.xmin = self.xbins[0]; self.xmax = self.xbins[-1];
      self.ymin = self.ybins[0]; self.ymax = self.ybins[-1];
      self.zmin = self.zbins[0]; self.zmax = self.zbins[-1];
      self.xstep = self.xbins[1] - self.xbins[0]
      self.ystep = self.ybins[1] - self.ybins[0]
      self.zstep = self.zbins[1] - self.zbins[0]
    else:
      self.nx = 1; self.ny = 1; self.nz = 1;
      self.xmin = 0; self.ymin = 0; self.zmin = 0;
      self.xmax = 0; self.ymax = 0; self.zmax = 0;
      self.xstep = 0; self.ystep = 0; self.zstep = 0; 
    self.noiseFactor = opt["MESH"]["ZNOISE"]["FACTOR"]
    self.noiseQvalue = opt["MESH"]["ZNOISE"]["PERCENTILE"]    

    # Initialize a sphere to be used later
    self.initSphere = self.InitializeSphere(opt["MESH"]["GRIDSIZE"])

  def InitializeSphere(self,stepSizeGrid):
    """
      Initialize a sphere or generating surface points on an unit sphere
      stepSizeGrid: Input grid size to generate surface points. Smaller
                    value generates more points on the sphere
      returns: set of vertices on the sphere surface
    """
    self.minTheta = 0.0; self.minPhi = 0.0;
    self.maxTheta = 360.0; self.maxPhi = 180.0;

    # Parameter calculation for sphere generation
    self.stepSizeGrid = stepSizeGrid;
    self.numThetaGrid = int((self.maxTheta - self.minTheta)/self.stepSizeGrid + 1)
    self.numPhiGrid = int((self.maxPhi - self.minPhi)/self.stepSizeGrid + 1)
    self.nvert = self.numThetaGrid*self.numPhiGrid

    # Generating grid
    tg = np.linspace(self.minTheta, self.maxTheta, self.numThetaGrid) * Mesh.DegToRadian
    pg = np.linspace(self.minPhi, self.maxPhi, self.numPhiGrid) * Mesh.DegToRadian

    # Calculating values on grid
    costg = np.cos(tg); sinpg = np.sin(pg);
    sintg = np.sin(tg); cospg = np.cos(pg);
    xv = np.outer(costg, sinpg).ravel()
    yv = np.outer(sintg, sinpg).ravel()
    zv = np.outer(np.ones_like(tg), cospg).ravel()

    # Pack values to vertices
    vertices = np.vstack([xv,yv,zv]).T

    # Place sphere (the vertex set) centered at origin (0,0,0)
    vertices = vertices - np.mean(vertices,axis=0)

    return vertices     

  def set_MeshPoints(self,MeshPoints):
    """
      Sets mesh-points
      MeshPoints: An array containing surface points of the ellipsoid
                  (N x M x 3 shaped, where M = number of points per ellipsoid,
                                           N = number of ellipsoid)
    """
    self.meshPoints = MeshPoints

  def buildMesh(self,ellipsoids):
    """
      The mesh generator using unit sphere (already initialized in class) 
      and input ellipsoids.
      ellipsoids: Input set of ellipsoids. An array of N ellipsoids, where
                  each ellipsoid is of class Ellipsoid
    """
    # Initialize variables for generation
    self.nGaussian = len(ellipsoids)
    self.meshPoints = np.empty([self.nGaussian,self.nvert,3])
    self.meshPointIndex = np.empty([self.nGaussian,self.nvert,3])
    self.radii = [None] * self.nGaussian; 
    self.scaled_radii = [None] * self.nGaussian;
    self.scaling_factor = [None] * self.nGaussian;
    self.rotation = [None] * self.nGaussian;
    self.centers = [None] * self.nGaussian; self.ScalledVertices = [None] * self.nGaussian;  
    self.covar = [None] * self.nGaussian;
    self.ellipWeight = [None] * self.nGaussian;

    # Loop over all ellipsoids taking ellipsoid parameters and 
    # generating vertices accordingly
    for ie,e in enumerate(ellipsoids):

      self.scaled_radii[ie] = e.get_radii(scaled=True); self.radii[ie] = e.get_radii(scaled=False);
      self.rotation[ie] = e.get_rotation_matrix()
      self.centers[ie] = e.mean; self.covar[ie] = e.covar;
      self.scaling_factor[ie] = e.scale
      self.ellipWeight[ie] = e.mass

      # Scale unit sphere to correct size and save scaled vertices per ellipsoid
      vertices = np.multiply(self.initSphere,self.scaled_radii[ie])      
      self.ScalledVertices[ie] = vertices.T
      print("debug:: shape = %d %d" % (self.ScalledVertices[ie].shape[0],self.ScalledVertices[ie].shape[1]))
      # Transform scaled vertices
      self.meshPoints[ie,:,:] = self.rotation[ie].dot(self.ScalledVertices[ie]).T + self.centers[ie] 
      self.meshPointIndex[ie,:,:] = ie

    tw = np.sum(self.ellipWeight)
    for ie in range(len(ellipsoids)):
      self.ellipWeight[ie] = self.ellipWeight[ie]/tw
    
    # Do a simple check on mesh and reshape the mesh-point array
    if self.checkMesh(self.radii):
      self.meshPoints = np.reshape(self.meshPoints,[self.nGaussian*self.nvert,3])
      self.meshPointIndex = np.reshape(self.meshPointIndex,[self.nGaussian*self.nvert,3])
      self.success = True
    else:
      self.success = False  

    # Model noise from saved parameters
    if self.noiseFactor is None:
      self.noiseFactor = np.percentile(self.meshPoints[:,2],self.noiseQvalue)
    self.bg = float(self.noiseFactor)*np.random.uniform(low=0.0,high=1.0,size=(self.nx,self.ny))

    # Initialize matrix to be projected
    self.zstat = np.full((self.nx,self.ny),self.zmin) + self.bg 
    self.zindex = np.zeros((self.nx,self.ny)); self.zindex.fill(None);

  def checkMesh(self,list_radii):
    """
      A simple check on mesh for radii values are positive
      list_radii: An array of radii values
    """    
    for radii in list_radii:
      for val in radii:
        if val <= 0 or val is None:
          return False
    return True

  def get_scaled_vertices(self,new_radii=False,radii_values=None):
    """
      This return scaled vertices using already used radii or use
      new_radii
      new_radii: Whether scaled vertices will be updated or not
      radii_values: new radii values (scaled radii)
      returns: updated scaled vertices from initial unit sphere or original
    """   
    if new_radii:
      for i in radii_values:
        rad = radii_values[i]
        if rad is None:
          continue
        vertices = np.multiply(self.initSphere,rad)      
        self.ScalledVertices[i] = vertices.T
        self.scaled_radii[i] = rad
        self.radii[i] = rad / self.scaling_factor[i]
    return self.ScalledVertices 

  def updateMesh(self,trs,new_radii=False,radii_values=None): 
    """
      This will update current mesh using new transformation. A mesh is formed
      from a set of scaled-vertices surface that define ellipsoids but not 
      transformed correctly. The ellipsoid parameters rotation and center are
      used to transform them correctly. Here, in stead of rotation and center 
      are retrieved from ellipsoids, they are given as class Transformation.
      trs: An array of Transformation class instances
      new_radii: Whether scaled vertices will be updated or not
      radii_values: new radii values (scaled radii)
    """          
    # Backup current mesh
    self.backup()    

    print("updating mesh")
    self.updated = True
    print("debug:: shape = %d %d" % (self.ScalledVertices[0].shape[0],self.ScalledVertices[0].shape[1]))
    ScalledVertices = self.get_scaled_vertices(new_radii=new_radii,radii_values=radii_values)

    # Initialize mesh-points array    
    meshPoints = np.full((self.nGaussian,self.nvert,3),0.0)
    meshPointIndex = np.empty([self.nGaussian,self.nvert,3])
    mp_prev = np.reshape(self.meshPoints,(self.nGaussian,self.nvert,3))    

    # Loop over each transformation instances and apply transformation to scaled
    # vertices for each ellipsoid already formed.
    for it in trs:
      tr = trs[it]
      meshPointIndex[it,:,:] = it
      if tr is None:
        meshPoints[it,:,:] = mp_prev[it,:,:]        
        print("Keep %d intact" % it)
        continue
      print("update %d-th ellipsoid,mesh" % it)
      curr_rotation = tr.change_rotation(fromtype='q',totype='m')
      curr_trans = tr.get_translation()
      meshPoints[it,:,:] = curr_rotation.dot(ScalledVertices[it]).T + curr_trans      
      self.centers[it] = curr_trans; self.rotation[it] = curr_rotation;
      # Using rotation and radii back-calculate covariance matrix
      self.covar[it] = invert_decomposition(curr_rotation,self.radii[it])

    # Save to class variables
    self.meshPoints = np.reshape(meshPoints,[self.nGaussian*self.nvert,3])
    self.meshPointIndex = np.reshape(meshPointIndex,[self.nGaussian*self.nvert,3])

    # Re-initialize matrix to be projected
    self.zstat = np.full((self.nx,self.ny),self.zmin) + self.bg
    self.zindex = np.zeros((self.nx,self.ny)); self.zindex.fill(None);
      
  def backup(self):
    """
      This is for backing up existing mesh class variables so that they can
      be restored if required.
    """
    print("backing-up mesh")    
    self.meshPoints_save = self.meshPoints[:,:]
    self.meshPointIndex_save = self.meshPointIndex[:,:]
    self.radii_save = self.radii[:]
    self.scaled_radii_save = self.scaled_radii[:]
    self.centers_save = self.centers[:]
    self.rotation_save = self.rotation[:]
    self.zstat_save = self.zstat[:,:]
    self.zindex_save = self.zindex[:,:]
    self.ScalledVertices_save = [None] * len(self.ScalledVertices)
    for i in range(len(self.ScalledVertices)):
      self.ScalledVertices_save[i] = np.zeros(self.ScalledVertices[i].shape)
      self.ScalledVertices_save[i][:,:] = self.ScalledVertices[i][:,:]

  def restore_backup(self):
    """
      This restores mesh (if it is updated) to the mesh prior to update
    """
    if not hasattr(self,"updated"): self.updated = False
    if self.updated:
      print("restoring mesh")
      self.meshPoints = self.meshPoints_save[:,:]
      self.meshPointIndex = self.meshPointIndex_save[:,:]
      self.radii = self.radii_save[:]
      self.scaled_radii = self.scaled_radii_save[:]
      self.centers = self.centers_save[:]
      self.rotation = self.rotation_save[:]
      self.zstat = self.zstat_save[:,:]      
      self.zindex = self.zindex_save[:,:]
      self.ScalledVertices = [None] * len(self.ScalledVertices_save)
      for i in range(len(self.ScalledVertices_save)):
        self.ScalledVertices[i] = np.zeros(self.ScalledVertices_save[i].shape)
        self.ScalledVertices[i][:,:] = self.ScalledVertices_save[i][:,:]      
    self.updated = False  

  def _do_help_projection(self,i,j):
    """
      A helper method for projection building. 
      Not for calling directly.      
    """
    k = np.intersect1d(self.wx[i],self.wy[j])
    if len(k) > 0:
      self.zz[i,j] = True
      self.xyindex[i*self.nuidy+j] = k      
  
  def buildProjection(self,num_thread = 8,timing=False):
    """
      The method to build AFM-like projection from mesh-points
      num_thread: Number of threads to use
      timing: Whether to use method with timing reports
    """
    print("building projection")

    if timing:
      # Despatch method with timing reports
      self.buildProjectionWithTiming(num_thread = num_thread)
    else:
      # Less verbose
      if not self.success: return

      # Initializations
      MP = self.meshPoints; MPi = self.meshPointIndex;
      xmin = self.xmin; xmax = self.xmax;
      ymin = self.ymin; ymax = self.ymax;
      
      # Restrict to xy-grids
      mpb = (MP[:,0] >= xmin) & (MP[:,0] <= xmax) & (MP[:,1] >= ymin) & (MP[:,1] <= ymax)
      MP = MP[mpb]; MPi = MPi[mpb];

      # Digitize each vertex value to a bin
      idX = np.digitize(MP[:,0],self.xbins) - 1
      idY = np.digitize(MP[:,1],self.ybins) - 1

      # Work with obtained bins only
      uidX = np.unique(idX)
      uidY = np.unique(idY)
      wx = []; wy = [];

      # Indexing unique bins X
      for ix in uidX: wx.append(np.where(ix == idX))
      # Indexing unique bins Y
      for iy in uidY: wy.append(np.where(iy == idY))

      # Saving parameters for thread call
      self.nuidx = len(uidX); self.nuidy = len(uidY); 
      self.xyindex = [None] * self.nuidx * self.nuidy; 
      self.zz = np.zeros(self.zstat.shape,dtype=bool)
      self.wx = wx; self.wy = wy;
      thread = [None] * self.nuidx * self.nuidy; k = 0;

      # Threaded loop to fill xyindex and zz
      for i in range(self.nuidx):
        for j in range(self.nuidy):
          thread[k] = Thread(target=self._do_help_projection, args=(i,j))
          k += 1    
      for i in range(0,len(thread),num_thread):
        thr = thread[i:i+num_thread]
        for th in thr: th.start()
        for th in thr: th.join()
      
      # Loop to take only representative vertex (with max(z)) from each grid
      for i,ix in enumerate(uidX):
        for j,iy in enumerate(uidY):
          k = self.xyindex[i*self.nuidy+j]
          if self.zz[i,j]: 
            self.zstat[ix,iy] = np.max(MP[k,2])
            self.zindex[ix,iy] = int(MPi[k[np.argmax(MP[k,2])],2])

  def buildProjectionWithTiming(self,num_thread = 8):
    """
      The method to build AFM-like projection from mesh-points
      num_thread: Number of threads to use
    """
    # Initializations
    start_clock = datetime.now()                
    MP = self.meshPoints; MPi = self.meshPointIndex;
    xmin = self.xmin; xmax = self.xmax; 
    ymin = self.ymin; ymax = self.ymax;

    # Restrict to xy-grids
    mpb = (MP[:,0] >= xmin) & (MP[:,0] <= xmax) & (MP[:,1] >= ymin) & (MP[:,1] <= ymax)
    MP = MP[mpb]; MPi = MPi[mpb]; 

    # Digitize each vertex value to a bin
    idX = np.digitize(MP[:,0],self.xbins) - 1
    idY = np.digitize(MP[:,1],self.ybins) - 1
    end_clock = datetime.now(); elap_clock = end_clock - start_clock; 
    elap_ms = int(elap_clock.total_seconds() * 1000000); 
    print("DEBUG: Time taken to build projection 1: %d" % elap_ms);

    # Work with obtained bins only
    start_clock = datetime.now();
    uidX = np.unique(idX)
    uidY = np.unique(idY)
    wx = []; wy = [];    
    end_clock = datetime.now(); elap_clock = end_clock - start_clock; 
    elap_ms = int(elap_clock.total_seconds() * 1000000); 
    print("DEBUG: Time taken to build projection 2: %d" % elap_ms); 

    # Indexing unique bins X
    start_clock = datetime.now();
    for ix in uidX: wx.append(np.where(ix == idX))
    end_clock = datetime.now(); elap_clock = end_clock - start_clock; 
    elap_ms = int(elap_clock.total_seconds() * 1000000); 
    print("DEBUG: Time taken to build projection 3: %d" % elap_ms); 

    # Indexing unique bins Y
    start_clock = datetime.now();
    for iy in uidY: wy.append(np.where(iy == idY))
    end_clock = datetime.now(); elap_clock = end_clock - start_clock; 
    elap_ms = int(elap_clock.total_seconds() * 1000000); 
    print("DEBUG: Time taken to build projection 4: %d" % elap_ms); 

    # Saving parameters for thread call
    start_clock = datetime.now();
    self.nuidx = len(uidX); self.nuidy = len(uidY); 
    self.xyindex = [None] * self.nuidx * self.nuidy; 
    self.zz = np.zeros(self.zstat.shape,dtype=bool)
    self.wx = wx; self.wy = wy;
    thread = [None] * self.nuidx * self.nuidy; k = 0;

    # Threaded loop to fill xyindex and zz
    for i in range(self.nuidx):
      for j in range(self.nuidy):
        thread[k] = Thread(target=self._do_help_projection, args=(i,j))
        k += 1    
    for i in range(0,len(thread),num_thread):
      thr = thread[i:i+num_thread]
      for th in thr: th.start()
      for th in thr: th.join()
    end_clock = datetime.now(); elap_clock = end_clock - start_clock; 
    elap_ms = int(elap_clock.total_seconds() * 1000000) / num_thread; 
    print("DEBUG: Time taken to build projection 5: %d" % elap_ms);

    # Loop to take only representative vertex (with max(z)) from each grid
    start_clock = datetime.now();
    for i,ix in enumerate(uidX):
      for j,iy in enumerate(uidY):
        k = self.xyindex[i*self.nuidy+j]
        if self.zz[i,j]: 
          self.zstat[ix,iy] = np.max(MP[k,2])
          self.zindex[ix,iy] = int(MPi[k[np.argmax(MP[k,2])],2])
    end_clock = datetime.now()
    elap_clock = end_clock - start_clock; elap_ms = int(elap_clock.total_seconds() * 1000000)
    print("DEBUG: Time taken to build projection 6: %d" % elap_ms)  
  
  def showProjection(self,rotate=False):
    """
      Plot the projected matrix
      rotate: Input degree to rotate the matrix for plotting
    """
    if not self.success: return
    vmin = self.zmin; vmax = self.zmax; mat = self.zstat;
    if rotate: plt.imshow(ndimage.rotate(mat,-90),cmap=plt.cm.gray_r,vmin=vmin,vmax=vmax)
    else:      plt.imshow(mat,cmap=plt.cm.gray_r,vmin=vmin,vmax=vmax)
    plt.show(); plt.close();

  def get_figure_dimension(self,dpi=None):
    """
      Get projection figure dimensions
      dpi: Input DPI
      returns: width, height and dpi of the figure
    """
    mag = self.imageMag
    if dpi is None: dpi = plt.rcParams['figure.dpi']     
    num_bin_horiz = self.ny
    num_bin_verti = self.nx
    width  = (num_bin_horiz * mag) / dpi
    height = (num_bin_verti * mag) / dpi
    return (width,height,dpi)

  def apply_image_rotation(self,mat,rotate=False):
    """
      Rotate (like an image) the input matrix
      mat: Input matrix (2D)
      rotate: Input degree or False for no rotation
      returns: Rotated matrix
    """
    if rotate: return ndimage.rotate(mat,-90)
    else:      return mat

  def apply_magnification(self,mat,param,proportion='constraint',direction=None):
    """
      Magnify the matrix (like an image)
      mat: Input matrix
      param: Image parameter dictionary, with keys 'width', 'height' and 'dpi'
      proportion: How to resize the matrix? Can be 'constraint' or 'constant'
                  If 'constraint' matrix (as image) is resized with resampling
                  method bi-linear in both the directions. 
                  If 'constant', then argument direction is used, the matrix
                  (as image) is resized along the direction
      direction: either 'width' or 'height'
      returns: magnified matrix
    """
    if self.imageMag == 1.0: return mat

    # Make image from matrix
    im = Image.fromarray(mat)

    # Get parameters
    w = param['width']; h = param['height']; dpi = param['dpi']
    wpixel = int(w * dpi)
    hpixel = int(h * dpi)

    # Magnify
    if proportion == 'constraint':
      # Magnify along both the direction similarly
      im2 = im.resize((wpixel,hpixel),resample=Image.BILINEAR)
    elif proportion == 'constant':
      # Magnify along a given direction width or height
      if direction is None:
        raise Exception('invalid direction')
      else:
        if direction == 'width':
          percent = (wpixel/float(im.size[0]))
          hpixel = int((float(im.size[1])*float(percent)))
          im2 = im.resize((wpixel,hpixel), resample=Image.BILINEAR)
        elif direction == 'height':
          percent = (hpixel/float(im.size[1]))
          wpixel = int((float(im.size[0])*float(percent)))
          im2 = im.resize((wpixel,hpixel), resample=Image.BILINEAR)
        else:
          raise Exception('invalid direction')
    else:
      raise Exception('invalid proportion')
    
    # Recreate matrix from image and return
    return np.asarray(im2.getdata(),dtype=np.float64).reshape((im2.size[1],im2.size[0]))

  def saveProjection(self,filename=None,rotate=False,dpi=None):
    """
      This saves the projected matrix to a filename and/or a internal class 
      variable MatrixProjected, which can be accessed to get the matrix
      filename: The output file to save the plotted matrix, or None to not
                save the matrix but save it in MatrixProjected. filename can
                be a buffer.
      rotate: Rotation in degrees to apply to the projection or False for no
              rotation
      dpi: The input DPI
      returns: The filename
    """
    if not self.success: return

    # Rotate if required
    mat = self.apply_image_rotation(self.zstat,rotate=rotate)

    if not self.has_image_param:
      # Get image parameters if not set already
      w,h,d = self.get_figure_dimension(dpi=dpi)
      self.imageParam['width'] = w
      self.imageParam['height']= h
      self.imageParam['dpi']   = d
      self.has_image_param = True
    else:
      w = self.imageParam['width']
      h = self.imageParam['height']
      d = self.imageParam['dpi']
    # Print some information about the figure
    print("Figure dimension (inches) [w x h]: %.3f %.3f" % (w,h))
    print("Figure dpi: %d" % d)
    
    if filename is None:      
      # If filename is not given save the magnified matrix to the
      # MatrixProjected variable
      self.MatrixProjected = self.apply_magnification(mat,self.imageParam)      

      # Print information corresponding to original and modified matrix
      s0 = mat.shape[0]; s1 = mat.shape[1]
      print("Original Projected Matrix dimension: %d x %d" % (s0,s1))
      s0 = self.MatrixProjected.shape[0]; s1 = self.MatrixProjected.shape[1];
      print("Magnified Projected Matrix dimension: %d x %d" % (s0,s1))
      return filename

    else:

      # If filename is given setup matplotlib figure
      figure = plt.figure(frameon=False)
      figure.set_size_inches([self.imageParam['width'],self.imageParam['height']],forward=True)         
      print("confirm: %f %f" % (self.imageParam['width'],self.imageParam['height']))
      plt.subplots_adjust(left=0,bottom=0,top=1,right=1)
      ax = plt.Axes(figure,[0.0,0.0,1.0,1.0])
      ax.set_axis_off()
      figure.add_axes(ax)
      bbox = figure.get_window_extent().transformed(figure.dpi_scale_trans.inverted())
      print("Content width x height: %.3f x %.3f" % (bbox.width,bbox.height))

      # make colored matrix and put it for printing
      mat_up_dark = 255 - (np.digitize(mat,self.zbins) - 1)
      #ax.imshow(mat,cmap=plt.cm.gray_r,aspect='auto',interpolation='none',vmin=self.zmin,vmax=self.zmax)
      ax.imshow(mat_up_dark,cmap=plt.cm.gray,aspect='auto',interpolation='none',vmin=0,vmax=255)
      ax.margins(0)         

      # Complete the figure
      figure.savefig(filename,dpi=dpi,pad_inches=0,bbox_inches=bbox)
      plt.close(figure)
      print("Original Projected Matrix dimension: %d x %d" % (mat.shape[0],mat.shape[1]))

      # If filename is a string save magnified matrix to the variable MatrixProjected
      # else return filename as a buffered io
      if isinstance(filename,str):

        self.MatrixProjected = self.apply_magnification(mat,self.imageParam)                    
        s0 = self.MatrixProjected.shape[0]; s1 = self.MatrixProjected.shape[1];
        print("Magnified Projected Matrix dimension: %d x %d" % (s0,s1))

        # Check the image dimension
        check = self.checkImage(filename)
        if check: 
          if check == 2: print("Image check: ............... [pass]")
          else:          print("Image check: ............... [pass with warning]")
        else: print("Image check: ............... [fail]")
        self.imageCheck = check
      else:
        return filename

  def get_projected_matrix(self,invert=False):
    """
      Returns projected matrix (AFM-like projection)
      invert: Whether to invert (return negative of the projected matrix)
              or not
      returns: projected matrix
    """
    if invert:      
      return -self.MatrixProjected
    else:
      return self.MatrixProjected

  def captureProjectionDirect(self):
    """
      To capture the matrix projected which would be printed to a figure
      withing printing it to a device
      returns: captured (photographed) matrix
    """
    vmin = self.zmin; vmax = self.zmax;
    Norm = plt.Normalize(vmin=vmin,vmax=vmax)
    pmat_rgba = plt.cm.gray_r(Norm(self.MatrixProjected),bytes=True)
    return np.dot(pmat_rgba[...,:3], [0.299, 0.587, 0.114])

  def checkImage(self,filename):
    """
      This checks size of the image given in filename
      filename: Input filename containing the image
      returns: integer (pass/fail test mark). 
               2 = all passed
               1 = one passed
               0 = failed
    """
    if self.checkSizeImage(filename):
      print("\tProduced Image Size ............... [pass]") 
      if self.checkReliable(filename): 
        print("\tSelf-Correlation ............... [pass]") 
        return 2
      else:
        print("\tSelf-Correlation ............... [fail]") 
        return 1
    else: 
      print("\tProduced Image Size ............... [fail]") 
      return 0

  def readImageAsMatrix(self,filename):
    """
      This reads the image in filename as a matrix
      filename: Input filename containing the image
      returns: Matrix representation of the image (RGB)
    """
    im = Image.open(filename)
    im = im.convert('L')
    mat = np.asarray(im.getdata(),dtype=np.float64).reshape((im.size[1],im.size[0]))
    return mat

  def checkSizeImage(self,filename):
    """
      Checks image size given in filename
      Note: due to roundoff error matplotlib saved figure dimenstion may not
            be identical to the matrix saved in the class
      filename: Input filename containing the image
      returns: True if image size as output in disk is identical to the matrix
               formed in the variable. If False, roundoff-error due to dpi 
               settings is present.
    """
    mat = self.readImageAsMatrix(filename)
    num_bin_verti, num_bin_horiz = mat.shape        
    if num_bin_horiz != self.ny*self.imageMag or num_bin_verti != self.nx*self.imageMag:
      print("\tRound-off error due to dpi settings detected")
      print("\tDisk Image dimension: w x h: %d x %d" % (num_bin_horiz,num_bin_verti))
      print("\tRequested Image dimension: w x h: %d x %d" % (self.ny,self.nx))
      return False
    else: return True

  def checkReliable(self,filename):
    """
      Check reliability of the color mapping made during saving.
      The image in the filename and image made internally should be nearly
      identical. Otherwise, self-correlation is less. However, such 
      self-correlation can be low due to error in output.
      filename: Input filename containing the image
      returns: True if self-correlation is nearly 1.0. False, lower 
               self-correlation, change is more than 0.1%.
    """
    pmat = self.captureProjectionDirect()
    imat = self.readImageAsMatrix(filename)
    cc = pearsonr(pmat.ravel(),imat.ravel())[0]       
    if (1.0 - cc) < 0.001: return True
    else:
      print("\tself correlation: %.6e" % cc)
      print("\tself correlation out of limit detected")
      return False

  def buildVRML(self,col=None):
    """
      Make VRML data from class variables
      Note: VRML files ellipsoid should use scaled radii, not unscaled radii.
      col: Input color for each ellipsoid
    """
    if not self.success: return
    self.vrmlData = VRML(center=self.centers,rotation=self.rotation,radii=self.scaled_radii,
                         covmat=self.covar,col=col,weight=np.array(self.ellipWeight))

  def buildGMData(self):
    """
      Make GMM-data from class variables
      Note: GMM files use covariance matrix formed from unscaled radii      
    """
    if not self.success: return
    if hasattr(self,"vrmlData"):
      # If vrmlData is present in class instance, use that to make GMM data
      self.gmData = GMData(fromVRMLData=True,VRMLData=self.vrmlData)
    else:
      # Make GMM-data from ellipsoid weights, centers, rotations and radii
      # w = np.repeat(1.0 / self.nGaussian,self.nGaussian)   
      w = np.array(self.ellipWeight)      
      self.gmData = GMData(weight=w,center=self.centers,rotation=self.rotation,radii=self.radii)

  def writeVRML(self,outFileName):
    """
      Write VRML data to file
      outFileName: Output filename
    """
    if not self.success: return
    data = self.vrmlData
    data.write(outFileName)

  def writeGMData(self,outFileName):
    """
      Write GMM-data to file
      outFileName: Output filename
    """
    if not self.success: return
    data = self.gmData
    data.write(outFileName)

  def plotMesh(self):
    """
      This plots mesh as a scatterplot. This can be very time-consuming and 
      should be used only during debugging
    """
    if not self.success: return
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(self.meshPoints[:,0],self.meshPoints[:,1],self.meshPoints[:,2],color='b',alpha=0.2)
    plt.show()    

class VRML:
  """
    A class to read, create, and write ellipsoids as VRML file
  """
  def __init__(self,center=None,rotation=None,radii=None,col=None,weight=None,scaleRadii=1.0,covmat=None,
               fromGMData=False,GMData=None,fromFile=False,filename=None):    
    """
      The constructor which either read VRML data from a file, or make VRML 
      data from inputs

      This class can be instantiated in three ways.
      1) vrml = VRML(center=Array_Of_Centers,rotation=Array_Of_Rotations,
                     radii=Array_Of_Radii,covmat=Array_Of_covmats)
      2) vrml = VRML(fromFile=True,filename=Input_VRML_File)
      3) vrml = VRML(fromGMData=True,GMData=GMData_Class_Instance)

      center: Array of centers of ellipsoids
      rotation: Array of rotation matrices of ellipsoids
      radii: Array of radii of ellipsoids
      col: Color of VRML file to be written
      scaleRadii: How the radii should be scaled
      covmat: Array of covariance matrices of ellipsoids
      fromGMData: Whether to make VRML data from GMM-data
      GMData: An instance of class GMData
      fromFile: Whether to read VRML data from a file
      filename: The VRML file from which VRML data will be read
    """
    # Set vrml header
    self.header = '#VRML V2.0 utf8'    

    # set input scale for radii
    self.scaleRadii = scaleRadii;

    # set color per ellipsoid, default to blue
    if col is None: self.col = [0.,1.,0.]
    else: self.col = col

    if fromFile and filename is not None:
      # Read from filename
      information = self.read(filename,verbose=0)
      self.translation = information[0]; self.scale = information[1]; 
      self.rotation = information[2]; self.rotationMat = information[3];
      self.covmat = information[4]; self.elmName = information[5];
      self.weight = information[6];
    else:
      if fromGMData and GMData is not None:     
        # Set from GMData
        self.translation = []; self.covmat = []; self.rotation = []; 
        self.scale = []; self.rotationMat = []; self.weight = [];
        chi2ppf = chi2(df=3); ppf_val = np.sqrt(chi2ppf.ppf(Ellipsoid.chi2CoverRatio))
        scaleFactor = self.scaleRadii * ppf_val
        for i in range(len(GMData.covmat)):
          m = GMData.center[i]
          self.weight.append(GMData.W[i])
          s,r = np.linalg.eig(GMData.covmat[i])
          r,s = sort_decomposition(r,s,decreasing=True); 
          r = checkRotation(r); s = np.sqrt(s);
          covmat_regain = invert_decomposition(r,s)
          namecov, minicov = self.getcovMinimal(covmat_regain)          
          tr = Transformation(trans=m,rot=r,rotationType='m')
          AxisAngle = tr.QuaternionToAxisAngle(tr.change_rotation(fromtype='m',totype='q'))
          self.translation.append(m); self.scale.append(s * scaleFactor);
          self.rotation.append(AxisAngle); self.rotationMat.append(r);
          self.covmat.append(minicov); self.elmName = namecov;          
      else:
        # set from input center, rotation, radii, covmat
        if (center is None) or (rotation is None) or (radii is None) or (covmat is None):
          raise Exception("invalid initialization")            
        self.translation = []; self.covmat = []; self.rotation = []; 
        self.rotationMat = []; self.scale = [];         
        if weight is None:
          self.weight = 1.0 / len(center)
        else: self.weight = weight
        for i in range(len(center)):          
          namecov, minicov = self.getcovMinimal(covmat[i]);
          self.covmat.append(minicov); self.elmName = namecov;
          self.rotationMat.append(rotation[i])
          tr = Transformation(trans=center[i],rot=rotation[i],rotationType='m')
          AxisAngle = tr.QuaternionToAxisAngle(tr.change_rotation(fromtype='m',totype='q'))
          self.translation.append(center[i])
          self.scale.append(radii[i])
          self.rotation.append(AxisAngle)
  
  def getcovMinimal(self,cmat):
    """
      Return only upper-triangular part of covmat (3x3)
      cmat: Input covariance matrix (3x3)
      returns: tuple of Array_of_element_name and Array_of_matrix_values
    """
    minicov = [cmat[0,0],cmat[0,1],cmat[0,2],cmat[1,1],cmat[1,2],cmat[2,2]]
    namecov = ["xx","xy","xz","yy","yz","zz"]
    return namecov,minicov

  def read(self,filename,verbose=0):
    """
      Reading VRML file
      filename: Input filename
      returns: tuple of center,radii,rotation,rotationMat,covmat,namecov,weight
               where, center = array of centers
                      radii = array of radii (scaled)
                      rotation = array of rotation (axis-angle)
                      rotationMat = array of rotation matrix (3x3)
                      covmat = array of covariance matrix (minimal)
                      namecov = indices of covariance matrix returned
                      weight = weight of Gaussian
    """
    if exists(filename):
      # chi2CoverRatio from Ellipsoid class is used to calculate cutoff. 
      # In this case, scaleRadii input to constructor is used to scale cutoff
      # obtained from chi2CoverRatio. This scaleRadii should be the same 
      # scaleRadii used to instantiate class Ellipsoid (if done before).
      chi2ppf = chi2(df=3); ppf_val = np.sqrt(chi2ppf.ppf(Ellipsoid.chi2CoverRatio))
      scaleFactor = self.scaleRadii * ppf_val
      center = []; radii = []; rotation = []; rotationMat = []; covmat = [];      
      # Read filename
      # Note: vrml file contains 'scale' which are scaled radii. But we need
      #       unscaled radii to make covariance matrix.
      Lines = []   
      with open(filename) as fp:
        for line in fp:
          line = line.strip()
          if not line[0:1] == "#":
            if line[0:9] == "Transform": line = " ".join([" ",line])              
            Lines.append(line)      
      Line1 = " ".join(Lines)      
      tpiter = re.finditer("Transform",Line1)
      indices = [m.start(0) for m in tpiter]
      n = len(indices)
      weight = np.repeat(1.0/n,n)
      with open(filename) as fp:
        for line in fp:
          line = line.strip()
          if line[0:6] == "#Block":
            block_words = line.split()
            index_block = int(block_words[1])
            if len(block_words) >=3:
              weight[index_block] = float(block_words[2])
      for i in range(0,n-1):
        tblock = Line1[indices[i]:indices[i+1]]
        m = np.asarray(self._read_helper_parse_tblock(tblock,"translation",3))
        center.append(m)
        radii_scaled = np.asarray(self._read_helper_parse_tblock(tblock,"scale",3))        
        AxisAngle = self._read_helper_parse_tblock(tblock,"rotation",4)
        rotation.append(AxisAngle); radii.append(radii_scaled);
        radii_unscaled = radii_scaled / scaleFactor
        quat = self.convert2q(AxisAngle)
        tr = Transformation(trans=m,rot=quat,rotationType='q')
        rmat = tr.change_rotation(fromtype='q',totype='m')
        rotationMat.append(rmat)                                 
        cv = invert_decomposition(rmat,radii_unscaled);        
        namecov, minicov = self.getcovMinimal(cv);
        covmat.append(minicov)
      tblock = Line1[indices[n-1]:len(Line1)]
      m = np.asarray(self._read_helper_parse_tblock(tblock,"translation",3))
      center.append(m)
      radii_scaled = np.asarray(self._read_helper_parse_tblock(tblock,"scale",3))            
      AxisAngle = self._read_helper_parse_tblock(tblock,"rotation",4)
      rotation.append(AxisAngle); radii.append(radii_scaled);
      radii_unscaled = radii_scaled/scaleFactor
      quat = self.convert2q(AxisAngle)
      tr = Transformation(trans=m,rot=quat,rotationType='q')
      rmat = tr.change_rotation(fromtype='q',totype='m')
      rotationMat.append(rmat)     
      cv = invert_decomposition(rmat,radii_unscaled);
      namecov, minicov = self.getcovMinimal(cv);
      covmat.append(minicov)
      return center, radii, rotation, rotationMat, covmat, namecov,weight
    else:
      return None, None, None, None, None, None

  def _read_helper_parse_tblock(self,tblock,string,n):
    """
      Helper function for parsing
    """
    pt = re.search(string,tblock).start() + len(string)
    bl = tblock[pt:len(tblock)].split()
    val = []
    for i in range(n):
      val.append(float(bl[i]))
    return val

  def convert2q(self,rot):
    """
      Convert rotational vector from axis-angle to quaternion
      rot: Rotational vector in axis-angle format
      returns: Rotational quaternion
    """
    ax = rot[0]; ay = rot[1]; az = rot[2]; angle = rot[3]
    qx = ax * math.sin(angle/2.0)
    qy = ay * math.sin(angle/2.0)
    qz = az * math.sin(angle/2.0)
    qw = math.cos(angle/2.0) 
    return [qw,qx,qy,qz]

  def write(self,outFileName):
    """
      Write vrml data to a output file
      outFileName: Output filename
    """
    with open(outFileName,"w") as File:
      File.write(self.header + "\n")
      redCh = self.col[0]
      bluCh = self.col[1]
      greCh = self.col[2]
      for i,tr in enumerate(self.translation):
        s = self.scale[i]; rot = self.rotation[i];
        celm = self.covmat[i]
        File.write("#Block " + str(i))
        File.write(" %.5f\nTransform{\n" % self.weight[i])        
        File.write("  translation ")
        for j in range(3): File.write(" %.5f " % (tr[j]))
        File.write("\n")
        File.write("  scale ")
        for j in range(3): File.write(" %.5f " % (s[j]))
        File.write("\n")
        File.write("  rotation ")
        for j in range(4): File.write(" %.5f " % (rot[j]))
        File.write("\n")
        File.write("# covar ")
        for j in range(6): File.write(" %2s %.5f" % (self.elmName[j],celm[j]))
        File.write("\n")
        File.write("# rotationMat\n")
        for j in range(3):
          File.write("#");
          for k in range(3):
            File.write(" [%d,%d] %.5f" % (j,k,self.rotationMat[i][j,k]))
          File.write("\n")         
        File.write("  children[\n")
        File.write("    Shape{\n")
        File.write("      appearance Appearance{\n")
        File.write("        material Material{\n")
        File.write("          diffuseColor %.7f %.7f %.7f\n" % (redCh,bluCh,greCh))        
        File.write("          transparency 0.500000\n")
        File.write("        }\n")
        File.write("      }\n")
        File.write("      geometry Sphere{\n")
        File.write("        radius 1.0\n")
        File.write("      }\n")
        File.write("    }\n")
        File.write("  ]\n}\n")
        File.write("#Block " + str(i) + " End\n")    

class GMData:
  """
    A class to read, create, and write ellipsoids as GMM file
  """
  tol = 1e-8
  def __init__(self,weight=None,center=None,rotation=None,
               radii=None,makeCov=True,
               fromFile=False,filename=None,
               fromVRMLData=False,VRMLData=None,
               scale=1.0,molData=None,is_mol_data=False): 
    """
      The constructor which either read GMM data from a file, or make GMM 
      data from inputs

      This class can be instantiated in 4 ways.
      1) gm = GMData(weight=Array_Of_Weight_per_ellipsoid,
                     center=Array_Of_Centers,
                     rotation=Array_Of_Rotations,
                     radii=Array_Of_Radii,makeCov=True)
      2) gm = GMData(weight=Array_Of_Weight_per_ellipsoid,
                     center=Array_Of_Centers,
                     rotation=Array_Of_Covmat,
                     makeCov=False)
      2) gm = GMData(fromFile=True,filename=Input_GMM_File)
      3) gm = GMData(fromVRMLData=True,VRMLData=VRML_Class_Instance)

      center: Array of centers of ellipsoids
      rotation: Array of rotation matrices of ellipsoids
      radii: Array of radii of ellipsoids      
      makeCov: Whether to make covariance matrix
      fromVRMLData: Whether to make GMM data from VRML-data
      VRMLData: An instance of class VRML
      fromFile: Whether to read GMM data from a file
      filename: The GMM file from which GMM data will be read
      scale: scale applied to covariance-matrix or radii
      molData: A dictionary of parameters for Gaussian Mixture (the summed
               distribution) with keys, weight: Total weight of all 
               Gaussians, mean: Center of all Gaussians, cov: Covariance
               matrix for the entire mixture, icov: Inverse of covariance
               matrix, constant: Constant factor for GDF, read: Whether to
               to read molData from file or not.
      is_mol_data: molData['read']
    """
    path = "clib/Rotate_lib.so"
    CFile = pkg_resources.resource_filename(__name__,path)
    self._Rotate_lib=ctypes.cdll.LoadLibrary(CFile)
    #self._Rotate_lib=ctypes.cdll.LoadLibrary(clib_dir+'./clib/Rotate_lib.so')
    path = "./clib/get_abs_struct_factor.so"
    CFile = pkg_resources.resource_filename(__name__,path)
    self._get_abs_struct_factor=ctypes.cdll.LoadLibrary(CFile)
    #self._get_abs_struct_factor=ctypes.cdll.LoadLibrary(clib_dir+'./clib/get_abs_struct_factor.so')
    self._get_abs_struct_factor.get_abs_struct_factor.restype=ctypes.c_double
    if fromFile:
      # Read from file, molData read from input file within read method if available
      header, comment, weight, center, covmat, detcov, detcons, molData = self.read(filename,scale=scale,verbose=0)
      if not detcov is None: self.detcov = detcov
      if not detcons is None: self.detcons = detcons
      is_mol_data = molData['read']
      self.header = header
      self.comment = comment 
      # No need to make covariance from rotation and radii
      makeCov = False   
    else:
      if fromVRMLData and VRMLData is not None:
        # set from VRMLData
        n = len(VRMLData.covmat)        
        covmat = [self.getcovExpand(VRMLData.covmat[i]) for i in range(n)] 
        center = [VRMLData.translation[i] for i in range(n)]
        weight = list(VRMLData.weight); # weight = np.repeat(1.0 / n,n)
        self.detcov = None; self.detcons = None; 
        # Make a dummy molData, as VRMLData does not have any molData, and
        # set read to False
        molData = {'weight': None, 'mean': [] * 3, 'cov': np.zeros((3,3)), 
                   'icov': np.zeros((3,3)),'constant': None, 'read': False}
        is_mol_data = molData['read']
        # No need to make covariance from rotation and radii
        makeCov = False
      else:
        # Set from inputs, then covmat initialized from rotation
        covmat = rotation
      self.header = 'HEADER 3D Gaussian Mixture Model'
      self.comment = 'REMARK COMMENT Created from GMData Class'
    self.N = len(center)
    self.covmat = []; self.center = center; self.W = weight;
    for i in range(self.N):
      if makeCov:
        # Make covariance matrix from rotation and radii, scale 
        # is scaling applied to covariance matrix as a whole
        # radii should be unscaled radii
        r = rotation[i]; s = radii[i];
        cv = invert_decomposition(r,s)
        cv = cv * scale
      else:
        # Make covariance matrix from saved variable before, scale 
        # is scaling applied to covariance matrix as a whole
        cv = covmat[i] * scale
      self.covmat.append(cv)

    # Generate Gaussian for Molecule (Mixture) or Set
    if not is_mol_data:
      self.generate_gauss_mol()
    else:
      self.set_molData(molData)
    self.celmName = ["xx","xy","xz","yy","yz","zz"]

  def getcovExpand(self,cmat):
    """
      Expand covariance from minimal representation to full (3x3)
      cmat: Minimal upper-triangular covariance matrix
      returns: Full covariance matrix (3x3)
    """
    cv = np.zeros((3,3))
    cv[0,0] = cmat[0]; cv[0,1] = cmat[1]; cv[0,2] = cmat[2];
    cv[1,1] = cmat[3]; cv[1,2] = cmat[4]; cv[2,2] = cmat[5];
    cv[1,0] = cv[0,1]; cv[2,0] = cv[0,2]; cv[2,1] = cv[1,2];
    return cv

  def read(self,filename,scale=1.0,verbose=0):
    """
      Read GMM-file
      filename: Input filename
      scale: scale applied to covariance matrix
      verbose: Verbosity
      returns: tuple of header, comment, weight, center, covmat, 
                        detcov, detcons, molData
         where, header = head string written to GMM-file
                comment = comment string written to GMM-file
                weight = weight of each Gaussian
                center = center of each Gaussian
                covmat = covariance matrix of each Gaussian
                detcov = Information followed by 'det' word
                         Determinant of covariance matrix
                detcons = Information followed by cons word
                         Constant related to GDF
                molData = A dictionary of parameters for Gaussian Mixture
                          See document from constructor
    """
    if exists(filename):
      header = []; remark = []; hetatm = []
      with open(filename) as fp:
        for line in fp:
          line = line.strip()
          lineType = line[0:6]
          if lineType == "HEADER": header.append(line)
          if lineType == "REMARK": remark.append(line)
          if lineType == "HETATM": hetatm.append(line)
      header = "\n".join(header)
      comment = []; gaussLine = []; metadata = []; gaussMolLine = [];
      n = None; gaussMolBool = False;
      for r in remark:
        words = r.split()
        remType = words[1]
        if remType == "COMMENT":
          comment.append(r)
        elif remType == "NGAUSS":
          n = int(words[2])
        elif remType == "GAUSS":
          gaussLine.append(r)
        elif remType == "GAUSS_MOLECULE":
          gaussMolLine.append(r)
          gaussMolBool = True
        else:
          metadata.append(r)
      if scale != 1.0: gaussMolBool = False
      comment = "\n".join(comment)
      metadata = "\n".join(metadata)
      comment = metadata + comment
      if n is None:
        raise Exception("file %s does not contain NGAUSS string" % filename)
      else:
        if verbose > 0:
          print("%d Gaussians are read" % n)
      weight = [None] * n
      center = [None] * n
      covmat = [np.zeros((3,3)) for i in range(n)]
      detcov = [None] * n
      detcons = [None] * n
      dd = {'x':0,'y':1,'z':2}
      for gl in gaussLine:
        words = gl.split()
        index = int(words[2])
        if index >= 1 and index <= n:
          valueType = words[3]
          if valueType == "W": 
            weight[index-1] = float(words[4])
          elif valueType == "det":
            detcov[index-1] = float(words[4])
          elif valueType == "cons":
            detcons[index-1] = float(words[4])
          elif valueType == "M":
            center[index-1] = [float(words[i+4]) for i in range(3)]
          elif valueType == "CovM":
            wa = np.array(words)
            index1 = [4,6,8]
            index2 = [5,7,9]
            elname = list(wa[index1])
            celm = [float(elm) for elm in list(wa[index2])]
            for ie,ename in enumerate(elname):
              epos1s, epos2s = list(ename)
              epos1 = dd[epos1s]; epos2 = dd[epos2s];
              covmat[index-1][epos1,epos2] = celm[ie]
      n1 = 0; n2 = 0;
      for i in range(n):
        if detcov[i] is None: n1 += 1
        if detcons[i] is None: n2 += 1
        covmat[i] = covmat[i] + covmat[i].T - np.diag(np.diag(covmat[i]))
      if n1 == n: detcov = None
      if n2 == n: detcons = None     
      if gaussMolBool:
        print("Generate gauss mol data by reading")  
        molData = {'weight': None, 'mean': [] * 3, 'cov': np.zeros((3,3)),
                   'icov': np.zeros((3,3)),'constant': None, 'read': True}
        for gl in gaussMolLine:
          words = gl.split()
          valueType = words[2]
          if valueType == "Weight":
            molData['weight'] = float(words[3])
          if valueType == "Gcen":
            molData['mean'] = [float(words[i+3]) for i in range(3)]
          if valueType == "CovM":
            wa = np.array(words)
            index1 = [3,5,7]
            index2 = [4,6,8]
            elname = list(wa[index1])
            celm = [float(elm) for elm in list(wa[index2])]
            for ie,ename in enumerate(elname):
              epos1s, epos2s = list(ename)
              epos1 = dd[epos1s]; epos2 = dd[epos2s];
              molData["cov"][epos1,epos2] = celm[ie]
        molData["cov"] = molData["cov"] + molData["cov"].T - np.diag(np.diag(molData["cov"]))
        dcovar = np.linalg.det(molData["cov"])
        cons_mol = 1.0 / ((2.0 * np.pi)**1.5 * np.sqrt(dcovar))
        if np.abs(dcovar) <= GMData.tol: 
          molData["icov"] = np.linalg.pinv(molData["cov"],rcond = GMData.tol)
        else:
          molData["icov"] = np.linalg.inv(molData["cov"])
        molData["constant"] = cons_mol
      else:
        molData = {'weight': None, 'mean': [] * 3, 'cov': np.zeros((3,3)),
                   'icov': np.zeros((3,3)),'constant': None, 'read': False}
      return header, comment, weight, center, covmat, detcov, detcons, molData
    else:
      raise Exception("file: %s not found" % filename)

  def set_molData(self,molData):
    """
      Set Gaussian Molecule (Mixture) parameters
      and calculate some parameters not available by reading GMM-file
      molData: A dictionary of parameters for Gaussian Mixture
               See document from constructor
    """
    self.mol_weight = molData['weight']
    self.mol_mean = molData['mean']
    self.mol_cov = molData['cov']
    self.imol_cov = molData['icov']
    self.cons_mol = molData['constant']

    # Generate secondary data not available in molData dictionary
    self.generate_gauss_mol(partial=True)

  def generate_gauss_mol(self,partial=False):
    """
      Generate parameters for Gaussian Molecule (Mixture)
      partial: Whether to generate all data or secondary data only                     
    """
    if not partial:

      print("generate gauss mol data fully")

      # Total weight
      self.mol_weight = np.sum(self.W)

      # Center of all Gaussians
      self.mol_mean = np.full((self.N,3),0.0)
      for i in range(self.N):
        self.mol_mean[i,:] = np.array(self.center[i]) * self.W[i]
      self.mol_mean = np.sum(self.mol_mean,axis=0) / self.mol_weight

      # Covariance matrix for the Molecule (mixture)
      self.mol_cov = np.full((3,3),0.0)
      for i in range(self.N):
        m = self.center[i]; w = self.W[i]; cv = self.covmat[i]
        for k in range(3):
          mk = m[k] - self.mol_mean[k]
          for l in range(3):
            ml = m[l] - self.mol_mean[l]
            self.mol_cov[k,l] = self.mol_cov[k,l] + w*mk*ml
            self.mol_cov[k,l] = self.mol_cov[k,l] + w*cv[k,l]        
      self.mol_cov /= self.mol_weight      

    # Generate secondary data
    if not hasattr(self,'pc'):

      # Generate principal components
      val,vec = np.linalg.eig(self.mol_cov)
      val_index = np.argsort(-val)
      val = val[val_index]; vec = vec[:,val_index]
      if sum(np.cross(vec[:,0],vec[:,1]) * vec[:,2]) < 0.0: vec[:,2] = -vec[:,2]
      self.pc = {}
      for i in range(3):
        self.pc[i] = {'value': val[i], 'sd' : np.sqrt(val[i]), 'axis': vec[:,i]}

    if not hasattr(self,'cons_mol') or self.cons_mol is None:

      # Calculate constant for GDF of molecule (Mixture)
      dcovar = np.linalg.det(self.mol_cov)
      self.cons_mol = 1.0 / ((2.0 * np.pi)**1.5 * np.sqrt(dcovar))

    if not hasattr(self,'imol_cov') or self.imol_cov is None:

      # Calculte inverse of covariance matrix of molecule (Mixture)
      if np.abs(dcovar) <= GMData.tol: 
        self.imol_cov = np.linalg.pinv(self.mol_cov,rcond = GMData.tol)
      else:
        self.imol_cov = np.linalg.inv(self.mol_cov)

  def write(self,outFileName,minimal=True):
    """
      Write GMM-data to file
      outFileName: Output GMM file
      minimal: Whether to write determinants and constants for 
               each Gaussians
    """
    with open(outFileName,"w") as File:
      File.write(self.header + "\n" + self.comment + "\n")
      File.write("REMARK NGAUSS %d\n" % self.N)
      File.write("REMARK GAUSS_MOLECULE Weight %.6f\n" % self.mol_weight)
      buf = StringIO.StringIO()
      buf.write("REMARK GAUSS_MOLECULE Gcen     ")
      for i in range(3): buf.write("%8.3f " % self.mol_mean[i])        
      File.write(buf.getvalue().strip() + "\n")
      buf = StringIO.StringIO()
      buf.write("REMARK GAUSS_MOLECULE CovM ")
      cov = self.get_covariances_uniq(self.mol_cov)
      for i in range(0,3): buf.write(" %s %15.10f" % (self.celmName[i],cov[i]))
      buf.write("\nREMARK GAUSS_MOLECULE CovM ")
      for i in range(3,6): buf.write(" %s %15.10f" % (self.celmName[i],cov[i]))
      File.write(buf.getvalue() + "\n")
      for i in range(3):
        File.write("REMARK GAUSS_MOLECULE PCvar%d   %10.3f PCsd%d   %8.3f\n" % 
                   (i,self.pc[i]['value'],i,self.pc[i]['sd']))
      for i in range(3):
        File.write("REMARK GAUSS_MOLECULE PCaxis%d  " % i)
        buf = StringIO.StringIO()
        for j in range(3): buf.write("%8.3f " % self.pc[i]['axis'][j])
        File.write(buf.getvalue().strip() + "\n")
      File.write("REMARK GAUSS_MOLECULE homo_id 1000000\n")
      File.write("REMARK GAUSS_MOLECULE Chain\n")
      for i in range(self.N):
        k = i + 1
        File.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n" % 
                   (k,"GAU","GAU"," ",k,self.center[i][0],self.center[i][1],self.center[i][2],
                   self.W[i],self.W[i]))
        File.write("REMARK GAUSS%4d W %.10f\n" % (k,self.W[i]))
        if not minimal:
          flag = False
          try: 
            self.detcov
            flag = True
          except NameError:
            detcov = np.linalg.det(self.covmat[i])
          if flag:
            try:
              detcov = self.detcov[i]
            except IndexError:
              detcov = np.linalg.det(self.covmat[i])
          flag = False
          try:
            self.detcons
            flag = True
          except NameError:
            detcons = 1.0 / ((2.0 * np.pi)**1.5 * np.sqrt(detcov))
          if flag:
            try:
              detcons = self.detcons[i]
            except IndexError:
              detcons = 1.0 / ((2.0 * np.pi)**1.5 * np.sqrt(detcov))
          File.write("REMARK GAUSS%4d det  %.10f\n" % (k,detcov))
          File.write("REMARK GAUSS%4d cons %.10f\n" % (k,detcons))
        File.write("REMARK GAUSS%4d M %.6f %.6f %.6f\n" % 
                    (k,self.center[i][0],self.center[i][1],self.center[i][2]))
        buf = StringIO.StringIO(); 
        buf.write("REMARK GAUSS%4d CovM " % k)
        cov = self.get_covariances_uniq(self.covmat[i])
        for j in range(0,3): buf.write(" %s %15.10f" % (self.celmName[j],cov[j]))
        buf.write("\nREMARK GAUSS%4d CovM " % k)
        for j in range(3,6): buf.write(" %s %15.10f" % (self.celmName[j],cov[j]))
        File.write(buf.getvalue() + "\n")
      File.write("TER\n")  
    
  def get_covariances_uniq(self,cmat):
    """
      Return minimal covariance matrix, upper-triangular only
      cmat: Input covariance matrix
      returns: Upper-triangular part only of covariance matrix
    """
    return [cmat[0,0],cmat[0,1],cmat[0,2],cmat[1,1],cmat[1,2],cmat[2,2]]

  def _create_input_gauss(self,index):
    """
      Helper function for other class methods, aim to make a simple
      dictionary with key 'mu' or center of a Gaussian and 's' or 
      covariance matrix of the Gaussian, given an index.
      Note: Index must be within 0 <= index < N, N = number of Gaussians
      index: Input index of Gaussian
      returns: A dictionary of with single Gaussian parameters
    """
    g = {'mu': np.array(self.center[index]), 's': self.covmat[index]}
    return g

  def create_ellipsoid_data(self,scale=1.0,npoint=1000,prune=True):  
    """
      Creates Ellipsoid data (Array of Ellipsoids) from Gaussian data
      scale: Input scale applied to Ellipsoid covariance
      npoint: Number of points resampled in total within all Gaussians
              Note: This number of points will be equally distributed among
                    all the Gaussians. So if there are many Gaussians, number
                    of points per Gaussian can be small
      prune: Whether resampled points to be pruned to exist only within 
             specified standard-deviation
      returns: Array of ellipsoids      
    """
    # sd_factor: A factor to be applied to the standard-deviation per 
    #             ellipsoid. The default value 5.0 is related to 
    #             Ellipsoid.chi2CoverRatio and maxsd in class Voxel
    sd_factor=5.0
    E = [None] * self.N
    if self.N > 1:
      PointPerEllipsoid = int(npoint / self.N)
      PointPerEllipsoid1 = npoint - (PointPerEllipsoid * (self.N - 1))
      N = np.concatenate((np.repeat(PointPerEllipsoid,self.N-1),np.array([PointPerEllipsoid1])))
    else:
      N = np.array([self.N])
    if prune:
      # First create Ellipsoids with small number of points (doesn't matter)      
      for i in range(self.N):      
        E[i] = Ellipsoid()
        E[i].set_fromGDF(mean=self.center[i],covmat=self.covmat[i],W=self.W[i],
                         scale=scale,resample_point=10,cutoff_pdf=None)
      # Using these primary data get the value of cutoff
      vox = Voxel()
      vox.buildVoxel(E,fromFile=False)
      cutoff = vox.getCutoff(sd_factor=sd_factor) 
      print("Cutoff = %e" % cutoff)
    else:
      cutoff = None
    # Now create ellipsoids with resampled points according to cutoff
    for i in range(self.N):        
      E[i] = Ellipsoid()
      E[i].set_fromGDF(mean=self.center[i],covmat=self.covmat[i],W=self.W[i],
                       scale=scale,resample_point=N[i],cutoff_pdf=cutoff)
      E[i].calculate_aux()
    return E
  def rotate(self,inputData):
    phi = inputData[0]
    theta = inputData[1]
    psi = inputData[2]
    Rz = self.construct_rotation_matrix(phi,'z')
    Rx = self.construct_rotation_matrix(theta,'x')
    Ry = self.construct_rotation_matrix(psi,'y')
    R = (Rz.dot(Ry)).dot(Rx)
    vert = np.zeros((self.N,3))
    for i in range(self.N):
      vert[i,:] = self.center[i]
      self.covmat[i] = (R.dot(self.covmat[i])).dot(R.T)
    cen = np.mean(vert,axis=0)
    vert = (R.dot((vert - cen).T)).T + cen
    for i in range(self.N):
      self.center[i] = vert[i,:]
    
  def construct_rotation_matrix(self,angle,axis):
    ca = np.cos(angle)
    sa = np.sin(angle)
    R = np.diag([1.0,1.0,1.0])
    if axis == 'z':      
      R[0,0] = ca; R[0,1] = sa; 
      R[1,0] =-sa; R[1,1] = ca;
    elif axis == 'x':
      R[1,1] = ca; R[1,2] = sa; 
      R[2,1] =-sa; R[2,2] = ca;
    elif axis == 'y':
      R[0,0] = ca; R[0,2] = sa; 
      R[2,0] =-sa; R[2,2] = ca;
    return R

#Nagai san's code
  

  def init_ctypes_of_gmm(self):
        """set ctypes object of gmm parameters"""
        self.ncenter=ctypes.c_int(self.N) # self.nmu is the ctypes of self.ng. Bad naming, which should be fixed in some future... #2

        """ initialize ctypes """
        self.weight_array_class =ctypes.c_double *self.N # 2
        self.weight_array=self.weight_array_class() # 1
        
        self.center_array_class=ctypes.c_double*(self.N*3) # 2
        self.center_array=self.center_array_class() #2

        self.covmat_array_class=ctypes.c_double*(self.N*3*3) # 2
        self.covmat_array=self.covmat_array_class() #2

        """ make ojbects that work as pointer """
        self.weight_ref=ctypes.byref(self.weight_array) # 2
        self.center_ref=ctypes.byref(self.center_array) #2
        self.covmat_ref=ctypes.byref(self.covmat_array) #2

        """copy information in numpy to the ctypes's object """
        self.put_data_ctypes_of_gmm()

        self._is_init_ctypes_of_gmm=True

  def put_data_ctypes_of_gmm(self):
        """put numpy's data into ctypes object"""
        self.ncenter.value=self.N # 2
        self.weight_array[:]=self.W[:] # 2(weight-W)
        c = np.array(self.center)
        self.center_array[:]=c.ravel()[:] #2
        d = np.array(self.covmat)
        self.covmat_array[:]=d.ravel()[:] #2
        
  def call_get_abs_struct_factor(self,kx,ky,kz):
        """ 
        YOU HAVE TO call "init_ctypes_of_gmm" before CALLING THIS!!
        the unit of k must be 1/A. 
        you can usally use this if you have already calculated diffraction image by calling xxxx
        As information won't be passed into ctypes objects in this function, you should care whether the ctypes objects contain intended information."""
        self._diff_one_point=self._get_abs_struct_factor.get_abs_struct_factor(ctypes.c_double(kx),ctypes.c_double(ky),ctypes.c_double(kz),self.weight_ref,self.center_ref,self.covmat_ref,self.ncenter) #4
        return self._diff_one_point

        
  def call_get_abs_struct_factor_all_mesh_on_Ewald_sphere(self, number_grid, det_width, det_dist, wave_length):
        """you can usally use this if you have calculated diffraction once with function like calc_diffraction_with_c. 
        As the information won't be passed in this function, you should care whether the ctypes objects contain intended information."""
        self._number_grid=number_grid
        self._diffraction_image_array_class=ctypes.c_double*(number_grid*number_grid)
        self._diffraction_image_array=self._diffraction_image_array_class()
        self.diff_ref=ctypes.byref(self._diffraction_image_array)

        self._get_abs_struct_factor.get_abs_struct_factor_all_grid_point\
         (self.diff_ref, \
         ctypes.c_int(number_grid), ctypes.c_double(det_width), ctypes.c_double(det_dist), ctypes.c_double(wave_length),\
         self.weight_ref,self.center_ref,self.covmat_ref,self.ncenter) #4
        self.diff=np.array(self._diffraction_image_array[:]).reshape(number_grid,number_grid)
        
        return self.diff


# Added 07/25/2018
class Voxel:
  """
    This class is for making Voxel data from Ellipsoids or GMData.
  """
  inv2pi = 1.0 / (2.0 * np.pi)**(1.5)
  def __init__(self,gw=4.0,maxsd=5.0,origin=None,numGrid=None): 
    """
      Constructor sets parameters to calculate probability density
      values within voxel
      gw: Grid-width
      maxsd: A factor to be multiplied to standard-deviations of
             covariances.
      origin: A tuple or list to set the origin
      numGrid: A tuple or list tp set the number of grids along xyz
    """  
    self.setOriginBool = False
    self.setNumGridBool = False
    self.setup = False
    self.setGridWidth(gw)
    self.maxsd = maxsd
    if origin is not None:
      origin = self._checkVector(origin,Name = "Origin")
      self.setOrigin(origin)    
    if numGrid is not None:
      nGrid = self._checkVector(numGrid,Name = "NumGrid")
      self.setNumGrid(nGrid)   

  def _checkVector(self,v,Name="Vector"):
    """
      Helper for constructor to check input origin and numgrid
      tuple or list and make numpy array
      v: Input vector
      Name: A string to print information
      returns: A numpy array
    """
    if isinstance(v,tuple): 
      v = list(v)
    elif isinstance(v,list):
      pass
    else:
      raise("%s (if given) must be a list or a tuple" % Name)
    if len(v) != 3:
      raise("%s (if given) must be of length 3" % Name)   
    return np.array(v)

  def buildVoxel(self,ellipsoids=None,
                 fromFile=False,filename=None):
    """
      Build Voxel data from ellipsoids or from GMM-file
      ellipsoids: Array of Ellipsoids
      fromFile: Whether to read from a GMM-file
      filename: Input GMM Filename 
    """
    if not self.setup:
      self.setupVoxel(ellipsoids=ellipsoids,fromFile=fromFile,filename=filename)
    self.createGrid()
    self.CalProbDensity()        

  def setupVoxel(self,ellipsoids=None,fromFile=False,filename=None,scale=1.0):
    """
      Sets up GMM data from ellipsoids or from GMM-file
      ellipsoids: Array of Ellipsoids
      fromFile: Whether to read from a GMM-file
      filename: Input GMM Filename 
    """
    if fromFile:

      # Read from GMM file using GMData class
      self.gm = GMData(fromFile=True,filename=filename)
      self.nGaussian = len(self.gm.center)
      self.weight = [None] * self.nGaussian
      self.pnorm = [None] * self.nGaussian
      for i in range(self.nGaussian):
        center = self.gm.center[i]
        covmat = self.gm.covmat[i] * scale
        self.weight[i] = self.gm.W[i]
        self.pnorm[i] = multivariate_normal(mean=center,cov=covmat,allow_singular=True)

    else:

      # Get data from input ellipsoids
      self.nGaussian = len(ellipsoids)
      center = [None] * self.nGaussian; covar = [None] * self.nGaussian;
      self.pnorm = [None] * self.nGaussian             
      self.weight = np.repeat(1.0 / self.nGaussian,self.nGaussian)
      for i,e in enumerate(ellipsoids):
        center[i] = e.mean
        covar[i] = e.covar * scale       
        self.pnorm[i] = multivariate_normal(mean=e.mean,cov=scale * e.covar,allow_singular=True)
      # Call GMData class using data obtained
      self.gm = GMData(weight=self.weight,center=center,rotation=covar,makeCov=False)

    # Get Gaussian for molecule (Mixture)
    self.mean = self.gm.mol_mean
    self.covar = self.gm.mol_cov
    self.sd = np.sqrt(np.diag(self.covar))

    # Set Origin, and Number of grids if not set already
    if not self.setOriginBool: self.setOrigin(self.getOrigin())
    if not self.setNumGridBool: self.setNumGrid(self.getNumGrid())
    self.setup = True

  def getOrigin(self):
    """
      Return Origin by calculating from mean, maxsd and sd of Gaussian
      molecule (Mixture) or return already set value  
      returns: origin of voxel
    """
    if self.setOriginBool: return self.Origin
    if hasattr(self,"mean") and hasattr(self,"maxsd") and hasattr(self,"sd"):
      return self.mean - (self.maxsd * self.sd)    
    else:
      print("Can't calculate Origin")
      return None

  def setOrigin(self,origin):
    """
      Set origin to the input
      origin: Input origin of voxel
    """
    self.Origin = origin
    self.setOriginBool = True

  def getNumGrid(self):
    """
      Return number of grids by calculating from maxsd, gw and sd (of Gaussian 
      molecule (Mixture)) or return already set value 
      returns: number of grids along xyz
    """
    if self.setNumGridBool: return self.numGrid
    if hasattr(self,"maxsd") and hasattr(self,"sd") and hasattr(self,"gw"):
      numgrid = np.ceil((2.0 * self.maxsd * self.sd) / self.gw)
      return np.asarray(numgrid,dtype=int)
    else:
      print("Can't calculate Number of grids")
      return None

  def setNumGrid(self,nGrid):
    """
      Set number of grids to the input
      nGrid: Input number of grids along xyz
    """
    self.numGrid = nGrid
    self.setNumGridBool = True

  def getNumGridbyGridWidth(self,gw):
    """
      Calculate number of grids along xyz using input grid-width and previous
      value of grid-width
      Note: If you set this number of grid in class then make sure to update
            grid-width also
      gw: Input grid-width
      returns: Updated number of grids along xyz
    """
    ng = ((self.numGrid * self.gw) / gw).astype(int)
    return ng

  def setGridWidth(self,gw):
    """
      Sets grid-width to the input
      gw: Input grid-width
    """
    self.gw = gw

  def createGrid(self):
    """
      Create grids along xyz using number of grids along xyz, grid-width
      and origin of voxel
    """
    nx =  self.numGrid[0]; ny = self.numGrid[1] ; nz = self.numGrid[2];
    self.rx = np.asarray(range(nx)) * self.gw + self.Origin[0]
    self.ry = np.asarray(range(ny)) * self.gw + self.Origin[1]
    self.rz = np.asarray(range(nz)) * self.gw + self.Origin[2]    

  def CalProbDensity(self):    
    """
      Calculate probability density at each grid-point
    """
    nx =  self.numGrid[0]; ny = self.numGrid[1] ; nz = self.numGrid[2];
    totalNumPoint = nx * ny * nz
    self.ProbDensity = np.zeros((nx,ny,nz))    
    rx = self.rx; ry = self.ry; rz = self.rz;
    point = np.array(np.meshgrid(rx, ry, rz)).reshape(3, totalNumPoint).T
    index = np.array(np.meshgrid(range(nx), range(ny), range(nz))).reshape(3, totalNumPoint).T
    pdf = np.zeros(totalNumPoint)
    for i,p in enumerate(self.pnorm):      
      pdf += self.weight[i] * p.pdf(point)
    for ip in range(index.shape[0]):
      ijk = index[ip,:]
      self.ProbDensity[ijk[0],ijk[1],ijk[2]] = pdf[ip]

  def WriteMRC(self,outFileName,format='map'):
    """
      Write MAP(ccp4) or MRC file
      outFileName: Output filename
      format: Either 'map' or 'mrc'
    """
    mrc = MRC(voxelData=self.ProbDensity,voxelOrigin=self.Origin,gridWidth=self.gw)
    mrc.write(outFileName,format=format)

  def getCutoff(self,sd_factor=None):
    """
      Get cutoff value to show data
      sd_factor: Scaling applied to standard-deviation of probability density
    """
    # Just use WriteMRC's functionality
    mrc = MRC(voxelData=self.ProbDensity,voxelOrigin=self.Origin,gridWidth=self.gw)
    vmin = mrc.AMIN; vmax = mrc.AMAX;
    self.vavg = mrc.AMEAN; self.vstd = mrc.ARMS;
    if sd_factor is None: sd_factor = self.maxsd
    return self.vavg + self.vstd * sd_factor

  def setCutoff(self,cutoff):
    """
      Set cutoff value to the input
      cutoff: Input cutoff value
    """
    self.cutoff = cutoff

  def CalMarchingCube(self):    
    """
      Performs marching cubes algorithm on probability density for 
      isosurface calculation
      Note: Uses skimage.measure
    """
    self.v,self.s,self.n,self.nv = measure.marching_cubes(self.ProbDensity,level=self.cutoff)

  def getIsoSurface(self):
    """
      Get Iso-surface on probability density using vertices and simplices 
      obtained from marching cubes algorithm
      returns: A set of vertices that defines iso-surface
    """
    vertices = self.v
    simplices = self.s
    return vertices[simplices]

  def getPrunedBox(self):
    """
      Get Pruned box limits based on iso-surface 
      returns: Limits along xyz based on iso-surface
    """
    vertices = self.v
    xyzMin = np.min(vertices,axis=0); xyzMax = np.max(vertices,axis=0)
    xlim = [xyzMin[0], xyzMax[0]]; ylim = [xyzMin[1], xyzMax[1]]; zlim = [xyzMin[2], xyzMax[2]];
    return xlim,ylim,zlim

class MRC:
  """
    This sets MRC or CCP4-MAP data and writes to a file
  """
  def __init__(self,voxelData,voxelOrigin,gridWidth):
    """
      Constructor based on voxel data, origin of voxel and grid-width
    """
    self.origin = voxelOrigin
    self.data = voxelData
    self.gw = gridWidth
    self.angle = 90.0

    # Calculate statistics for writing
    self.getStatistics()

    # Set variables for writing
    self.setParam()

  def getStatistics(self):
    """
      Calculate different statistics based on voxel data
    """
    self.AMIN  = np.asscalar(np.min(self.data))
    self.AMAX  = np.asscalar(np.max(self.data))
    self.AMEAN = np.asscalar(np.mean(self.data))
    delta = self.data - self.AMEAN
    self.ARMS  = sqrt(np.mean((delta*delta).ravel()))   
 
  def setParam(self):
    """
      Sets different variables which will be used for writing
      map(ccp4) or mrc file
    """
    mrcOrigin = np.asarray(self.origin / self.gw,dtype=int)
    boxLen = np.asarray(self.data.shape) * self.gw
    self.NCSTART = np.asscalar(mrcOrigin[0])
    self.NRSTART = np.asscalar(mrcOrigin[1])
    self.NSSTART = np.asscalar(mrcOrigin[2])
    print("NCSTART, NRSTART, NSSTART = %d %d %d" % (self.NCSTART,self.NRSTART,self.NSSTART))
    self.NC = self.data.shape[0]; self.NR = self.data.shape[1]; self.NS = self.data.shape[2];
    self.NX = self.NC; self.NY = self.NR; self.NZ = self.NS;
    self.XLEN = np.asscalar(boxLen[0])
    self.YLEN = np.asscalar(boxLen[1])
    self.ZLEN = np.asscalar(boxLen[2]); 
    self.CELLALPHA = self.angle; self.CELLBETA = self.angle; self.CELLGAMMA = self.angle;
    self.MAPC = 1; self.MAPR = 2; self.MAPS = 3;
    self.ISPG = 1; self.NSYMBT = 0;
    self.MODE = 2
    self.ORIGINX = np.asscalar(self.origin[0])
    self.ORIGINY = np.asscalar(self.origin[1])
    self.ORIGINZ = np.asscalar(self.origin[2])    
    print("ORIGIN = %f %f %f" % (self.ORIGINX,self.ORIGINY,self.ORIGINZ))
    self.MAP = "MAP "; self.MACHST = "DA"
    self.NLABEL = 10
    self.LABELNULL = " "
    self.LABEL0 = "This is made by AFM Module"
    self.LABEL1 = "This follows GMCONVERT implementation"

  def write(self,outFileName,format='map'):
    """
      Writes data to file follwing given format
      outFileName: Output filename
      format: Either 'map' or 'mrc'
    """
    with open(outFileName, "wb") as File:
      File.write(struct.pack('@i', self.NC))
      File.write(struct.pack('@i', self.NR))
      File.write(struct.pack('@i', self.NS))
      File.write(struct.pack('@i', self.MODE))
      File.write(struct.pack('@i', self.NCSTART))
      File.write(struct.pack('@i', self.NRSTART))
      File.write(struct.pack('@i', self.NSSTART))
      File.write(struct.pack('@i', self.NX))
      File.write(struct.pack('@i', self.NY))
      File.write(struct.pack('@i', self.NZ))
      File.write(struct.pack('@f', self.XLEN))
      File.write(struct.pack('@f', self.YLEN))
      File.write(struct.pack('@f', self.ZLEN))
      File.write(struct.pack('@f', self.CELLALPHA))
      File.write(struct.pack('@f', self.CELLBETA))
      File.write(struct.pack('@f', self.CELLGAMMA))
      File.write(struct.pack('@i', self.MAPC))
      File.write(struct.pack('@i', self.MAPR))
      File.write(struct.pack('@i', self.MAPS))
      File.write(struct.pack('@f', self.AMIN))
      File.write(struct.pack('@f', self.AMAX))
      File.write(struct.pack('@f', self.AMEAN))
      File.write(struct.pack('@i', self.ISPG))
      File.write(struct.pack('@i', self.NSYMBT))
      zero = 0; zero_float = float(0.0);
      if format == 'map':
        File.write(struct.pack('@i', zero))
        for x in range(9): File.write(struct.pack('@f', zero_float))
        for x in range(3): File.write(struct.pack('@f', zero_float))
        for x in range(38,53): File.write(struct.pack('@i', zero))
      elif format == 'mrc': 
        for x in range(25,50): File.write(struct.pack('@i', zero))
        File.write(struct.pack('@f', self.ORIGINX))
        File.write(struct.pack('@f', self.ORIGINY))
        File.write(struct.pack('@f', self.ORIGINZ))
      File.write(struct.pack('@4s', self.MAP))
      File.write(struct.pack('@4s', self.MACHST))
      File.write(struct.pack('@f', self.ARMS))
      File.write(struct.pack('@i', self.NLABEL))
      File.write(struct.pack('@80s', self.LABEL0))
      File.write(struct.pack('@80s', self.LABEL1))
      for x in range(2,10): File.write(struct.pack('@80s', self.LABELNULL))
      for k in range(self.NZ):
        for j in range(self.NY):
          for i in range(self.NX):
            File.write(struct.pack('@f', np.asscalar(self.data[i,j,k])))

# End Added 07/25/2018

class Overlap:
  """
    Calculation of overlap between two Gaussians
  """
  tol = 1e-8

  def __init__(self,g1=None,g2=None):
    """
      Constructor to set input Gaussians
      g1: First Gaussian
      g2: Second Gaussian
    """
    self.g1 = g1
    self.g2 = g2

  def set_input_gaussian(self,g,datatype="1"):
    """
      Set first or second Gaussian data
      g: Input Gaussian
      datatype: Either '1' or '2', if '1' sets first Gaussian,
                if '2' sets second Gaussian
    """
    if datatype == "1":
      self.g1 = g
    elif datatype == "2":
      self.g2 = g

  def set_input_gaussians(self,g1,g2):
    """
      Set first and second Gaussians
      g1: First Gaussian
      g2: Second Gaussian
    """
    self.g1 = g1
    self.g2 = g2

  def _get_mu(self,g):
    """
      A helper method to return center of input Gaussian
      g: Input Gaussian
      return: center of the input Gaussian
    """
    return g['mu']

  def _get_sigma(self,g):
    """
      A helper method to return covariance matrix of input Gaussian
      g: Input Gaussian
      return: covariance matrix of the input Gaussian
    """
    return g['s']

  def get_overlap_gaussian(self):
    """
      Calculate overlap between Gaussians
      returns: The value of the overlap
    """
    S1 = self._get_sigma(self.g1)
    S2 = self._get_sigma(self.g2)
    add_S1_S2 = S1 + S2
    det_add_S1_S2 = np.linalg.det(add_S1_S2)
    if np.abs(det_add_S1_S2) <= Overlap.tol:
      inv_add_S1_S2 = np.linalg.pinv(add_S1_S2,rcond = Overlap.tol)
    else:
      inv_add_S1_S2 = np.linalg.inv(add_S1_S2)
    M1 = self._get_mu(self.g1)
    M2 = self._get_mu(self.g2)
    substract_M1_M2 = M1 - M2    
    q = self._Quadratic_Form_3D(inv_add_S1_S2,substract_M1_M2)
    #print("qform = %e\n" % q)
    ov = exp(-0.5*q) * (1.0 / ((2.0 * np.pi)**1.5)) / np.sqrt(det_add_S1_S2)    
    return ov

  def _Quadratic_Form_3D(self,mat,mu):
    """
      A helper method to calculate Quadratic term from given covariance matrix
      and center
      mat: Covariance matrix
      mu: Center
      returns: The quadratic term (mu^T mat mu)
    """
    qform  = mat[0,0]*mu[0]*mu[0] + mat[1,1]*mu[1]*mu[1] + mat[2,2]*mu[2]*mu[2];
    qform += 2.0*(mat[0,1]*mu[0]*mu[1] + mat[0,2]*mu[0]*mu[2] + mat[1,2]*mu[1]*mu[2]);
    return qform

class VolCompare:
  """
    Calculate correlation coefficient between two Gaussian mixtures
  """
  def __init__(self,gmdata1,gmdata2,base=None,algo=None):
    """
      Constructor to set first and second Gaussian mixture model
      gmdata1: First Gaussian mixture model of class GMData
      gmdata2: Second Gaussian mixture model of class GMData
      base: Either None or a reference Gaussian mixture model of class GMData
            This reference data can be used to fit (rigid-body transformation)
            second GMM to the reference GMM. Note that, in this fitting only
            centers of Gaussians are used. It is assumed that reference GMM
            and second GMM contains same number of GDFs with the similar
            correspondance, i.e. i-th GDF from reference GMM corresponds
            to i-th GDF of second GMM
      algo: The algorithm to fit. The choices are, 'trans' to fit by
            translation only, 'rot' to fit rotation only, 'both' to fit by
            both translation and rotation, 'transx' to fit by translation
            along x, 'transy' to fit by translation along y, 'transz' to fit
            by translation along z
    """
    self.gm1 = gmdata1; 
    self.cc11 = None; self.cc12 = None; self.cc22 = None;
    if algo is not None:
      if base is not None:
        self.gm2 = self.remove_drift_gm(gmdata2,base,algo=algo)
    else: self.gm2 = gmdata2

  def set_data(self,data,datatype="1"):
    """
      Set first or second Gaussian data
      data: Input Gaussian
      datatype: Either '1' or '2', if '1' sets first Gaussian,
                if '2' sets second Gaussian
    """
    if datatype == "1":
      self.gm1 = data
    if datatype == "2":
      self.gm2 = data

  def remove_drift_gm(self,source,target,algo='both'):
    """
      This method removes drift (i.e. fit) source to target by moving
      source.
      source: The source gmm
      target: The target gmm
      algo: The algorithm to fit. See documentation of constructor
      returns: The modified source gmm
    """
    print("removing drift")

    # Get centers
    source_pc = np.array(self.get_centers_gm(source))
    target_pc = np.array(self.get_centers_gm(target))

    # Treats centers as list of corresponding points fit by using
    # class pcFit
    pc_fit = pcFit(target=target_pc,source=source_pc,algo=algo)
    param = pc_fit.get_fit_param()
    source_pc2 = pc_fit.apply_fit_param(source_pc,param)

    # Set source centers
    for ic in range(len(source.center)):
      source.center[ic] = source_pc2[ic,:]
    source.centers = pc_fit.fit()    

    return source

  def get_centers_gm(self,gm):
    """
      Get centers of the GMM (the centers of the ellipsoids)
      gm: Input GMM of class GMData
      returns: Centers of GMM
    """
    return gm.center

  def getCC(self):
    """
      Calculate correlation coefficient
    """
    if self.cc11 is None:
      self.cc11 = self._overlap_gmm(self.gm1,self.gm1)
    if self.cc12 is None:
      self.cc12 = self._overlap_gmm(self.gm1,self.gm2)
    if self.cc22 is None:    
      self.cc22 = self._overlap_gmm(self.gm2,self.gm2)
    return self.cc12 / np.sqrt(self.cc11 * self.cc22)

  def setCC(self,cc,cctype="self1"):
    """
      Set sum of overlaps
      cc: The sum of overlaps
      cctype: Either 'self1' for sum of overlaps among first
              GMM, or 'self2' for sum of overlaps among second
              GMM, or 'cross' for sum of overlaps between first
              and second GMMs
    """
    if cctype == "self1": self.cc11 = cc
    if cctype == "self2": self.cc22 = cc
    if cctype == "cross": self.cc12 = cc

  def _overlap_gmm(self,gm1,gm2):
    """
      Calculate overlap values and sum them using class Overlap
      gm1: First GMM
      gm2: Seconf GMM
      returns: sum of overlap values
    """
    n = gm1.N
    m = gm2.N
    cc = np.zeros((n,m))
    for i in range(n):
      g1 = gm1._create_input_gauss(i)      
      for j in range(m):        
        g2 = gm2._create_input_gauss(j)        
        ovObj = Overlap(g1,g2);
        ov_ij = ovObj.get_overlap_gaussian()
        cc[i,j] = gm1.W[i] * gm2.W[j] * ov_ij
    return np.sum(cc) 

class AFMInputImage:
  """
    For working with AFM-like image
  """
  def __init__(self,data,fromarray=False,mode='up_is_dark',normalize=True,bg=255,verbose=0):
    """
      Constructor to set AFM-like grayscale image
      data: Input data from which class variables are set
            if data = matrix, then set fromarray = True
            if data = filename, then set fromarray = False
      fromarray: Whether to treat input data as a matrix containing image
                 or a string of filename
      mode: Either 'up_is_dark' (higher z-values are darker than lower 
            z-value) or 'up_is_light' (higher z-values are lighter than 
            lower z-value)
      normalize: Whether to normalize the matrix data to be within 0 and 
                 background value
      bg: The background value of color
      verbose: Verbosity
    """    
    self.bg = bg

    # Set mode
    if mode == 'up_is_dark':    self.mode = 1
    elif mode == 'up_is_light': self.mode = 2
    else: self.mode = 0
    if self.mode == 0:
      raise Exception("Invalid mode")

    # Read
    if fromarray:
      # Read from array
      if normalize: data = self.normalize_data(data)
              
      im = Image.fromarray(data)
      if verbose > 0: print("Fetching image from data")
    else: 
      # Read from file
      im = Image.open(data,'r')
      if verbose > 0: print("Fetching image from file")
    # Convert to grayscale data and convert mode
    im = im.convert('L')
    mat = np.asarray(im.getdata(),dtype=np.float64)
    mat = mat.reshape((im.size[1],im.size[0]))        
    if self.mode == 2:
      im, mat = self.convert_mode(mat)    

    # Save to class variables
    self.im = im; self.imat = mat
    self.imClean = self.im; self.imatClean = self.imat;   
  
  def normalize_data(self,data):    
    """
      Minmax centering of input data
      data: Input data matrix
      returns: normalized matrix
    """
    data = (data - np.min(data)) / (np.max(data) - np.min(data))
    return data*self.bg

  def convert_mode(self,imat,frm=2,to=1):
    """
      convert between modes
      imat: input image matrix
      frm: The mode from convert
      to: The mode to convert
      returns: tuple of converted image and corresponding matrix
    """
    if frm == 2 and to == 1:
      fill = np.full(imat.shape,self.bg)
      imatConv = np.asarray(np.abs(imat-fill),dtype=np.uint8)
      imConv = Image.fromarray(imatConv)
    else:
      raise Exception("Incomplete support")
    return (imConv,imatConv)

  def rm_background(self,hilo=(None,3),reverse=False,nbin=50):
    """ 
      Removes background from image matrix
      hilo: A tuple of high and low pixel value to auto-detect
            background
      reverse: Either True to make an image with mode up_is_light
               or False to make an image with mode up_is_dark
      nbin: Number of bins to create color histogram
      Note: Produce an image (a matrix with type unit8), where
      entries range from 0-255. Background is substracted.
      If reverse is True, then image value is proportional
      to the height. That is, image mode is up_is_light
      This mode is good to compare to a normalized height
      warp projection from xyz data
      If reverse is False, then image value is inversely
      proportional to the height, or it is up_is_dark.
      By default, reverse is false, so we get background
      substracted afm-image.
    """ 
    self.reverse = reverse
    sobel = filters.sobel(self.im)
    blurred = filters.gaussian(sobel,sigma=2.0)
    hi = hilo[0]; lo = hilo[1];
    if hi is None: hi = self.getHicut(nbin=nbin)
    light_spots = np.array((self.imat > hi).nonzero()).T            
    dark_spots = np.array((self.imat < lo).nonzero()).T
    bool_mask = np.zeros(self.imat.shape,dtype=bool)
    bool_mask[tuple(light_spots.T)] = True
    bool_mask[tuple(dark_spots.T)] = True
    seed_mask,num_seeds = ndimage.label(bool_mask)
    ws = morphology.watershed(blurred,seed_mask)
    background = max(set(ws.ravel()), key=lambda g: np.sum(ws == g))
    background_mask = (ws == background)
    cleaned = self.im * ~background_mask
    if reverse:
      reverse_fill = np.full(cleaned.shape,self.bg)
      self.imatClean = np.asarray(np.abs(cleaned-reverse_fill),dtype=np.uint8)
    else:
      self.imatClean = cleaned - background_mask
    self.imClean = Image.fromarray(self.imatClean)

  def getHicut(self,nbin=50):
    """
      Helper method to remove background
      nbin: Number of bins to make color histogram
      returns: cutoff for higher pixel
    """
    hist, bin_edges = np.histogram(self.imat,bins=nbin,density=True)
    x = np.mean(np.dstack([bin_edges[:-1],bin_edges[1:]]).reshape((nbin,2)),axis=1)
    y = hist
    xy = np.dstack([x,y]).reshape((nbin,2))
    which = np.argmax(np.abs(np.diff(y)))
    if which >= 1: which = which - 1
    return x[which]

  def plotImage(self,plotClean=True):
    """
      To plot Image
      plotClean: True to plot clean image (background substracted) and 
                 False to plot origin image
    """
    if plotClean:
      plt.imshow(self.imatClean,cmap=plt.cm.gray)     
    else:
      plt.imshow(self.imat,cmap=plt.cm.gray)
    plt.show()

  def setAFMobject(self):
    """
      To set object embedded in image
    """
    self.ao = self.getAFMobject(self.imatClean)

  def getAFMobject(self,inputMat):
    """
      Get molecular object embedded in the image
      inputMat: Input matrix (preferrably background substracted)
      returns: An object of class AFMObject
    """
    ao = AFMObject(inputMat)    
    ao.findObjectInImage(background=self.bg)
    ao.setObjectCenter()
    ao.setObjectSize()
    return ao

  def crop2object(self):
    """
       Crop image only to molecular object
    """
    ao = self.ao
    self.imatCrop = ao.crop2object(self.imatClean)
    self.imCrop = Image.fromarray(self.imatCrop)

  def getAxis(self,crop=True):
    """
      Get 2D principal component axis along molecular
      object in the image
      crop: True use cropped image to get molecular
            object, else use uncropped molecular object
    """
    if crop:
      mo = self.getAFMobject(self.imatCrop)
    else:
      mo = self.ao
    mo.setObjectShape()
    mo.setAxis()
    self.co = mo
    self.major = mo.major
    self.minor = mo.minor
    self.cropOption = crop

  def showObjectWithAxis(self,crop=True):
    """
      Plot object along with 2D axes
      crop: True use image from cropped molecular object, else
            use image from uncropped molecular object
    """
    if self.cropOption != crop:
      raise Exception("Argument crop mismatch to the previous option")
    if crop:      
      self.co.showObjectWithAxis(self.imatCrop)
    else:
      self.ao.showObjectWithAxis(self.imatClean)

class AFMObject:
  """
    To work with molecular object embedded in AFM image
  """
  def __init__(self,imageMat):
    """
      Constructor to set input image matrix
      imageMat: input image matrix
    """
    self.imageMat = imageMat
    self.setShape()

  def setShape(self):
    """
      Set height and width of the image
    """
    # N.B. in matrix representation shape gives height by width
    self.height, self.width = self.imageMat.shape

  def findObjectInImage(self,background=255):
    """
      To find molecular object in image. The molecular object must have
      pixels less than background. So image mode up_is_dark is assumed.
      background: The background value
    """
    w = self.width; h = self.height
    omat = np.zeros((h,w))
    for i in range(h):
      for j in range(w):
        omat[i,j] = self.imageMat[(i, j)] < background
    self.ObjectMat = omat

  def setObjectCenter(self):
    """
      Set center of the object
    """
    self.center = self.findObjectCenter(self.ObjectMat)

  def findObjectCenter(self,mat):
    """
      To find center of the molecular object in image
      mat: input matrix
      returns: center along width and height
    """
    h = mat.shape[0]
    w = mat.shape[1]
    mat = mat / np.sum(np.sum(mat))
    dh = np.sum(mat, axis=1)
    dw = np.sum(mat, axis=0)
    cw = np.sum(dw * np.arange(w))
    ch = np.sum(dh * np.arange(h))    
    return (cw,ch)

  def setObjectSize(self):
    """
      Sets size and corners values of the object
    """
    self.size, self.corners = self.findObjectSize(self.ObjectMat)

  def findObjectSize(self,ObjectMat):    
    """
      To find size of the object and corners
      ObjectMat: input matrix containing the object
      returns: a tuple of size and corners (topleft, topright, 
               bottomright,bottomleft)
    """
    ry,rx = np.nonzero(ObjectMat)
    topleft = (min(ry),min(rx))
    botrigh = (max(ry),max(rx))
    toprigh = (min(ry),max(rx))
    botleft = (max(ry),min(rx))
    corners = np.array([topleft, botleft, toprigh, botrigh]).reshape(4,2)
    dy = np.max(corners[:,0]) - np.min(corners[:,0])
    dx = np.max(corners[:,1]) - np.min(corners[:,1])
    size = np.array([dy,dx])
    return (size,corners)

  def setObjectShape(self):
    """
      Set shape of the object
    """
    self.shape = self.findObjectShape(self.ObjectMat,self.center)

  def findObjectShape(self,ObjectMat,center):
    """
      To find shape of the object. Shape indicates a 2D covariance
      matrix, radii along width and height, and rotation matrix 
      corresponding to the covariance matrix. This 3-length tuple is
      indicated as shape.
      ObjectMat: input matrix containing the object
      center: Center of the object
    """
    rx,ry = np.nonzero(ObjectMat)
    n = len(rx)
    xy = np.dstack([rx,ry])
    cx = center[0]; cy = center[1];
    mean = [cy,cx]
    xy = xy.reshape((n,2))    
    dd = (xy - mean).T
    covar = np.dot(dd, dd.T) / n
    u,s,rotation = np.linalg.svd(covar)
    radii = np.sqrt(s)
    return (covar,radii,rotation)

  def crop2object(self,inputMat):
    """
      Crop to object
      inputMat: input matrix
      returns: modified form of input matrix with modified shape
    """
    corners = self.corners
    fromI = corners[0,0]; fromJ = corners[1,1];
    toI = np.min([corners[3,0]+1,inputMat.shape[0]])
    toJ = np.min([corners[2,1]+1,inputMat.shape[1]])
    return inputMat[fromI:toI,fromJ:toJ]

  def setAxis(self):
    """
      Set 2D principal axis of object in image
    """
    self.major, self.minor = self.findAxis(self.shape,self.center)

  def findAxis(self,shape,center,scaleFactor=1.0):
    """
      To find 2D principal axis along object in image
      shape: input shape. See documentation for findObjectShape
      center: center of object
      scaleFactor: scaling along 2D plane to set length of major and 
                   minor axes
      returns: major and minor axes
    """
    radii = shape[1]
    rotation = shape[2]
    xy1 = self.findLine(rotation[:,0],radii[0],center,scaleFactor=scaleFactor)
    xy2 = self.findLine(rotation[:,1],radii[1],center,scaleFactor=scaleFactor)
    major = np.asarray([xy1[1],xy1[0]])
    minor = np.asarray([xy2[1],xy2[0]])
    return (major,minor)

  def findLine(self,vec,radii,center,scaleFactor=1.0):
    """
      A helper method for findAxis 
      vec: 2D vector (principal component axis)
      radii: a principal component value
      center: center of the object
      scaleFactor: scaling along 2D plane to set length of major and 
                   minor axes
    """
    vec = vec*radii*scaleFactor;
    cx = center[0]; cy = center[1];
    mean = [cy,cx]
    pt = np.dstack([mean - vec,mean + vec])
    x = pt[:,0].ravel(); y = pt[:,1].ravel();
    return (x,y)

  def showObjectWithAxis(self,inputMat):
    """
      Plot object along with 2D axes
      inputMat: input matrix containing object
    """
    plt.imshow(inputMat)
    plt.plot(self.center[0],self.center[1],'go');
    plt.plot(self.major[0,:],self.major[1,:]); plt.plot(self.minor[0,:],self.minor[1,:])
    plt.show()

class ImageProcessOption:
  """
    A container for different options for image processing
  """
  def __init__(self,fromarray=False,mode='up_is_dark',
               hilo=(None,3),reverse=False,nbin=50,bgcorr=True,normalize=True,
               plotClean=True,crop=True,
               bg=255,bgfill=255,measure='pearson'):
    """
      Constructor to set different options
      fromarray: Whether to read image from matrix or file
      mode: Either 'up_is_dark' (higher z-values are darker than lower 
            z-value) or 'up_is_light' (higher z-values are lighter than 
            lower z-value)
      hilo: A tuple of high and low pixel value to auto-detect
            background
      reverse: Either True to make an image with mode up_is_light
               or False to make an image with mode up_is_dark
      nbin: Number of bins to create color histogram
      bgcorr: Whether to perform background corrections (removal)
      normalize: Whether to normalize the matrix data to be within 0 and 
                 background value
      plotClean: Whether to plot cleaned matrix or original matrix
      crop: Whether to use cropped matrix or original matrix
      bg: The background used for correction or normalize
      bgfill: A filling value to the background to bias correlation 
              coefficient to use more object detail
      measure: A measure of similarity to compare two images. Can be
               'pearson', 'l1norm' or 'l2norm'
    """
    self.fromarray = fromarray
    self.mode = mode
    self.hilo = hilo
    self.reverse = reverse
    self.nbin = nbin
    self.plotClean = plotClean
    self.crop = crop
    self.normalize = normalize
    self.bg = bg
    self.bgcorr = bgcorr
    self.bgfill = bgfill
    self.measure = measure

class ImageCompare:
  """
    For comparison of two images
  """
  verbose = 0
  def __init__(self,target,reference,               
               targetProcessOption=ImageProcessOption(),
               referenceProcessOption=ImageProcessOption(),
               postprocess=True,superpose=False,
               verbose=0):
    """
      Constructor method to set target and reference data; and
      also sets different options to process target and reference
      data
      target: target data
      reference: reference data
      targetProcessOption: Set of options to process target data 
                           of class ImageProcessOption
      referenceProcessOption: Set of options to process reference data
                              of class ImageProcessOption
      postprocess: Complete postprocessing using ImageProcessOption instances
      superpose: Whether to superpose two images before finding similarity
      verbose: Verbosity
    """
    self.verbose = verbose

    # Save options
    self.targetProcessOption = targetProcessOption
    self.referenceProcessOption = referenceProcessOption

    # Use options to create AFMInputImage instances
    targetMode = targetProcessOption.mode
    refernMode = referenceProcessOption.mode
    targetFromArray = targetProcessOption.fromarray
    refernFromArray = referenceProcessOption.fromarray
    targetBG = targetProcessOption.bg
    refernBG = referenceProcessOption.bg
    if self.targetProcessOption.crop != self.targetProcessOption.crop:
      raise Exception("Incompatible options")
    targetNormalize = targetProcessOption.normalize
    refernNormalize = referenceProcessOption.normalize

    print(target); 
    print(targetFromArray); 
    print(targetMode); 
    print(targetBG);
    print(targetNormalize)
    print(reference); 
    print(refernFromArray); 
    print(refernMode); 
    print(refernBG);
    print(refernNormalize)
    print("All good here")

    targetImage = AFMInputImage(target,fromarray=targetFromArray,mode=targetMode,
                                bg=targetBG,normalize=targetNormalize,verbose=verbose)
    refernImage = AFMInputImage(reference,fromarray=refernFromArray,mode=refernMode,
                                bg=refernBG,normalize=refernNormalize,verbose=verbose)   
    self.targetImage = targetImage
    self.refernImage = refernImage            

    # If required postprocess here directly
    if postprocess: self.do_postprocess(superpose=superpose)

  def do_postprocess(self,superpose=False):
    """
      To perform all postprecessing on input images
      superpose: Whether to superpose two images before finding similarity
    """

    # Background correction of target image if asked
    if self.targetProcessOption.bgcorr:
      reverse = self.targetProcessOption.reverse
      hilo = self.targetProcessOption.hilo
      nbin = self.targetProcessOption.nbin
      self.targetImage.rm_background(reverse=reverse,hilo=hilo,nbin=nbin)
    # Crop target image if asked
    crop = self.targetProcessOption.crop
    if crop:    
      self.targetImage.setAFMobject(); self.targetImage.crop2object();
      self.targetImage.getAxis(crop=crop)
      if not superpose:
        imT = self.targetImage.imCrop
    else:
       imT = self.targetImage.imClean

    
    # Background correction of reference image if asked
    if self.referenceProcessOption.bgcorr:
      reverse = self.referenceProcessOption.reverse
      hilo = self.referenceProcessOption.hilo
      nbin = self.referenceProcessOption.nbin
      self.refernImage.rm_background(reverse=reverse,hilo=hilo,nbin=nbin)

    # Crop reference image if asked
    crop = self.referenceProcessOption.crop
    if crop:
      if self.verbose > 0: print("Crop is used")
      self.refernImage.setAFMobject(); self.refernImage.crop2object();
      self.refernImage.getAxis(crop=crop)      
      
      if superpose:

        # superpose images if asked
        if self.verbose > 0: print("Correlation with scaling and superposition")            
        self.targetMat,self.referenceMat = self.superposeObject(self.targetImage,self.refernImage)

      else:
        # scale objects using their shape information
        if self.verbose > 0: print("Correlation without superposition but with scaling")
        self.targetMat,self.referenceMat = self.smartscaleObject(self.targetImage,self.refernImage)          

    else:

      if self.verbose > 0: print("Crop not used")
      if self.verbose > 0: print("Correlation with simple resizing")
      imR = self.refernImage.imClean 
           
      # naive scaling
      self.targetMat,self.referenceMat = self.scaleObject(imT,imR)

    # Fill bakcground with some color if required    
    tbgfill = self.targetProcessOption.bgfill
    rbgfill = self.referenceProcessOption.bgfill
    if tbgfill != 255.0:
      self.targetMat = self.BGfiller(self.targetMat,tbgfill)
    if rbgfill != 255.0:
      self.referenceMat = self.BGfiller(self.referenceMat,rbgfill)    

  def BGfiller(self,mat,fill):
    """
      Helper method of fill input matrix background part. 
      Note: For this, either the matrix needs to background removed or matrix
            contains an object with no background. mode 'up_is_dark' is 
            assumed
      mat: Input matrix
      fill: Fill value
      returns: modified matrix
    """
    mat[mat == 255.0] = fill
    return mat
 
  def smartscaleObject(self,imTObj,imRObj):
    """
      Scale target and reference object based on their shape
      imTObj: Target AFM object of class AFMInputImage
      imRObj: Reference AFM object of class AFMInputImage
      returns: a tuple of matrices from target and reference object
    """
    Tobj = imTObj.co; Robj = imRObj.co;
    tcen = np.array(Tobj.center); rcen = np.array(Robj.center);
    th,tw = Tobj.size; rh,rw = Robj.size;
    maxH = max((th,rh)); maxW = max((tw,rw));    
    icen = np.array((maxW / 2.0, maxH / 2.0))
    tcorn_shifted = np.round((Tobj.corners - (tcen - icen)))
    rcorn_shifted = np.round((Robj.corners - (rcen - icen)))
    all_corners = np.vstack((tcorn_shifted,rcorn_shifted))    
    all_corners = all_corners - np.min(all_corners,axis=0)
    size_new = np.max(all_corners,axis=0) - np.min(all_corners,axis=0)
    new_maxH, new_maxW = np.asarray(np.round(size_new),dtype=np.int)
    im1 = Image.new('L',(new_maxW,new_maxH),color=255)
    im2 = Image.new('L',(new_maxW,new_maxH),color=255)
    xoff,yoff = np.asarray(np.round(all_corners[0,:]),dtype=np.int)
    im1.paste(imTObj.imCrop,(xoff,yoff))
    xoff,yoff = np.asarray(np.round(all_corners[4,:]),dtype=np.int)
    im2.paste(imRObj.imCrop,(xoff,yoff))
    mat1 = np.asarray(im1.getdata(),dtype=np.float64).reshape((im1.size[1],im1.size[0]))
    mat2 = np.asarray(im2.getdata(),dtype=np.float64).reshape((im2.size[1],im2.size[0]))
    return (mat1,mat2)

  def scaleObject(self,imT,imR):
    """
      Naive scaling of input matrices
      imT: target input matrix
      imR: reference input matrix
      returns: a tuple of matrices from modified target image matrix 
               and reference image matrix
    """
    sx,sy = imR.size   
    imTemp = imT.resize((sx,sy))
    matTemp = np.asarray(imTemp.getdata(),dtype=np.float64).reshape((imTemp.size[1],imTemp.size[0]))
    matR = np.asarray(imR.getdata(),dtype=np.float64).reshape((imR.size[1],imR.size[0]))    
    return (matTemp,matR)

  def rotate_with_expand(self,im,theta):
    """
      Helper method for superposing target and reference object
      im: Input image
      theta: Degree to rotate
      returns: modified image rotated and expanded
    """
    im2 = im.convert('RGBA')
    rot = im2.rotate(theta, expand=True)
    fff = Image.new('RGBA', rot.size, (255,)*4)
    out = Image.composite(rot, fff, rot)
    return out.convert(im.mode)

  def superposeObject(self,imTObj,imRObj):
    """
      Method to superpose target image to reference image
      imTObj: target image object of class AFMInputImage
      imRObj: reference image object of class AFMInputImage
      returns: a tuple with modified matrices
    """
    Tobj = imTObj.co; Robj = imRObj.co;
    tcen = np.array(Tobj.center); rcen = np.array(Robj.center);
    th,tw = Tobj.size; rh,rw = Robj.size;
    maxH = max((th,rh)); maxW = max((tw,rw)); icen = np.array((maxW / 2.0, maxH / 2.0))
    tvec = rcen - tcen; Tw = tvec[0]; Th = tvec[1];
    mg = Robj.shape[1][0] / Tobj.shape[1][0]
    sw = tw * mg; sh = th * mg;
    tv1 = Tobj.shape[2][:,0]; rv1 = Robj.shape[2][:,0];
    theta = np.arccos(np.dot(tv1,rv1))
    im_rot = self.rotate_with_expand(imTObj.imCrop,theta * (180.0 / np.pi))        
    matTemp = np.asarray(im_rot.getdata(),dtype=np.float64).reshape((im_rot.size[1],im_rot.size[0]))
    imtarget = AFMInputImage(matTemp,fromarray=True,mode='up_is_dark',bg=255,normalize=False)
    imtarget.setAFMobject(); imtarget.crop2object();
    imtarget.getAxis(crop=True)
    mat1, mat2 = self.smartscaleObject(imtarget,imRObj)
    return (mat1,mat2)

  def getCorrelation(self,overlapOnly=False,unionOnly=False,q=50):
    """
      Calculate similarity measure between target and reference matrices 
      (by subsetting if asked)
      Note: The method of subsettings requires background is 255.0 in both the
            images of target and reference
      overlapOnly: A method of subsetting that takes intersection of target
                   and reference matrices                   
      unionOnly: A method of subsetting that takes union of target and 
                 reference matrices
                 Note: unionOnly takes precedence over overlapOnly 
      q: Percentile value used to fill missing values (backgrounds) in 
         unionOnly case
      If both overlapOnly and unionOnly are False, no subsetting is used
      returns: similarity measure
    """
    if q is None: q = 50.0
    if overlapOnly and unionOnly: overlapOnly = False
    A = []
    B = []
    if overlapOnly:
      for i in range(self.targetMat.shape[0]):
        for j in range(self.targetMat.shape[1]):
          if self.targetMat[i,j] < 255.0 and self.referenceMat[i,j] < 255.0:
             A.append(self.targetMat[i,j]); B.append(self.referenceMat[i,j]);
      if len(A) > 0: 
        cc = self.get_measure(A,B)
      else: cc = -1.0
    elif unionOnly:
      if q < 0. or q > 100.: q = 50.
      tb = self.targetMat < 255.0
      rb = self.referenceMat < 255.0
      nA = len(self.targetMat[tb]); nB = len(self.referenceMat[rb]);
      if q != 50.:
        tmedian = np.percentile(self.targetMat[tb],q)
        rmedian = np.percentile(self.referenceMat[rb],q)
      else:
        tmedian = np.median(self.targetMat[tb])
        rmedian = np.median(self.referenceMat[rb])     
      for i in range(self.targetMat.shape[0]):
        for j in range(self.targetMat.shape[1]):
          if self.targetMat[i,j] < 255.0 or self.referenceMat[i,j] < 255.0:
             if self.targetMat[i,j] < 255.0 and self.referenceMat[i,j] == 255.0:
               A.append(self.targetMat[i,j]); B.append(rmedian);
             elif self.referenceMat[i,j] < 255.0 and self.targetMat[i,j] == 255.0:
               A.append(tmedian); B.append(self.referenceMat[i,j]);
             else:
               A.append(self.targetMat[i,j]); B.append(self.referenceMat[i,j]);
      if (nA + nB) > len(A):
        cc = self.get_measure(A,B)
      else: cc = -1.0
    else:  
      A = self.targetMat.ravel()
      B = self.referenceMat.ravel()            
      cc = self.get_measure(A,B)
    return cc

  def get_measure(self,v1,v2):
    """
      Use different similarity measure metric to compare input vectors
      v1: First input vector
      v2: Second input vector
      returns: similarity measure
    """
    meth = self.referenceProcessOption.measure
    if meth == 'pearson':
      val,_ = pearsonr(v1,v2)
    elif meth == 'l1norm':
      v12 = v1 - v2
      val = np.sum(np.abs(v12))
    elif meth == 'l2norm':
      v12 = v1 - v2
      val = np.sqrt(np.sum(v12*v12))
    else:
      val,_ = pearsonr(v1,v2)
    return val

  def plotCompareImage(self,filename,overlapOnly=False,unionOnly=False,q=50):
    """
      Plot the images that were compared side-by-side
      filename: Output image file
      overlapOnly, unionOnly, q: See documentation for getCorrelation
    """
    if q is None: q = 50.0
    if overlapOnly and unionOnly: overlapOnly = False      
    if overlapOnly:
      A = np.full(self.targetMat.shape,255.0)
      B = np.full(self.referenceMat.shape,255.0)
      for i in range(self.targetMat.shape[0]):
        for j in range(self.targetMat.shape[1]):
          if self.targetMat[i,j] < 255.0 and self.referenceMat[i,j] < 255.0:
            A[i,j] = self.targetMat[i,j]; B[i,j] = self.referenceMat[i,j];
    elif unionOnly:
      if q < 0. or q > 100.: q = 50.
      tb = self.targetMat < 255.0
      rb = self.referenceMat < 255.0
      if q != 50.:
        tmedian = np.percentile(self.targetMat[tb],q)
        rmedian = np.percentile(self.referenceMat[rb],q)
      else:
        tmedian = np.median(self.targetMat[tb])
        rmedian = np.median(self.referenceMat[rb]) 
      A = np.full(self.targetMat.shape,255.0)
      B = np.full(self.referenceMat.shape,255.0)
      for i in range(self.targetMat.shape[0]):
        for j in range(self.targetMat.shape[1]):
          if self.targetMat[i,j] < 255.0 or self.referenceMat[i,j] < 255.0:
            if self.targetMat[i,j] < 255.0 and self.referenceMat[i,j] == 255.0:
               A[i,j] = self.targetMat[i,j]; B[i,j] = rmedian;
            elif self.referenceMat[i,j] < 255.0 and self.targetMat[i,j] == 255.0:
               A[i,j] = tmedian; B[i,j] = self.referenceMat[i,j];
            else:
               A[i,j] = self.targetMat[i,j]; B[i,j] = self.referenceMat[i,j];
    else:
      A = self.targetMat
      B = self.referenceMat
    f, (ax1, ax2) = plt.subplots(1, 2)    
    ax1.imshow(A,cmap=plt.cm.gray,vmin=0,vmax=255);
    ax2.imshow(B,cmap=plt.cm.gray,vmin=0,vmax=255);
    f.savefig(filename)
    plt.close()

class EllipsoidPairDistRestraint:
  """
    Ellipsoid pairwise center to center distance restraint
  """
  def __init__(self,ellipsoids,Pair,strength=1000.0):            
    """
      Constructor to set initial pairwise distances from set
      of ellipsoids and all pairs
      ellipsoids: An array of ellipsoids
      Pair: The list of pairs 
      strength: Force-constant
    """
    self.ok = np.full(len(Pair),0)
    self.ok[:] = False
    self.native_distances = self.get_pairs_dist(ellipsoids,Pair)  
    self.Pair = Pair
    self.strength = strength

  def get_pair_dist(self,e1,e2):
    """
      Helper for distance between centers calculation
      e1,e2: First and second ellipsoids
      returns: distance between centers of ellipsoids
    """
    coord = np.zeros((2,3))
    coord[0,:] = e1.mean; coord[1,:] = e2.mean;
    return sqrt(np.sum((coord[0,:] - coord[1,:])**2))

  def get_pairs_dist(self,ellipsoids,Pair):
    """
      Get list of distances between centers for input pair lists of ellipsoids
      ellipsoids: Array of ellipsoids
      Pair: list of pairs
    """
    ne = len(ellipsoids); distances = []
    for ip,p in enumerate(Pair):
      i = p[0]; j = p[1];
      if self.ok[ip]:
        dist = self.get_pair_dist(ellipsoids[i],ellipsoids[j])
      else:
        if i < ne and j < ne:
          dist = self.get_pair_dist(ellipsoids[i],ellipsoids[j])
          self.ok[ip] = True
        else:
          dist = None
      distances.append(dist)
    return np.array(distances)

  def set_input(self,input):
    """
      Set ellipsoids to the input
      input: Input array of ellipsoids
    """
    self.ellipsoids = input

  def get_measure(self):
    """
      Get the value of restraint score
      returns: restraint score
    """
    curr_distances = self.get_pairs_dist(self.ellipsoids,self.Pair)
    delta = (curr_distances - self.native_distances)
    delta = delta * delta * self.strength    
    return np.sum(delta)

  def evaluate(self,inp):
    """
      Evaluate restraint score given input ellipsoids
      inp: Input array of ellipsoids
    """
    self.set_input(inp)
    return self.get_measure()

class EllipsoidExcludedVolumeRestraint:
  """ 
    Ellipsoid excluded volume restraint
  """
  def __init__(self,ellipsoids,strength=np.inf,cutoff=0.01):
    """
      Constructor to set initial correlation from set
      of ellipsoids
      ellipsoids: An array of ellipsoids
      strength: Force-constant per pair
      cutoff: A threshold to apply input strength
    """    
    ne = len(ellipsoids)
    self.native_correlation = []
    self.native_correlation = self.get_corr_pairs(ellipsoids,set_self=True)
    self.strength = strength
    self.cutoff = cutoff

  def set_cc_self(self,cc):
    """
      Set self correlations
      cc: self-correlation
    """
    self.cc_self = cc

  def get_cc_self(self):
    """
      Get self correlation
      returns: self-correlation
    """
    return self.cc_self

  def _create_input_gauss(self,e):
    """
      Helper to make a simple dictionary with center and
      covariance matrix of an ellipsoid
      e: input ellipsoid
      returns: a dictionary with keys 'mu' including center
               and 's' including covariance matrix
    """
    center = np.array(e.mean); covmat = e.covar;
    return {'mu': center, 's': covmat}  

  def get_corr_pair(self,g1,g2):
    """
      Calculate overlap between input Gaussians
      g1,g2: First and second Gaussian input to class Overlap
    """
    ov = Overlap(g1,g2)
    return ov.get_overlap_gaussian()

  def get_corr_pairs(self,ellipsoids,set_self=True):
    """
      Calculate pairwise overlaps
      ellipsoids: Input array of ellipsoids
      set_self: Whether to save current ellipsoids self correlations
      returns: an array of correlations between all pairs of Gaussians
    """
    ne = len(ellipsoids)
    if set_self:
      cc_self = []
      for i in range(ne):
        g = self._create_input_gauss(ellipsoids[i])
        cc_self.append(np.sqrt(self.get_corr_pair(g,g)))
      self.set_cc_self(cc_self)
    else:
      cc_self = self.get_cc_self()
    cc = []
    for i in range(ne):
      g1 = self._create_input_gauss(ellipsoids[i])
      cc_ii_deno1 = cc_self[i]
      for j in range(i+1,ne):
        g2 = self._create_input_gauss(ellipsoids[j])
        cc_jj_deno2 = cc_self[j]
        cc_ij = self.get_corr_pair(g1,g2)
        cc.append(cc_ij/(cc_ii_deno1*cc_jj_deno2))
    return np.array(cc)

  def set_input(self,input):
    """
      Set input ellipsoids
      input: Input array of ellipsoids
    """
    self.ellipsoids = input

  def get_measure(self):
    """
      Get excluded volume restraint score
      returns: Excluded volume restraint score
    """
    curr_correlation = self.get_corr_pairs(self.ellipsoids)
    #delta = 100.0 * (np.abs(curr_correlation - self.native_correlation)/self.native_correlation)
    #print(delta)
    #value = int(np.any(delta > self.cutoff)) * self.strength
    for i,c in enumerate(curr_correlation): print("Current Correlation: %d %e" % (i,c))
    print("Max overlap correlation: %e" % (np.max(curr_correlation)))   
    value = int(np.any(curr_correlation > self.cutoff)) * self.strength
    if np.isnan(value): value = 0
    return value

  def evaluate(self,inp):
    """
      Evaluate excluded volume restraint score
      inp: Input array of ellipsoids
      returns: Excluded volume restraint score
    """
    self.set_input(inp)
    return self.get_measure()

class ConstantRestraint(object):
  """
    Simple restraint container to return constant score
  """
  def __init__(self,value=1.0):
    """
      Constructor to set constant score
    """
    self.value = value

  def evaluate(self,inp):
    """
      Evaluate constant score
    """
    return self.value

class pcFit:
  """
    For fitting point-cloud data with known correspondance
  """
  def __init__(self,target,source,index=None,algo='both'):
    """
      Sets target and source so that by fitting source moves towars target
      target: Target matrix with N x 3 shape including N atoms
      source: Source matrix with M x 3 shape including M atoms
      index: An pairlist to give correspondance between N atoms to M atoms
             By default, it is assumed that N = M and i-th atom from target
             matrix corresponds to i-th atom from source matrix
      algo: The algorithm to fit. The choices are, 'trans' to fit by
            translation only, 'rot' to fit rotation only, 'both' to fit by
            both translation and rotation, 'transx' to fit by translation
            along x, 'transy' to fit by translation along y, 'transz' to fit
            by translation along z
    """
    self.sourceFull = source
    if index is None:
      n = np.min([len(target),len(source)])
      seqI = np.array([i for i in range(n)])
      index = np.vstack([seqI,seqI]).T        
    self.target = []; self.source = [];
    for i in range(index.shape[0]):
      tindex = index[i,0]; sindex = index[i,1];
      if tindex >= 0 and tindex < len(target) and sindex >= 0 and sindex < len(source):
        self.target.append(target[tindex]); self.source.append(source[sindex]);
    self.target = np.array(self.target)
    self.source = np.array(self.source)
    self.algo = algo

  def fit(self):
    """
      Perform fitting
      returns: Modified source matrix
    """
    param = self.get_fit_param()
    return self.apply_fit_param(self.sourceFull,param)

  def get_fit_param(self):
    """
      Get fitting parameter
      returns: An instance of class Tranformation including optimum rotation 
               and translation depending on the algo
    """
    if self.algo == 'both':
      translation = self.get_translation()
      rotation = self.get_rotation()      
    elif self.algo == 'trans':
      translation = self.get_translation()
      rotation = np.identity(3)
    elif self.algo == 'rota':
      translation = np.array([0.,0.,0.])
      rotation = self.get_rotation()
    elif self.algo == 'transx':
      translation = self.get_translation()
      translation[1] = 0.0; translation[2] = 0.0;
      rotation = np.identity(3)
    elif self.algo == 'transy':
      translation = self.get_translation()
      translation[0] = 0.0; translation[2] = 0.0;
      rotation = np.identity(3)
    elif self.algo == 'transz':
      translation = self.get_translation()
      translation[0] = 0.0; translation[1] = 0.0;
      rotation = np.identity(3)
    tr = Transformation(trans=translation,rot=rotation,rotationType='m')
    return tr

  def get_translation(self):
    """
      Get translation for target matrix
      returns: center of target matrix
    """
    return np.mean(self.target,axis=0)

  def get_rotation(self):    
    """
      Get rotation matrix to rotate source to target after putting both
      target and source at the origin (0,0,0)
      returns: rotation matrix (3x3)
    """
    rotation = np.identity(3)
    translation = np.array([0.,0.,0.])
    dummy_tr = Transformation(trans=translation,rot=rotation,rotationType='m')
    target_origin = self.apply_fit_param(self.target,dummy_tr)    
    source_origin = self.apply_fit_param(self.source,dummy_tr)
    rotmat = self.get_rotation_matrix(target_origin,source_origin)
    return rotmat

  def get_rotation_matrix(self,mat1,mat2):
    """
      Given two matrix use Diamonds method based on quaternions to 
      superpose 2nd matrix to 1st matrix
      mat1,mat2: First and second input matrix
      returns: Optimum rotation matrix (3x3)
    """
    M = np.dot(mat2.T,mat1)
    TraceM = np.sum(np.diag(M))
    I = np.diag(np.repeat(2*TraceM,3))
    Q = M + M.T - I
    V = np.array([M[1,2] - M[2,1],M[2,0] - M[0,2],M[0,1] - M[1,0]])
    P = np.zeros((4,4))
    P[0:3,0:3] = Q; P[0:3,3] = V; P[3,0:3] = V
    val,vec = np.linalg.eig(P)
    val_index = np.argsort(-val)
    val = val[val_index]; vec = vec[:,val_index];
    l,m,v,s = vec[:,0]
    params = []
    params.append( l*l - m*m - v*v + s*s); params.append(2.0*(l*m - v*s)); params.append(2.0*(l*v + m*s));
    params.append(2.0*(l*m + v*s)); params.append(-l*l + m*m - v*v + s*s); params.append(2.0*(m*v - l*s));
    params.append(2.0*(l*v - m*s)); params.append(2.0*(m*v + l*s)); params.append(-l*l - m*m + v*v + s*s);
    return np.array(params).reshape(3,3).T

  def apply_fit_param(self,mat,tr):
    """
      Apply fitting parameters on input matrix
      mat: Input matrix
      tr: An instance of class Transformation including optimum rotation and
          translation depending on the algo
      returns: fitted matrix
    """
    trans = tr.get_translation()
    rot = tr.get_rotation()
    mat_cen = np.mean(mat,axis=0)
    if self.algo == 'transx': 
      mat_cen[1] = 0.0; mat_cen[2] = 0.0;
    elif self.algo == 'transy':
      mat_cen[0] = 0.0; mat_cen[2] = 0.0;      
    elif self.algo == 'transz':
      mat_cen[0] = 0.0; mat_cen[1] = 0.0;
    mat_origin = mat - mat_cen
    return mat_origin.dot(rot) + trans

class Mover:
  """
    For moving an ellipsoid by translation, rotation
  """
  def __init__(self,ellipsoid,name,index,translation,rotation,scale=0.0):
    """
      Constructor to set ellipsoids and translation and rotation
      parameters
      ellipsoid: An instance of the class ellipsoid
      name: A name string
      translation: Maximum value of translation
      rotation: Maximum value of rotation
      scale: Maximum value of scaling
    """
    self.ellipsoid = ellipsoid
    self.name = name
    self.index = index
    # translation in Angstrom and rotation in radian
    # include "high" exclude "low"
    self.tval = -np.random.uniform(low=-translation,high=0)
    self.rval = -np.random.uniform(low=-rotation,high=0)
    if scale == 0.0:
      self.sval = scale
    else:  self.sval = -np.random.uniform(low=-scale,high=0)
    self.backup()

  def _get_ellipsoid(self):
    """
      Helper method returning saved ellipsoid
      returns: saved ellipsoid
    """
    return self.ellipsoid

  def propose(self,tpropose=None,rpropose=None,spropose=None):
    """
      Proposal to move method. It proposes a translational move to the
      ellipsoid center and a rotational move to the axis of the ellipsoid
      tpropose: Translation proposal. If None generated automatically.
      rpropose: Rotation proposal. If None generated automatically.
      spropose: Scaling proposal. If None generated automatically.
      returns: proposed moved ellipsoid
    """
    e = copy.copy(self.ellipsoid)
    if tpropose is None:
      tpropose = self.propose_center(self.tval)
    if self.tval > 0.0:
      center = e.mean + tpropose
      self.dist = np.linalg.norm(center-e.mean)
    else:
      center = e.mean
      self.dist = 0.0
    if rpropose is None:
      rpropose = self.propose_rotation(self.rval)
    if self.rval > 0.0:
      q1 = e.get_rotation4D()
      q2 = rpropose
      rotatedQ = self.qmult(q1,q2)    
      trf = Transformation(rot=rotatedQ,rotationType='q',trans=0.0)
      rotatedM = trf.change_rotation(fromtype='q',totype='m')
    else:
      rotatedM = e.rotation
    if spropose is None:
      spropose = self.propose_scale(self.sval)
    if self.sval > 0.0:
      radii = e.radii + spropose
      e.set_radii(radii,scaled=False)      
    covar = invert_decomposition(rotatedM,e.radii)
    e.set_mean(center); e.set_covar(covar)
    # Following fixes covariance matix if there is any error
    e.calculate_aux()
    self.ellipsoid = copy.copy(e)
    return self._get_ellipsoid()

  def generate_point_on_unit_sphere(self):
    """
      Helper for random translation and rotation generator
      returns: a point on unit sphere
    """
    vec = np.random.randn(3,1).flatten()
    vec /= np.linalg.norm(vec)
    return vec
  def qmult(self,q1,q2):
    """
      Quaternion multiplication
      q1,q2: Input quaternions
      returns: multiplied quaternion
    """
    a1 = q1[0]; b1 = q1[1]; c1 = q1[2]; d1 = q1[3];
    a2 = q2[0]; b2 = q2[1]; c2 = q2[2]; d2 = q2[3];
    qq12 = np.array([a1*a2 - b1*b2 - c1*c2 - d1*d2,
                     a1*b2 + b1*a2 + c1*d2 - d1*c2,
                     a1*c2 + c1*a2 - b1*d2 + d1*b2,
                     a1*d2 + d1*a2 + b1*c2 - c1*b2])
    if qq12[0] < 0.0: qq12 = -qq12
    return qq12

  def propose_center(self,t):    
    """
      Generating random translation of given length
      t: length of the translation
    """
    vec = self.generate_point_on_unit_sphere()    
    return vec * t

  def propose_rotation(self,max_radian):
    """
      Generating random rotation of input magnitude
      max_radian: max degree of rotation
    """
    axis = self.generate_point_on_unit_sphere()
    halfangle = np.random.uniform(low=-max_radian,high=max_radian) * 0.5
    s = np.sin(halfangle)
    a = np.cos(halfangle)
    bcd = axis * s
    self.angle = 2 * halfangle
    return np.insert(bcd,0,a)

  def propose_scale(self,s):
    """
      Gerenating a random scaling along a random direction of ellipsoid
      s: input scaling magnitude
    """
    which_dim = random.sample(xrange(3),1)
    svec = np.zeros(3)
    svec[which_dim] = np.random.uniform(low=-s,high=s)
    return svec

  def accept(self):
    """
      Accept the move
      returns: moved ellipsoid
    """
    return self._get_ellipsoid()

  def reject(self):
    """
      Reject the move
      returns: restored ellipsoid
    """
    self.restore_backup()
    return self._get_ellipsoid()

  def backup(self):       
    """
      Backup ellipsoid instance
    """
    self.ellipsoid_save = copy.copy(self.ellipsoid)
  def restore_backup(self):    
    """
      Restore ellipsoid instance from backup
    """
    self.ellipsoid = copy.copy(self.ellipsoid_save)

class Trajectory:
  """
    To save trajectory of ellipsoids in HDF5 format
  """
  def __init__(self,filename):
    """
      Initialize trajectory by setting filename and file-pointer to write
    """
    self.filename = filename
    print("Creating file: %s" % filename)   
    self.blockCounter = 0

  def add_data(self,data,datatype,first=False):
    if first: mode = 'w'      
    else:     mode = 'a'    
    hf = h5py.File(self.filename,mode) 
    hfpath = "/%s/i%d" % (datatype,self.blockCounter)
    block = hf.create_group(hfpath)
    if datatype == 'Ellipsoid': 
      data_dic = self.extract_ellipsoid_data(data)
    for k in data_dic.keys():
      block.create_dataset(k,data=data_dic[k])
    hf.close()
    self.blockCounter += 1        

  def extract_ellipsoid_data(self,data):
    n = len(data); weight = [None] * n;
    mean = [None] * n; covar = [None] * n;
    radii = [None] * n; scaled_radii = [None] * n;
    rotation = [None] * n; quat = [None] * n
    weight = np.repeat(1.0/n,n);
    mean = np.zeros((n,3)); covar = np.zeros((n,3,3));
    radii = np.zeros((n,3)); scaled_radii = np.zeros((n,3));
    rotation = np.zeros((n,3,3)); quat = np.zeros((n,4));
    for ie,e in enumerate(data):
      mean[ie,:] = e.mean; covar[ie,:,:] = e.covar;
      radii[ie,:] = e.get_radii(scaled=False)
      scaled_radii[ie,:] = e.get_radii(scaled=True)
      rotation[ie,:,:] = e.get_rotation_matrix()
      quat[ie,:] = e.quat 
      weight[ie] = e.mass  
    tw = np.sum(weight)
    for ie in range(len(data)):
      weight[ie] = weight[ie]/tw
    return {'mean':mean,'covar':covar,'radii':radii, 'w': weight,
            'scaled_radii':scaled_radii,'rot': rotation,'q': quat}


class MCSampler:
  """
    For Monte-Carlo sampling by randomly moving ellipsoids and judging 
    the proposal by input scoring functions
  """
  def __init__(self,Options,Data,timing=True):
    """
      Constructor to set misc. information from Option dictionary
      and Data dictionary. The Option dictionary contains extended
      information from input yaml options and default options. Most
      importantly, here options under OPTIMIZE is used. The Data
      dictionary includes ellipsoid representation of the molecule,
      image options and image itset for the target, restraint scoring
      class instance etc. For calling this constructor see documentation
      of method 'sample' from 'Project' class.
      Options: Input Options
      Data: Input Data
      timing: Whether to write timing reports
    """
    self.opt = Options
    self.timing = timing

    # Set random state
    seed = Options['SEED']
    np.random.seed(seed); random.seed(seed);
    self.random_state = random.getstate();
    self.random_state_np = np.random.get_state();

    # Set general private settings
    self.limit_translation = OptionAFM._LIMIT_TRANSLATION_
    self.limit_rotation = OptionAFM._LIMIT_ROTATION_ 
    self.limit_scale = OptionAFM._LIMIT_SCALE_ 

    # copy options
    self.max_translation = Options['OPTIMIZE']['MAX_MOVEMENT']['TRANSLATION']
    self.max_rotation = Options['OPTIMIZE']['MAX_MOVEMENT']['ROTATION']
    self.max_scale = Options['OPTIMIZE']['MAX_MOVEMENT']['SCALE']
    self.kt_mech = Options['OPTIMIZE']['TEMP_MECH']
    self.kt_cc = Options['OPTIMIZE']['TEMP']
    self.wt_factor =  Options['OPTIMIZE']['WEIGHT']
    self.movers_opt = Options['OPTIMIZE']['MOVERS']
    self.npropose = Options['OPTIMIZE']['MOVERS']['NPROPOSE']
    self.total_optimization_steps = Options['OPTIMIZE']['STEPS']
    self.period = Options['OPTIMIZE']['PERIOD_TRJ']           

    # Copy data
    self.ellipsoid_opt = Data['ellipsoid']
    #self.ellipsoid_backup = Data['ellipsoid']
    self.restraints = Data['restr']
    self.refMatClean = Data['refernMat']
    self.f2list = Data['f2list']    

    # Do rest of the processing
    self._process()

  def _process(self):
    """
      Helper for constructor to complete misc. secondary settings
      based on options
    """
    # Set rotation and translation (max. values)
    self.rotation = self.max_rotation
    self.translation = self.max_translation
    self.scale = self.max_scale   

    # Set beta-value for Monte-Carlo based on 3D scoring functions
    self.beta_mech = 1.0 / (self.kt_mech * 8.314)

    # Set anneal schedule if asked
    self.start_step = 0
    self.set_anneal_schedule()
    
    # Initialize movers
    self.initMovers(self.translation, self.rotation, self.scale, 
                    useproj = self.opt["OPTIMIZE"]["MOVERS"]["USEPROJ"])
    
    self.initial_ellipsoid_instance = copy.deepcopy(self.ellipsoid_opt)

    # Calculate initial score
    restr = self.evaluate_restraints()
    if self.opt['MODE'] == 'Debug':
      self.score_mech_curr = self.get_evaluated_scores(restr)
    else: self.score_mech_curr = self.print_evaluated_scores(restr)
    # Reset initial mech score = 0, even if it is an invalid conformation
    if self.score_mech_curr == np.inf: self.score_mech_curr = 0.0      
    print("Initial Score: %15.4f" % self.score_mech_curr)     
        
    # Calculate correlation and score based on image similarity
    self.cc = np.zeros(self.total_optimization_steps)
    self.cc0 = self.get_correlation_from_diffractionImage(verbose=0)
    self.pcc_val_curr = self.cc0
    print("Initial pcc: %.6f" % self.cc0)
    method = self.opt['OPTIMIZE']['BIAS']['METHOD']
    param = self.opt['OPTIMIZE']['BIAS']['PARAM']
    w = self.wt_factor
    self.score_img_curr = self.get_score_on_cc(self.cc0,method=method,param=param,weight=w)

    # Set loop parameters
    self.firstStep = True            
    self.tol = 0.05

    # Initialize drift removal
    self.init_drift_method()

    # Initialize all scores and saves
    self.score_mech_prev = self.score_mech_curr; self.score_img_prev = self.score_img_curr;
    self.score_mech_temp = self.score_mech_curr; self.score_img_temp = self.score_img_curr;
    self.pcc_val_prev = self.pcc_val_curr; self.pcc_val_temp = self.pcc_val_curr;

    # Start writing 'STATFILE'
    op = self.opt['OPTIMIZE']['OUTPUT']
    if op['STAT']:
      baseName = op['STATPATH'] + "/" + op['STATFILE']
      statFile = baseName + "_stat.out"
      with open(statFile,'w') as fp:
        fp.write("%d %.6f %d\n" % (0,self.pcc_val_curr,0))

    # Initialize Trajectory of Ellipsoids
    self.traj = Trajectory(filename=self.opt['OPTIMIZE']['TRJFILE'])
    self.traj.add_data(data=self.ellipsoid_opt[1],
                       datatype='Ellipsoid',first=True)

  def set_anneal_schedule(self):
    """
      Set anneal schedule based on options OPTIMIZE->ANNEAL or simply set
      beta against 3D scoring functions
    """
    if self.opt["OPTIMIZE"]["ANNEAL"]:
      start_step = self.start_step
      self.anneal_schedule = self.read_schedule(self.opt["OPTIMIZE"]["ANNEAL_SCHEDULE_FILE"])
      if self.anneal_schedule is None:
        anneal_from = float(self.opt["OPTIMIZE"]["TEMP_ANNEAL_FROM"])
        anneal_to = float(self.opt["OPTIMIZE"]["TEMP_ANNEAL_TO"])
        anneal_dt = float(self.opt["OPTIMIZE"]["TEMP_ANNEAL_DT"])
        anneal_dn = float(self.opt["OPTIMIZE"]["TEMP_ANNEAL_DN"])
        flag = True; t = anneal_from; step = 0 + start_step;
        self.anneal_schedule = {};
        while flag:        
          self.anneal_schedule[step] = 1.0 / (t  * 8.314)        
          t -= anneal_dt;
          step += anneal_dn
          if t < anneal_to:
            flag = False 
        self.anneal_schedule[self.total_optimization_steps] = 1.0 / (anneal_to  * 8.314)
    else:
      self.beta_cc = 1.0 / (self.kt_cc * 8.314)

  def init_drift_method(self):
    """
      Initialize drift removal based on OPTIMIZE->DRIFT
    """
    self.drift_method = self.opt["OPTIMIZE"]["DRIFT"]    
    if self.drift_method != "none":
      self.remove_drift_level = 1     
      ellipsoids = self.ellipsoid_opt[self.remove_drift_level]
      ecen = []
      for e in ellipsoids:
        ecen.append(np.array(e.mean))
      crd = np.asarray(ecen)
      self.EllipCenInit = crd      
      self.drift_method = self.reset_drift_method(coord=crd,method=self.drift_method)
    print("Removal of drift, method : %s" % self.drift_method)

  def initMovers(self,trans, rot, scale, useproj=False):
    """
      Initilize movers based on maximum values of translation,
      rotation and scaling
      trans: Max translation
      rot: Max rotation in degree
      scale: Max scaling value added
      useproj: Whether to use projection matrix and reference matrix
               difference to update translation and rotation permitted
    """
    
    rotation = (np.pi / 180.0) * rot
    elevels = list(reversed(sorted(self.ellipsoid_opt.keys())))
    self.movers = []
    for elevel in elevels:
      ells = self.ellipsoid_opt[elevel]
      mlist = self.movers_opt[elevel]
      for mindex in mlist:
        m = mindex - 1 # debug
        e = ells[m]  
        name = "EINDEX: " + str(m)         
        self.movers.append(Mover(e,name=name, index = m,
                                   translation=trans,
                                   rotation=rotation,
                                   scale=scale))
      trans = trans / 2
      rotation = rotation / 2 
      scale = scale / 2  
    self.nmover = len(self.movers) 
      
  def updateMovers(self,nreject):
    """
      update mover's movement limits if required
      nreject: number of rejected steps after last accept
    """
    # This is now an useless feature as limit_reject is very high
    limit_nreject = 10**6
    if nreject > limit_nreject:      
      if self.max_translation < self.limit_translation:
        slope = 0.005; y0 = self.max_translation      
        trans = (nreject * slope) + y0
        if trans > self.limit_translation: trans = self.limit_translation
      else:
        trans = self.max_translation
      if self.max_rotation < self.limit_rotation:
        slope = 0.01; y0 = self.max_rotation
        rot = (nreject * slope) + y0
        if rot > self.limit_rotation: rot = self.limit_rotation
      else:
        rot = self.max_rotation
    else:
      trans = self.max_translation; rot = self.max_rotation;
    if self.max_scale > self.limit_scale:
      scale = self.limit_scale
    else:
      scale = self.max_scale    
    self.initMovers(trans, rot, scale,
                    useproj = self.opt["OPTIMIZE"]["MOVERS"]["USEPROJ"])
    self.rotation = rot
    self.translation = trans
    self.scale = scale

  def move_selector(self,method='random',verbose=1):
    """
      Selecting which of the mover to move
      method: 'random' to select mover randomly based on class variable
              npropose. npropose is number of movers to move in one step
              Note: Only random is supported now
      verbose: Verbosity
    """
    if method == 'random':
      self.which = random.sample(xrange(self.nmover),self.npropose)
    else:
      print("\tOnly random selection is implemented, Falling back")
      self.which = random.sample(xrange(self.nmover),self.npropose)
    if verbose > 0:
      print("\tMoving:",end='')
      for w in self.which:
        mv = self.movers[w]        
        print("[%d,%s]" % (w,mv.name),end='')
      print()
      print("\trotation: %.2f translation: %.2f" % (self.rotation,self.translation))
    
  def start(self):
    """
      Start sampling
    """
    print("Performing MC sampling for %d steps" % (self.total_optimization_steps))
    cent = self.opt["OPTIMIZE"]["LIMIT_PROPOSAL_PERCENT"]
    if self.opt['OPTIMIZE']['STOPBY'] == 'accept':
      self.num_propose_limit = (self.total_optimization_steps * 100.0) / cent
    elif self.opt['OPTIMIZE']['STOPBY'] == 'propose':
      self.num_propose_limit = self.total_optimization_steps
    print("  MC number of proposal is restricted to %d" % self.num_propose_limit)
    self.continueSample = True
    self.reset_loop()
    print("Step: %5d [%10d]" % (self.loop,self.stat['number_of_proposed_steps']),end="\n")
    print("Current State: Mec(prev,post): %15.6f,%15.6f Img(prev,post): %15.6f,%15.6f" %
          (self.score_mech_prev,self.score_mech_curr,self.score_img_prev,self.score_img_curr))
    print("---------------------------------------------------------------------------------")
    self.nreject = 0
    self.do_sample()
    self.print_sample_statistics()

  def restart(self,opt):
    """
      Restart sampling after reading STATFILE 
    """
    print("restarting")
    np.random.set_state(self.random_state_np); random.setstate(self.random_state);
    
    statblock = self.readStatFile()
    self.writeStatBlock(statblock,summary=False)
        
    # reset some class variable, according to the new options
    self.max_translation = opt['OPTIMIZE']['MAX_MOVEMENT']['TRANSLATION']
    self.max_rotation = opt['OPTIMIZE']['MAX_MOVEMENT']['ROTATION']
    self.max_scale = opt['OPTIMIZE']['MAX_MOVEMENT']['SCALE']
    self.kt_mech = opt['OPTIMIZE']['TEMP_MECH']
    self.kt_cc = opt['OPTIMIZE']['TEMP']
    self.wt_factor =  opt['OPTIMIZE']['WEIGHT']
    self.movers_opt = opt['OPTIMIZE']['MOVERS']
    self.npropose = opt['OPTIMIZE']['MOVERS']['NPROPOSE']
    self.period = opt['OPTIMIZE']['PERIOD_TRJ'] 
    self.rotation = self.max_rotation
    self.translation = self.max_translation
    self.scale = self.max_scale  
    print("scale = %f" % self.scale) 
    self.beta_mech = 1.0 / (self.kt_mech * 8.314)
    self.start_step = self.total_optimization_steps
    self.total_optimization_steps += opt['OPTIMIZE']['STEPS']
    self.set_anneal_schedule()
    self.opt = opt # Simply set for now, but only few values should be reconfigurable
    self.init_drift_method()
    self.cc = np.append(self.cc,np.zeros(opt['OPTIMIZE']['STEPS']))

    print("Performing MC sampling for %d steps" % (self.total_optimization_steps))
    cent = self.opt["OPTIMIZE"]["LIMIT_PROPOSAL_PERCENT"]
    if self.opt['OPTIMIZE']['STOPBY'] == 'accept':
      self.num_propose_limit = (self.total_optimization_steps * 100.0) / cent
    elif self.opt['OPTIMIZE']['STOPBY'] == 'propose':
      self.num_propose_limit = self.total_optimization_steps
    print("  MC number of proposal is restricted to %d" % self.num_propose_limit)
    self.continueSample = True    
    print("Step: %5d [%10d]" % (self.loop,self.stat['number_of_proposed_steps']),end="\n")
    print("Current State: Mec(prev,post): %15.6f,%15.6f Img(prev,post): %15.6f,%15.6f" %
          (self.score_mech_prev,self.score_mech_curr,self.score_img_prev,self.score_img_curr))
    print("---------------------------------------------------------------------------------")
    self.nreject = 0
    self.check_continue()
    self.do_sample()  
    self.print_sample_statistics()  

  def finish(self,dumpFile=".dump"): 
    """
      Finish sampling by writing dumpFile
      dumpFile: Name of the dump file
    """   
    self.random_state = random.getstate(); self.random_state_np = np.random.get_state()
    with open(dumpFile, 'wb') as fp: 
      dill.dump(self, fp)

  def do_sample(self): 
    """
      This method dictates how sampling will be performed
    """   
    while self.continueSample:
      loop_count = self.loop+1; nprop = self.stat['number_of_proposed_steps'] + 1;
      print("Step: %5d [%10d]" % (loop_count,nprop),end="\n")
      self.set_beta()
      self.optimize_single_step() # core: performing one monte carlo step, accept or reject based on energy calculated from correlation coefficient between current and reference diffraction data
      print("\tScore1 : %15.6f  Score2: %15.6f pcc: %15.6f" % \
            (self.score_mech_curr,self.score_img_curr,self.pcc_val_curr))
      self.check_continue()#check steps (optimization steps)
      self.firstStep = False
      print("Current State: Mec(prev,post): %15.6f,%15.6f Img(prev,post): %15.6f,%15.6f" % 
            (self.score_mech_prev,self.score_mech_curr,self.score_img_prev,self.score_img_curr))
      print("---------------------------------------------------------------------------")

  def check_continue(self):
    """
      A check to see wheter to continue sample or not
    """
    nprop = self.stat['number_of_proposed_steps']
    cent = self.opt["OPTIMIZE"]["LIMIT_PROPOSAL_PERCENT"]    
    if self.opt['OPTIMIZE']['STOPBY'] == 'accept':
      if self.loop >= self.total_optimization_steps:
        print("Completed MC sampling\n")
        self.continueSample = False
    elif self.opt['OPTIMIZE']['STOPBY'] == 'propose':
      if nprop >= self.total_optimization_steps:
        print("Completed MC sampling\n")
        self.continueSample = False
      if nprop >=  self.num_propose_limit:
        print("Completed MC sampling, limited to maximum number of proposal\n")
        self.continueSample = False
      else:
        nacc = self.stat['number_of_accepted_steps']
        acceptance_cent = ((nacc*1.0) / nprop) * 100.0
        if acceptance_cent < cent and acceptance_cent > 0.:
          print("Completed MC sampling, limited to maximum proposal percent\n")
          self.continueSample = False
        else:
          print("Continue MC sampling, current rate: %.3f" % acceptance_cent)

  def set_beta(self):
    """
      Set beta against 3D scoring functions from anneal schedule or step count
    """
    if self.opt["OPTIMIZE"]["ANNEAL"]:
      if self.opt['OPTIMIZE']['STOPBY'] == 'accept':
        current_step_count = self.stat['number_of_accepted_steps']      
      elif self.opt['OPTIMIZE']['STOPBY'] == 'propose':
        current_step_count = self.stat['number_of_proposed_steps']
      for k in sorted(self.anneal_schedule.keys()):        
        if k <= current_step_count:
          beta = self.anneal_schedule[k]
      self.beta_cc = beta
      print("\tCurrent beta: %.6f" % beta)

  def read_schedule(self,filename):
    """
      To read anneal schedule file
      Filename: name of the anneal schedule file
      returns: schedule of annealling steps with corresponding temperature
    """
    if isfile(filename):
      anneal_schedule = {}
      with open(filename) as fp:        
        for line in fp:
          line = line.strip()
          first,second = line.split(" ")
          step = int(first)
          temp = float(second)
          anneal_schedule[step] = 1.0 / (temp  * 8.314)
      max_step = self.total_optimization_steps     
      if np.max(anneal_schedule.keys()) < max_step:
        anneal_schedule[max_step] = 1.0 / (self.opt["OPTIMIZE"]["TEMP"] * 8.314)
      return anneal_schedule
    else: return None

  def optimize_single_step(self):# core :performing one monte carlo step, accept or reject based on energy calculated from correlation coefficient between current and reference diffraction data
    """
      Helper method to optimize a single step. The method do_sample calls this
      method iteratively
    """
    verdict = False
    if not self.firstStep: 
      self.updateMovers(self.nreject)
    
    # take backup
    #ell = self.ellipsoid_opt[1]
    #for eindex in range(len(ell)):
    #  self.ellipsoid_backup[1][eindex] = self.ellipsoid_opt[1][eindex]		        

    self.propose(verbose=1) #no need to change propose
    if self.get_accept_mech(): # no need to change                  
      if self.get_accept_img():                        
        self.accept(); verdict = True;           
      else: 
        self.reject(roll_back_from=2)        
    else:
      self.reject(roll_back_from=1)    

    if verdict:
      self.stat['reject'] = self.nreject
      self.nreject = 0
      print("\tVerdict: Accept")
      self.stat['number_of_accepted_steps'] += 1
      if self.score_mech_prev >  self.score_mech_curr and self.score_img_prev >  self.score_img_curr:
        self.stat['number_of_down_steps'] += 1
      if self.score_mech_prev <= self.score_mech_curr and self.score_img_prev <= self.score_img_curr:
        self.stat['number_of_up_steps'] += 1
      if self.score_mech_prev >  self.score_mech_curr and self.score_img_prev <= self.score_img_curr:
        self.stat['number_of_neutral_steps'] += 1
      if self.score_mech_prev <=  self.score_mech_curr and self.score_img_prev >  self.score_img_curr:
        self.stat['number_of_neutral_steps'] += 1         
      self.cc[self.loop] = self.pcc_val_curr     
      self.outputBlock()      
      self.loop += 1;
    else:
      print("\tVerdict: Reject")
      self.nreject += 1

  def evaluate_restraints(self):
    """
      Evaluate 3D restrains (scoring functions)
      returns: Restraint object instances
    """
    Restr = copy.copy(self.restraints)
    for keyIndex in sorted(Restr.keys()):
      robj = Restr[keyIndex]
      if robj['EVAL']:
        inputkey = robj['INPUTKEY']
        if inputkey == 'ellipsoid':
          if robj['NAME'] == 'DOMAIN_PAIR' or robj['NAME'] == 'DOMAIN_EXCL':
            ellipsoids = self.ellipsoid_opt[1]
            robj['SCORE'] = robj['RESTR'].evaluate(ellipsoids)            
          else: robj['SCORE'] = 0.0
        else:
          robj['SCORE'] = 0.0
      else: robj['SCORE'] = 0.0
      Restr[keyIndex] = robj
    return Restr

  def get_evaluated_scores(self,Restr):
    """
      Calculates scores based on restrains objects
      Restr: Basically a dictionary (here a hash-of-hash) including
             its name, score value
    """
    scores = []
    for keyIndex in sorted(Restr.keys()):
      keyName = Restr[keyIndex]['NAME']
      sval = Restr[keyIndex]['SCORE']
      if Restr[keyIndex]['EVAL']:      
        scores.append(sval)
    return np.sum(scores)

  def print_evaluated_scores(self,Restr):
    """
      Prints scores based on restrains objects
      Restr: Basically a dictionary (here a hash-of-hash) including
             its name, score value
    """
    scores = []
    print("\t\tRESTRAINT TYPE               EVALUATED VALUE")
    print("\t\t--------------------------------------------")
    for keyIndex in sorted(Restr.keys()):
      keyName = Restr[keyIndex]['NAME']
      sval = Restr[keyIndex]['SCORE']
      if Restr[keyIndex]['EVAL']:
        print("\t\t%-28s %-15.4f" % (keyName,sval))
        scores.append(sval)
    print("\t\t%-28s %-15.4f" % ("Total Score:",np.sum(scores)))    
    return np.sum(scores)

  def get_accept_mech(self): 
    """
      Calculate whether to accept current move in terms of 3D scoring
      functions only.
      returns: True or False depending on acceptance or rejection
    """                   
    restr = self.evaluate_restraints()
    if self.opt['MODE'] == 'Debug':            
      self.score_mech_temp = self.print_evaluated_scores(restr)
    else: self.score_mech_temp = self.get_evaluated_scores(restr)    
    curr = self.score_mech_curr; next = self.score_mech_temp
    if curr == np.inf or next == np.inf:
      return False   
    return self.do_metropolis(curr,next,self.beta_mech)
   
  def get_accept_img(self):  # change only this and everything related 
    """
      Calculate whether to accept current move in terms of 2D image
      comparison
      returns: True or False depending on acceptance or rejection
    """  
    cc = self.get_correlation_from_diffractionImage(verbose=0)# change
    self.pcc_val_temp = cc
    method = self.opt['OPTIMIZE']['BIAS']['METHOD']
    param = self.opt['OPTIMIZE']['BIAS']['PARAM']    
    w = self.wt_factor
    self.score_img_temp = self.get_score_on_cc(cc,method=method,param=param,weight=w)# no need to change
    curr = self.score_img_curr; next = self.score_img_temp;
    return self.do_metropolis(curr,next,self.beta_cc) # no need to change

  def propose(self,method='random',verbose=1):
    """
      Propose a move
    """
    self.move_selector(method=method,verbose=verbose)
    if self.opt["OPTIMIZE"]["MOVERS"]["PROPOSALTYPE"] == 'SYNC':
      self.propose_as_rigid(verbose=verbose)
    else:
      for w in self.which:
        mv = self.movers[w]
        eindex = mv.index
        print("propose a move for %d[%s]" % (eindex,mv.name))
        if verbose > 0:
          c = self.ellipsoid_opt[1][eindex].mean
          print("Before move: center = %.6f %.6f %.6f" % (c[0],c[1],c[2]))
        self.ellipsoid_opt[1][eindex] = mv.propose()#important for understanding		
        if verbose > 0:
          c = self.ellipsoid_opt[1][eindex].mean
          print("After move: center = %.6f %.6f %.6f" % (c[0],c[1],c[2]))
      self.stat['number_of_proposed_steps'] += 1    
      self.remove_drift()

  def propose_as_rigid(self,verbose=1):
    """
      Propose a move identical to all ellipsoids
    """
    w = self.which[0]
    mv = self.movers[w]
    eindex = mv.index
    if mv.tval > 0.0: 
      tpropose = mv.propose_center(mv.tval)
    else:
      tpropose = np.zeros(3)
    if mv.rval > 0.0:
      rpropose = mv.propose_rotation(mv.rval)
    else:
      rpropose = np.array([1.0,0.0,0.0,0.0])
    if mv.sval > 0.0:
      spropose = mv.propose_scale(mv.sval)
    else:
      spropose = np.zeros(3)
    print("propose a move for %d[%s]" % (eindex,mv.name))    
    for w in self.which:
      mv = self.movers[w]
      eindex = mv.index
      print("propose a move for %d[%s]" % (eindex,mv.name))
      if verbose > 0:
        c = self.ellipsoid_opt[1][eindex].mean
        print("Before move: center = %.6f %.6f %.6f" % (c[0],c[1],c[2]))
      self.ellipsoid_opt[1][eindex] = mv.propose(tpropose=tpropose,rpropose=rpropose,spropose=spropose)
      if verbose > 0:
        c = self.ellipsoid_opt[1][eindex].mean
        print("After move: center = %.6f %.6f %.6f" % (c[0],c[1],c[2]))  
    self.stat['number_of_proposed_steps'] += 1    
    self.remove_drift()  

  def accept(self): #no need to change
    """ 
      Accept the proposed move
    """
    self.score_mech_prev = self.score_mech_curr
    self.score_mech_curr = self.score_mech_temp
    self.score_img_prev = self.score_img_curr
    self.score_img_curr = self.score_img_temp
    self.pcc_val_prev = self.pcc_val_curr
    self.pcc_val_curr = self.pcc_val_temp
    for w in self.which:
      mv = self.movers[w]
      eindex = mv.index
      self.ellipsoid_opt[1][eindex] = mv.accept()
      e = self.ellipsoid_opt[1][eindex]
      print("accepted radii: %f %f %f" % (e.scaled_radii[0],e.scaled_radii[1],e.scaled_radii[2]))
    
  def reject(self,roll_back_from,check=True):
    """
      Reject the proposed move
      roll_back_from: An interger to decide roll back from which algorithmic
                      step. 2 means roll back from 2nd Metropolis step, 1
                      means roll back from 1st Metropolis step. 1st Metropolis
                      step is based on 3D scoring function only, whereas the
                      2nd Metropolis step is based on 2D image similarity, 
                      which is evaluated only when 1st Metropolis step is 
                      accepted.
      check: Check ellipsoid state whether they comes back successfully
    """    
    for w in self.which:
      mv = self.movers[w]
      eindex = mv.index
      self.ellipsoid_opt[1][eindex] = mv.reject()    
    self.roll_back(roll_back_from=roll_back_from,check=check)

  def roll_back(self,roll_back_from,check=True):
    """
      Helper for reject method.
      roll_back_from: An interger to decide roll back from which algorithmic
                      step.
      check: Check ellipsoid state whether they comes back successfully
    """
    print("\trolling back from %d" % roll_back_from)       
    if check:
      flag1 = 0; flag2 = 0;
      restr = self.evaluate_restraints()    
      cscore = self.get_evaluated_scores(restr)
      if roll_back_from == 2:
        cc = self.get_correlation_from_diffractionImage(verbose=0)
        if np.abs(cc - self.pcc_val_curr) <= 1e-06:
          flag1 = 1
          print("\tChecking pcc  : %15.6f (pass)" % (cc),end="\n")
        else:
          print("\tMismatch : %15.6f vs %15.6f" % (cc,self.pcc_val_curr))
      else: flag1 = 1
      if cscore == np.inf and self.score_mech_curr == np.inf:
        flag2 = 1
        print("\tChecking score: %15.6f (pass)" % (cscore),end="\n")
      elif np.abs(cscore - self.score_mech_curr) <= 1e-06: 
        flag2 = 1
        print("\tChecking score: %15.6f (pass)" % (cscore),end="\n")
      if flag1 != 1 and flag2 != 1:
        raise Exception("Check Failed in rolling back, %d %d" % (flag1,flag2))

  def do_metropolis(self,E1,E2,beta,verbose=1):#no change
    """
      Perform Metropolis criterion
      E1,E2: Energy (score) values
      beta: 1/KT
      verbose: Verbosity
    """
    accept = True
    if verbose>0: 
      print("\tMC Metropolis: prev=%15.6f post=%15.6f " % (E1,E2),end="")
      if E2 > E1:
        bp = exp(-beta*(E2  - E1))
        rp = np.random.random_sample()      
        if bp < rp: accept = False
        print("BoltzFac:: %15.6f,rp=%15.6f,accept=%d" % (bp,rp,accept),end="\n")
      else:
        print("BoltzFac:: %15.6f,rp=%15.6f,accept=%d" % (1.0,0.0,accept),end="\n")
    else:
      if E2 > E1:
        bp = exp(-beta*(E2  - E1))
        rp = np.random.random_sample()      
        if bp < rp: accept = False
    return accept
  
  def reset_drift_method(self,coord,method='both'):
    """
      To reset drift removal method depending in meth requested
      coord: The coordinate of points from which drift will be removed
      method: Input method which will be checked
      returns: updated method name
    """
    m = "none"
    if method == 'both' or method == 'trans' or method == 'rota'              \
       or method == 'transx' or method == 'transy' or method == 'transz':
      m = method
    else:
      return m
    if coord.shape[0] < 3:
      if method == 'both':
        m = 'trans'
      if method == 'rota':
        m = 'none'
      if coord.shape[0] < 2:
        if method == 'trans':
          m = 'none'      
    return m

  def remove_drift(self):
    """
      Drift removal process
    """
    if self.drift_method == 'none':
      return            
    ecen = []
    ellipsoids = self.ellipsoid_opt[self.remove_drift_level][:]
    for e in ellipsoids:
      ecen.append(np.array(e.mean))
    print("\tTry to remove drift by %s\n" % self.drift_method)
    EllipCenCurr = np.asarray(ecen)          
    pc_fit = pcFit(target=self.EllipCenInit,source=EllipCenCurr,algo=self.drift_method)
    crdTransformed = pc_fit.fit()   
    for k in range(len(ellipsoids)):
      ellipsoids[k].mean = crdTransformed[k,:]
    self.ellipsoid_opt[self.remove_drift_level] = ellipsoids[:]

  def outputBlock(self):
    """
      Output various files as set in the options. Particularly, 
      OPTIMIZE->OUTPUT is used. OPTIMIZE->OUTPUT is itself a dictionary
      containing information about different outputs.
    """
    op = self.opt['OPTIMIZE']['OUTPUT']

    period = int(op['IMAGEPERIOD'])
    if not (self.loop + 1) % self.period:
      self.traj.add_data(data=self.ellipsoid_opt[1],datatype='Ellipsoid')

    if not (self.loop + 1) % period:           
      if op['IMAGE']:
        baseName = op['IMAGEPATH'] + "/" + op['IMAGEFILE']
        imgFile = baseName + str(self.loop + 1) + ".png"
        f, (ax1, ax2) = plt.subplots(1, 2)
        dw = float(self.opt["det_width"])  # change k: 2019.4.1
        ds = float(self.opt["det_dist"])
        wave = float(self.opt["wave_length"])
        k_inc = 1.0/wave
        k= (k_inc*dw/sqrt(ds*ds + 2*dw*dw))
        extent = [-k,k,-k,k]
        ax1.imshow(self.targetMat.T,origin = 'lower',extent = extent); #2019.9.2 Miyashita san chadded '.T'
        ax2.imshow(self.refMatClean.T,origin = 'lower',extent = extent);
	im1=ax1.imshow(self.targetMat.T,origin = 'lower',extent = extent);	
	im2=ax2.imshow(self.refMatClean.T,origin = 'lower',extent = extent);
	cbar1=plt.colorbar(im1,ax=ax1,pad=0.07,shrink=0.42) #2021.10.30 added colorbar
	cbar2=plt.colorbar(im2,ax=ax2,pad=0.07,shrink=0.42) #2021.10.30 added colorbar
	#cbar1.set_ticks([-4,-3,-2,-1,0]) #2021.11.1 added for EF2_7
	#cbar2.set_ticks([-4,-3,-2,-1,0]) #2021.11.1 added for EF2_7
	#im1.set_clim(-4,0) # added for EF2_7
	#im2.set_clim(-4,0) # added for EF2_7
	#cbar1.set_ticks([-3.0,-2.5,-2,-1.5,-1,-0.5,0]) #2021.11.1 added for AK
	#cbar2.set_ticks([-3.0,-2.5,-2,-1.5,-1,-0.5,0]) #2021.11.1 added for AK
	#im1.set_clim(-3,0) # added for AK
	#im2.set_clim(-3,0) # added for AK
	ax1.xaxis.set_major_locator(plt.MaxNLocator(7)) #2021.11.1 changing number of ticks
        ax2.xaxis.set_major_locator(plt.MaxNLocator(7)) #2021.11.1 changing number of ticks
	ax1.tick_params(axis='both',which='major',labelsize=8,direction='out') #2021.10.31 changed tick size
	ax2.tick_params(axis='both',which='major',labelsize=8,direction='out') #2021.10.31 changed tick size
	plt.subplots_adjust(wspace=0.35) #2021.10.30 ASI HAN changed the  subplots space
	#f.tight_layout 2021.10.30 this function adjust the subplots to fit into the figure area 		perfectly
        f.savefig(imgFile)
        plt.close()        

    period = int(op['GMDATPERIOD'])
    if not (self.loop + 1) % period:
      if op['GMDAT']:
        baseName = op['GMDATPATH'] + "/" + op['GMDATFILE']
        gmDatFile = baseName + str(self.loop + 1) + ".gmm"
        gmData = self.gmData
        gmData.write(gmDatFile)

    period = int(op['STATPERIOD'])
    if period != 1: 
       print("\t!!!WARNING!!! It is not recommended to skip writting in stat files")
    if not (self.loop + 1) % period:
      if op['STAT']:
        baseName = op['STATPATH'] + "/" + op['STATFILE']
        statFile = baseName + "_stat.out"
        with open(statFile,'a') as fp:
          fp.write("%d %.6f %d\n" % (self.loop+1,self.pcc_val_curr,self.stat['reject']+1))
    else:
      print("\tSkipping writting in stat: %d %.6f %d" % (self.loop+1,self.pcc_val_curr,
                                                         self.stat['reject']+1))

  def readStatFile(self):
    """
      To read OPTIMIZE->OUTPUT->{STATPATH/STATFILE}__stat.out
      returns: The read data
    """
    op = self.opt['OPTIMIZE']['OUTPUT']
    baseName = op['STATPATH'] + "/" + op['STATFILE']
    statFile = baseName + "_stat.out"
    StatBlock = {'trj': [], 'summary': []}
    with open(statFile,'r') as fp:
      lines = fp.readlines()
      for line in lines:
        line = line.strip()
        if line[0] == '#':  StatBlock['summary'].append(line)
        else: StatBlock['trj'].append(line)
    return StatBlock

  def writeStatBlock(self,StatBlock,summary=False):
    """
      To write OPTIMIZE->OUTPUT->{STATPATH/STATFILE}__stat.out
      StatBlock: The data present in the above file
      summary: Whether to write the summary data at the end
    """
    op = self.opt['OPTIMIZE']['OUTPUT']
    if op['STAT']:
      baseName = op['STATPATH'] + "/" + op['STATFILE']
      statFile = baseName + "_stat.out"
      trjStat = StatBlock['trj']; summStat = StatBlock['summary']
      with open(statFile,'w') as fp:
        for i in range(len(trjStat)):
          fp.write("%s\n" % trjStat[i])
        if summary:
          for i in range(len(summStat)):
            fp.write("%s\n" % summStat[i])

  def get_correlation_from_diffractionImage(self,verbose=1): #change / rewrite
    """
      verbose: verbosity
    """
    if self.timing: start_clock = datetime.now() 

    det_width = float(self.opt['det_width'])
    det_dist  = float(self.opt['det_dist']) 
    wave_length = float(self.opt['wave_length']) 
          
    ell = self.ellipsoid_opt[1]
    eweight = [None] * len(ell)
    for i in range(len(ell)):
      eweight[i] = ell[i].mass
    sumw = np.sum(eweight)       
    for i in range(len(ell)):
      eweight[i] = eweight[i]/sumw
    ecen = []
    ecov = []
    for i in range(len(ell)):      
      ecen.append(ell[i].mean)
      ecov.append(ell[i].covar)
    gmnew = GMData(weight=eweight,center=ecen,rotation=ecov,makeCov=False) #1
    gmnew.init_ctypes_of_gmm()
    tmat = gmnew.call_get_abs_struct_factor_all_mesh_on_Ewald_sphere(200,det_width,det_dist,wave_length)
    targetMat= np.log10(tmat)   
    
    targetMat = targetMat.astype(np.float64)
    if self.opt["correlation"] == "pearson": 
      A = targetMat.ravel()
      B = self.refMatClean.ravel()
      cc,pval = pearsonr(A,B)  #2
    elif self.opt["correlation"] == "nagai":
      dw = float(self.opt["det_width"])  # change k: 2019.4.1
      ds = float(self.opt["det_dist"])
      wave = float(self.opt["wave_length"])
      k_inc = 1.0/wave
      k= (k_inc*dw/sqrt(ds*ds + 2*dw*dw))
      cd = float(self.opt["circle_width"])
      lowest_k = float(self.opt["lowest_k"])
      ntheta=360
      rlist = np.arange(lowest_k,k+cd,cd)
      dtheta=2*np.pi/ntheta
      theta_list=np.arange(0.0,2*np.pi,dtheta)
      coslist,sinlist=np.cos(theta_list),np.sin(theta_list)
      f1list=np.zeros((len(rlist),len(theta_list)))
      for i_r in range(len(rlist)):
        for each in range(len(theta_list)):
          kx,ky,kz=TOOLS.get_coord_on_ewald_sphere_from_zero_curvature_coord(rlist[i_r]*coslist[each],rlist[i_r]*sinlist[each],k_inc)
          f1list[i_r][each]=gmnew.call_get_abs_struct_factor(kx,ky,kz)
      f1list=np.log10(f1list)*2.0
      cc = TOOLS.get_overlap_with_jit(f1list,self.f2list)

    print("diffraction gives correlation: %.6f" % cc)
    #col = self.map_correlation_to_color(cc)    
    col = np.array([0.0,1.0,0.0])
    self.targetMat = targetMat
    self.gmData = gmnew       
    if verbose > 0: print("\tpcc   : %15.6f" % (cc),end="\n")    
    return cc

  def reset_loop(self):
    """
      Reset loop parameters for optimization and reset other
      summary statistics
    """
    self.loop = 0    
    self.cc = np.zeros(self.total_optimization_steps)
    self.stat = {}
    self.stat['number_of_accepted_steps'] = 0
    self.stat['number_of_down_steps'] = 0
    self.stat['number_of_up_steps'] = 0
    self.stat['number_of_proposed_steps'] = 0
    self.stat['number_of_neutral_steps'] = 0
    self.stat['reject'] = 0
      
  def print_sample_statistics(self):
    """
      Append summary about sampling at the end of STATFILE
    """
    self.cc = self.cc[0:self.loop]
    cc = np.append(self.cc0,self.cc)    
    cc_min = np.min(cc)
    cc_max = np.max(cc)
    cc_mean = np.mean(cc)
    op = self.opt['OPTIMIZE']['OUTPUT']
    if op['STAT']:
      baseName = op['STATPATH'] + "/" + op['STATFILE']
      statFile = baseName + "_stat.out"
      fp = open(statFile,'a')
      fp.write("# proposed steps %d\n" % (self.stat['number_of_proposed_steps']))
      fp.write("# accepted steps %d\n" % (self.stat['number_of_accepted_steps']))
      fp.write("# downward steps %d\n" % (self.stat['number_of_down_steps']))
      fp.write("#  upward steps %d\n" % (self.stat['number_of_up_steps']))
      fp.write("# neutral steps %d\n" % (self.stat['number_of_neutral_steps']))
      fp.write("# Correlation: Min,Max,Mean = %.4f %.4f %.4f\n" % (cc_min,cc_max,cc_mean))            
      fp.close()
    else:
      print("proposed steps %d" % (self.stat['number_of_proposed_steps']))
      print("accepted steps %d" % (self.stat['number_of_accepted_steps']))
      print("downward steps %d" % (self.stat['number_of_down_steps']))
      print("  upward steps %d" % (self.stat['number_of_up_steps'])) 
      print("# neutral steps %d\n" % (self.stat['number_of_neutral_steps']))   
      print("Correlation: Min,Max,Mean = %.4f %.4f %.4f" % (cc_min,cc_max,cc_mean))
 
  def get_score_on_cc(self,cc,method='naive',param=1.0,weight=1.0):
    """
      Calculate score based on 2D similarity
      cc: 2D similarity
      method: One of 'naive', 'harmonic' and 'linear'
      param,weight: Score calculation parameters
      returns: The value of the score
    """
    if cc == -2.0:
      return 0.0
    if method == 'naive':
      rv = -cc * weight
    elif method == 'harmonic':
      delta = cc - param
      rv = delta * delta * weight
    elif method == 'linear':
      delta = cc - param
      rv = delta * weight
    return rv

  def map_correlation_to_color(self,cc):
    """
      Convert correlation (2D) to a color
    """
    if len(self.cc) > 0:
      cc0 = self.cc[0]
    else: cc0 = cc
    colrange = np.linspace(0.,1.,100)
    crange = np.linspace(cc0,1.0,100)
    cid = np.digitize(cc,crange)    
    if cc < cc0:      
      red  = colrange[cid]
      green= 1.0 - colrange[cid]
      blue = 0.0
    elif cc > cc0:      
      red  = 0.0
      green= 1.0 - colrange[cid]
      blue = colrange[cid]
    else:      
      red  = 0.0
      green= 1.0
      blue = 0.0
    return np.array([red,green,blue])
