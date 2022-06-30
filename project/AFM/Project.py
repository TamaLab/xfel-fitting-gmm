from AFM_Header import *
from .AFM import *

class Project:
  """
    This is like a soup. This class uses different AFM classes to perform
    specified project. Two important projects are 1) sampling 2) emulation
    Executables are written (afmSampler and afmEmulator) which are wrapper
    over this class. Currently class 'Project' has methods that are build
    on AFM classes, in order to simplify calls to AFM class.

    Currently Project (as well as AFM class) supports two types of restraints
    1) Domain-to-domain (center-to-center) pairwise interaction restrained
       to the native distances
    2) Excluded volume score that depends on pairwise correlation between 
       ellipsoids that define domains. Currently excluded volume score can be
       infinity or zero, depending on a cutoff on correlation, if correlation
       is higher than the cutoff score is infinity.

  """
  RESTR_TYPES = {};
  RESTR_NAME = [{'NAME': 'DOMAIN_PAIR'},
                {'NAME': 'DOMAIN_EXCL'}]
  def __init__(self,opt):    
    """
      Set options depending on input option dictionary
      opt: Input option
    """
    self.opt = opt.opt
    self.RESTR_TYPES = Project.RESTR_TYPES
    self.RESTR_NAME = Project.RESTR_NAME
    for i in range(len(self.RESTR_NAME)): self.RESTR_TYPES[i] = self.RESTR_NAME[i]
    self.pdb = None

  def set_input(self,asis=True,set_representation=True):
    """
      Set input object depending on its format
      Depending project type (EMULATE or SAMPLE) a transformation to input 
      data is made. However, asis option suppresses any transformation.
      For project type SAMPLE, if asis is False a transformation is guessed
      (a recommended workflow, see afmSampler executable)
    """
    if self.opt["PROJECT"] == "EMULATE":
      buffx = float(self.opt["MESH"]["BUFFER"]["X"])
      buffy = float(self.opt["MESH"]["BUFFER"]["Y"])
      buffz = float(self.opt["MESH"]["BUFFER"]["Z"])
      transf = Transformation(rot=np.diag([1.0,1.0,1.0]),
                              trans=np.array([buffx,buffy,buffz]),rotationType='m')
    elif self.opt["PROJECT"] == "SAMPLE":      
      transf = None
    input_fmt = self.opt['INPUTFILE']['FORMAT']
    input_objfile = self.opt['INPUTFILE']['NAME']
    if input_fmt == 'pdb':
      self.set_pdb(input_objfile,asis=asis,transf=transf)
      if set_representation:
        self.representation = Representation(self.pdb,self.opt)
    elif input_fmt == 'gmm':
      self.pdb = None
      if set_representation:
        self.representation = Representation(self.pdb,self.opt)        
        self.representation.apply_transformation(asis,transf)

  def set_pdb(self,filename,asis=True,transf=None):
    """
      Sets pdb object as read from option.
      Depending project type (EMULATE or SAMPLE) a transformation to input 
      coordinate is made. However, asis option suppresses any transformation.
      For project type SAMPLE, if asis is False a transformation is guessed
      (a recommended workflow, see afmSampler executable)
      filename: Name of the pdb file
      asis: Whether to transform or not
      transf: transformation details (shift and rotation)
    """
    self.pdb = PDB(filename,asis=asis,transf=transf)

  def write_pdb(self):
    """
      Writes pdb object to a file named BaseName_trs.pdb, where BaseName is
      derived from input pdb filename
    """
    pdbInputRoot = os.path.basename(self.opt['INPUTFILE']['NAME'])    
    self.pdb.write_pdb(self.pdb.object,pdbInputRoot.split(".")[0] + "_trs.pdb")

  def set_representation(self):
    """
      Make ellipsoidal representation from PDB data
    """
    self.representation = Representation(self.pdb,self.opt)

  def set_mesh(self):
    """
      Make mesh from ellipsoidal representation
    """
    self.representation.SetupMesh()
    self.opt['dpi'] = self.representation.opt['dpi']

  def set_image(self,expandProperties=False):
    """
      Read input AFM-like image and optionally calculates different settings
      required later
      expandProperties: Whether to expand image properties for further 
                        processing of coordinate data
      Note: set expandProperties = True when reading in sampling work-flow
            see afmSampler executable
    """
    imgName = self.opt["IMGFILE"]
    self.img = AFMInputImage(imgName,fromarray=False,normalize=False)
    if expandProperties:
      self.set_resolutions()
      self.set_image_property()

  def set_resolutions(self):
    """
      Helper method for set_image to expand properties. This sets
      x,y,z resolutions. The user is expected to know x,y-resolution
      (as set in MESH->STEP->{X or Y}). For z resolution the user
      is expected to know dynamic range of z-values in Angstrom, but 
      this can be approximate.
    """
    self.xres = self.get_resolution('X')
    self.yres = self.get_resolution('Y')
    self.xlen = self.get_LengthFromStep(self.xres,'X')
    self.ylen = self.get_LengthFromStep(self.yres,'Y')
    if self.opt['MESH']['STEP']['Z']['FORMAT'] == 'LENGTH':
      zlen = self.opt['MESH']['STEP']['Z']['VALUE']
      self.zres = self.get_StepFromLength(zlen,'Z')
    elif self.opt['MESH']['STEP']['Z']['FORMAT'] == 'STEP': 
      self.zres = self.opt['MESH']['STEP']['Z']['VALUE']
      zlen = self.get_LengthFromStep(self.zres,'Z')
    self.zlen = zlen
    self.ncolor = zlen

  def set_image_property(self):
    """
      Helper method for set_image to expand properties. This sets
      scale along xyz, which will be used later as grid-points
      When project type is SAMPLE is sets up scales, see documentation
      for setup_scales
    """
    if self.opt["PROJECT"] == "SAMPLE":
      self.img_height = self.img.imatClean.shape[0] # height = number of row
      self.img_width  = self.img.imatClean.shape[1] # width  = number of col
      self.imgObjectOffsetHW = {0: {"pixel": None, "A": None}, 1: {"pixel": None, "A": None}}
      nxpixel = np.min(np.where(self.img.imatClean != 255.0)[0])
      nypixel = np.min(np.where(self.img.imatClean != 255.0)[1])
      self.darkest = np.min(self.img.imatClean)
      self.lightest = np.max(self.img.imatClean[self.img.imatClean != 255.0])
      self.imgObjectOffsetHW[0]["pixel"] = nxpixel
      self.imgObjectOffsetHW[1]["pixel"] = nypixel
      self.imgObjectOffsetHW[0]["A"] = self.pixel2A(nxpixel,self.xres)
      self.imgObjectOffsetHW[1]["A"] = self.pixel2A(nypixel,self.yres)
      self.opt['MESH']['XSCALE'] = {'FROM': 0.0, 'TO': (self.img_height-1) * self.xres, 'BY': self.xres}
      self.opt['MESH']['YSCALE'] = {'FROM': 0.0, 'TO': (self.img_width-1)  * self.yres, 'BY': self.yres}
      self.opt['MESH']['ZSCALE'] = {'FROM': 0.0, 'LENGTH': self.ncolor, 'BY': self.zres}
      self.setup_scales()
    elif self.opt["PROJECT"] == "EMULATE":
      self.imgObjectOffsetHW = {0: {"pixel": None, "A": None}, 1: {"pixel": None, "A": None}}
      if self.opt['INPUTFILE']['FORMAT'] == 'pdb':
        xmax = self.pdb.objMax[0] + float(self.opt["MESH"]["BUFFER"]["X"])
        ymax = self.pdb.objMax[1] + float(self.opt["MESH"]["BUFFER"]["Y"])
        zmax = self.pdb.objMax[2] + float(self.opt["MESH"]["BUFFER"]["Z"])
      elif self.opt['INPUTFILE']['FORMAT'] == 'gmm':
        xmax = self.representation.repMax[0] + float(self.opt["MESH"]["BUFFER"]["X"])
        ymax = self.representation.repMax[1] + float(self.opt["MESH"]["BUFFER"]["Y"])
        zmax = self.representation.repMax[2] + float(self.opt["MESH"]["BUFFER"]["Z"])
      xlim = self.opt["MESH"]["RANGE"]["X"]
      ylim = self.opt["MESH"]["RANGE"]["Y"]
      zlim = self.opt["MESH"]["RANGE"]["Z"]
      if xmax <  xlim and ymax < ylim and zmax < zlim:
        self.opt['MESH']['RANGE']['X'] = xmax
        self.opt['MESH']['RANGE']['Y'] = ymax
        self.opt['MESH']['RANGE']['Z'] = zmax
        self.set_resolutions()
        self.opt['MESH']['XSCALE'] = {'FROM': 0.0, 'TO': xmax, 'BY': self.xres}
        self.opt['MESH']['YSCALE'] = {'FROM': 0.0, 'TO': ymax, 'BY': self.yres}
        self.opt['MESH']['ZSCALE'] = {'FROM': 0.0, 'TO': zmax, 'BY': self.zres}        
        self.setup_scales()
        print("Range along X(updated): %8.3f" % self.opt['MESH']['XSCALE'][1])
        print("Range along Y(updated): %8.3f" % self.opt['MESH']['YSCALE'][1])
        print("Range along Z(updated): %8.3f" % self.opt['MESH']['ZSCALE'][1])
        self.opt['MESH']['OFFSET'] = [0.,0.,0.]
      else:
        raise Exception("axis ranges too short")

  def check_zcolor(self,ncolor,zres,dark,light):
    """
      Check color mapping using z-coordinate
      ncolor: Number of color divisions along z (usually 256)
      zres: Resolution along z-direction
      dark: Darkest pixel in image
      light: Lightest pixel in image
      return: Check result
    """
    err = 0
    if light > ncolor: err += 1
    pdb_present = True
    try:
      self.pdb
    except NameError:
      pdb_present = False
    if pdb_present:
      zmax = self.pdb.objMax[2]
      ztrans = (ncolor - light) * zres
      colorMax = 255
      darkest_possible = colorMax - int((zmax + ztrans) / zres)
      if darkest_possible < 0: err += 1
    if err > 0: return False
    else:       return True

  def get_resolution(self,whichAxis):
    """
      This calculates resolution (typically for xy direction)
      using input option already saved
    """
    value = self.opt['MESH']['STEP'][whichAxis]['VALUE']
    format = self.opt['MESH']['STEP'][whichAxis]['FORMAT']
    if format == 'STEP':
      return value
    elif format == 'LENGTH':
      return self.get_StepFromLength(value,whichAxis)

  def get_StepFromLength(self,value,whichAxis):
    """
      From dynamic range along a direction axis calculate a resolution
      in Angstrom using total length along the direction in Angstrom
      returns: resolution along a axis
    """
    vrange = self.opt['MESH']['RANGE'][whichAxis]
    if vrange != 'unk':
      return vrange / (value - 1.0)
    else:
      raise Exception("can't set stepsize along axis %s" % whichAxis)  
    
  def get_LengthFromStep(self,value,whichAxis):
    """
      From dynamic range along a direction axis calculate a total
      length along the direction using resolution in Angstrom
    """
    vrange = self.opt['MESH']['RANGE'][whichAxis]
    if vrange != 'unk':
      return int(ceil((vrange/value)+1));
    else:
      raise Exception("can't set length along axis %s" % whichAxis) 

  def pixel2A(self,npixel,resolution,beg=0):
    """
      Convert number of pixel and resolution information to a length
      in Angstrom
    """
    return beg + resolution*npixel

  def create_sequence(self,sdict):
    """
      This is an utility method to convert a dictionary with keys 'FROM',
      'TO' and 'BY' or 'LENGTH' to a sequence of numbers
      sdict: Input dictionary with aforesaid keywords
      returns: the sequence as described above
    """
    algo = 0
    if sdict.has_key("FROM"): 
      try:
        Beg = float(sdict["FROM"])
      except ValueError:
        sdict.pop("FROM",None)
    if sdict.has_key("TO"):   
      try:
        End = float(sdict["TO"])
      except ValueError:
        sdict.pop("TO",None)
    if sdict.has_key("BY"):
      try:
        By = float(sdict["BY"])
      except ValueError:
        sdict.pop("BY",None)
    if sdict.has_key("LENGTH"):   
      try:
        Len = int(sdict["LENGTH"])
      except ValueError:
        sdict.pop("LENGTH",None)    
    if sdict.has_key("FROM") and sdict.has_key("TO"):
      if sdict.has_key("BY"):
        if End >= Beg: algo = 1          
        else: algo = 2        
      if sdict.has_key("LENGTH") and not sdict.has_key("BY"):
        if End >= Beg: algo = 3
        else: algo = 4
    if sdict.has_key("FROM") and not sdict.has_key("TO"):
      if sdict.has_key("BY") and sdict.has_key("LENGTH"): algo = 5
    if sdict.has_key("TO") and not sdict.has_key("FROM"):
      if sdict.has_key("BY") and sdict.has_key("LENGTH"): algo = 6    
    if not algo: return None
    else:
      if algo == 1:        
        nelm = int(ceil(((End - Beg)/By)+1)); step = By; elm = [None] * nelm;
        for i in range(nelm): elm[i] = Beg + step*i
      elif algo == 2:
        nelm = int(ceil(((Beg - End)/By)+1)); step = By; elm = [None] * nelm;
        for i in range(nelm): elm[i] = End - step*i
      elif algo == 3:
        nelm = Len; step = (End - Beg)/(nelm - 1.0); elm = [None] * nelm;
        for i in range(nelm): elm[i] = Beg + step*i
      elif algo == 4:
        nelm = Len; step = (Beg - End)/(nelm - 1.0); elm = [None] * nelm;
        for i in range(nelm): elm[i] = End - step*i
      elif algo == 5:
        nelm = Len; step = By; elm = [None] * nelm;        
        for i in range(nelm): elm[i] = Beg + step*i
      elif algo == 6:
        nelm = Len; step = By; elm = [None] * nelm;
        for i in range(nelm): elm[i] = End - step*i
      return (elm,step)

  def setup_scales(self):
    """
      Using the range of values along x,y,z the grid-points are
      calculated here
    """
    opt = self.opt
    if opt.has_key("MESH"):
      if opt["MESH"].has_key("XSCALE"):
        g,step = self.create_sequence(opt["MESH"]["XSCALE"])
        opt["MESH"]["XGRID"] = g
        opt["MESH"]["XSCALE"] = [g[0],g[-1]]
        opt["MESH"]["XRESOLUTION"] = step
      else:
        print("Error: can't set x-scales")
      if opt["MESH"].has_key("YSCALE"):
        g,step = self.create_sequence(opt["MESH"]["YSCALE"])
        opt["MESH"]["YGRID"] = g
        opt["MESH"]["YSCALE"] = [g[0],g[-1]]
        opt["MESH"]["YRESOLUTION"] = step
      else:
        print("Error: can't set y-scales")
      if opt["MESH"].has_key("ZSCALE"):
        g,step = self.create_sequence(opt["MESH"]["ZSCALE"])
        opt["MESH"]["ZGRID"] = g
        opt["MESH"]["ZSCALE"] = [g[0],g[-1]]
        opt["MESH"]["ZRESOLUTION"] = step
      else:
        print("Error: can't set z-scales")
    else:
      print("Error: can't set scales")
    self.opt = opt

  def get_initialCorrelation(self,shift=True,plotFile=None): 
    """
      The method that used IMGCOMPARE options to calculate image similarity
      shift: Whether to shift the projected matrix
      plotFile: Output filename for resultant image
    """   
    self.set_image_compare_options(shift=shift)
    to = self.targetOpt
    ro = self.referenceOpt
    rmat = self.referenceMat
    tmat = self.targetMat
    superpose = self.opt["IMGCOMPARE"]["SUPERPOSE"]
    compare = ImageCompare(target=tmat,reference=rmat,
                           targetProcessOption=to,referenceProcessOption=ro,
                           postprocess=True,superpose=superpose)
    corrOpt1 = self.opt["IMGCOMPARE"]["CORRELATION"]["OVERLAP"]
    corrOpt2 = self.opt["IMGCOMPARE"]["CORRELATION"]["UNION"]
    qval = self.opt["IMGCOMPARE"]["CORRELATION"]['Q']
    self.initCorr = compare.getCorrelation(overlapOnly=corrOpt1,unionOnly=corrOpt2,q=qval)
    print("Correlation with OVERLAP(%r) and UNION(%r) = %.3f" % (corrOpt1,corrOpt2,self.initCorr))
    if not plotFile is None:
      corrOpt1 = self.opt["IMGCOMPARE"]["PLOTCORR"]["OVERLAP"]
      corrOpt2 = self.opt["IMGCOMPARE"]["PLOTCORR"]["UNION"]
      qval = self.opt["IMGCOMPARE"]["PLOTCORR"]['Q']
      compare.plotCompareImage(plotFile,overlapOnly=corrOpt1,unionOnly=corrOpt2,q=qval)
      print(" Corresponding compared image is in %s" % plotFile)    

  def set_image_compare_options(self,shift=True):
    """
      Helper method for get_initialCorrelation that optionally shift projected
      image using 'ShiftProjectedImage' method
      shift: Whether to shift the projected image
    """
    crop = self.opt["IMGCOMPARE"]["CROP"]["PROJECTION"]
    normalize = False
    mode = self.opt["IMGCOMPARE"]["MODE"]
    bgcorr = self.opt["IMGCOMPARE"]["BGCORRECTION"]["PROJECTION"]
    bgfill = self.opt["IMGCOMPARE"]["BGFILL"]
    method = self.opt["IMGCOMPARE"]["CORRELATION_MEASURE"]
    self.targetOpt = ImageProcessOption(fromarray=True,crop=crop,bgcorr=bgcorr,
                                        mode=mode,normalize=False,bgfill=bgfill,
                                        measure=method)
    crop = self.opt["IMGCOMPARE"]["CROP"]["IMAGE"]
    bgcorr = self.opt["IMGCOMPARE"]["BGCORRECTION"]["IMAGE"]
    self.referenceOpt = ImageProcessOption(fromarray=True,crop=crop,bgcorr=bgcorr,
                                           mode=mode,normalize=False,bgfill=bgfill,
                                           measure=method)
    self.referenceMat = self.img.imatClean
    if shift:
      self.targetMat = self.ShiftProjectedImage(self.representation.colorMat)
    else:
      self.targetMat = copy.deepcopy(self.representation.colorMat).astype(np.float64)

  def ShiftProjectedImage(self,pmat):
    """
      To shift projected image from coordinate data so as to match to the
      reference image. The reference image containing a molecule is normaly
      recorded as shifted towards the middle of the image. However, projected
      image as obtained from coordinate data contain a molecule whose topleft
      corner is usually at the topleft corner of the image (i.e. no offset).
      Here such offset will be defined to be used in 'MoveStructureByImage'
      method.
      pmat: Projected image matrix
      returns: Shifted image matrix
    """
    rmat = self.img.imatClean
    pmat = copy.deepcopy(pmat).astype(np.float64)
    rb = rmat != 255.; pb = pmat != 255.;    
    whereX,whereY = np.where(rmat != 255.0)
    nxpixel = np.min(whereX)
    nypixel = np.min(whereY)
    nxA = self.pixel2A(nxpixel,self.xres)
    nyA = self.pixel2A(nypixel,self.yres)
    zmaxRef,zminRef = np.max(rmat[rb]),np.min(rmat[rb])
    zmax,zmin = np.max(pmat[pb]),np.min(pmat[pb])
    zmove = zmin - zminRef
    pmat[pb] = pmat[pb] - zmove
    nzA = self.pixel2A(zmove,self.zres)
    omat = np.full(pmat.shape,255.0)
    for i in range(rmat.shape[0]-nxpixel):
      for j in range(rmat.shape[1]-nypixel):
        omat[i+nxpixel,j+nypixel] = pmat[i,j]
    self.opt['MESH']['OFFSET'] = [nxA,nyA,nzA]
    return omat

  def MoveStructureByImage(self):
    """
      In method 'ShiftProjectedImage' offsets are defined. Here those offsets
      are used to move pdb-coordinate data. Later, such transformed pdb data
      will be used.
    """

    rot=np.diag([1.0,1.0,1.0])
    trans = np.array(self.opt['MESH']['OFFSET'])
    if self.opt['INPUTFILE']['FORMAT'] == 'pdb':
      coord = self.pdb.get_coordinates(self.pdb.object)
      coordNew = self.pdb.transform_coord(coord,trans,rot)
      self.pdb.object,coord = self.pdb.set_coordinates(self.pdb.object,coordNew)
    elif self.opt['INPUTFILE']['FORMAT'] == 'gmm':
      transf = Transformation(rot=rot,trans=trans,rotationType='m')
      self.representation.apply_transformation(asis=False,transf=transf)

  def set_movers(self):
    """
      Set mover list according to the options defined in OPTIMIZE->MOVERS
    """
    movers_opt = self.opt["OPTIMIZE"]["MOVERS"]
    mover_level = 0
    dom_exists = True
    try:
      self.representation.DomEllipsoids
    except:
      dom_exists = False
    res_exists = True
    try:
      self.representation.ResEllipsoids
    except:
      res_exists = False
    if dom_exists:
      ndom = len(self.representation.DomEllipsoids)
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
    if res_exists:
      nres = len(self.representation.ResEllipsoids)
      if movers_opt['RESIDUE'] == 'All':
        movers_opt['RESIDUE'] = range(1,nres+1)
        mover_level += 1
      elif movers_opt['RESIDUE'] == 'None':
        movers_opt['RESIDUE'] = []  
      else:
        mlist = []
        for m in movers_opt['RESIDUE'].split(","):
          try:
            im = int(m)
          except ValueError:
            raise Exception("Mover List Residue includes non-integer")
          if im > nres or im < 1:
            raise Exception("Mover List Residue improper")
          mlist.append(im)
        movers_opt['RESIDUE'] = mlist
        mover_level += 1
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
    self.opt["OPTIMIZE"]["MOVERS"] = movers_opt
    if self.opt["MESH"]["LEVEL"] == "domain":
      self.opt["MESH"]["LEVEL"] = 1
    elif self.opt["MESH"]["LEVEL"] == "residue":
      self.opt["MESH"]["LEVEL"] = 0
    else:
      raise Exception("Mesh level setting wrong")

  def set_restraints(self):
    """
      Set restraint calculations (scores based on 3D coordinate data)
      Particularly important settings are under the RESTRAINTS key
    """
    closePairList = []
    RestraintDict = self.RESTR_TYPES

    keyName = 'DOMAIN_PAIR'
    keyIndex = self.get_key_index(keyName,RestraintDict)
    do_eval = self.opt['RESTRAINTS'][keyName]['EVALUATE']
    if do_eval:
      try:
        pairs_read = self.opt['DOMAIN'][keyName]
        pairs = []
        for p in pairs_read: 
          ps = p.split(","); 
          ps[0] = ps[0].strip(); ps[1] = ps[1].strip(); 
          pairs.append(ps)
      except KeyError:
        self.generate_domain_pairs()
        pairs = self.opt['DOMAIN'][keyName]
      strength = self.opt['RESTRAINTS'][keyName]['STRENGTH']
      rs = EllipsoidPairDistRestraint(self.representation.DomEllipsoids,pairs,
                                 strength=strength)
    else: rs = ConstantRestraint(value=0.0)
    RestraintDict[keyIndex]['RESTR'] = rs
    RestraintDict[keyIndex]['EVAL'] = do_eval
    RestraintDict[keyIndex]['INPUTKEY'] = "ellipsoid"
    RestraintDict[keyIndex]['VALUE'] = None

    keyName = 'DOMAIN_EXCL'
    keyIndex = self.get_key_index(keyName,RestraintDict)
    do_eval = self.opt['RESTRAINTS'][keyName]['EVALUATE']
    if do_eval: 
      cutoff = self.opt['RESTRAINTS'][keyName]['STRENGTH']      
      rs = EllipsoidExcludedVolumeRestraint(self.representation.DomEllipsoids,
                                       cutoff=cutoff)
    else: rs = ConstantRestraint(value=0.0)
    RestraintDict[keyIndex]['RESTR'] = rs
    RestraintDict[keyIndex]['EVAL'] = do_eval
    RestraintDict[keyIndex]['INPUTKEY'] = "ellipsoid"
    RestraintDict[keyIndex]['SCORE'] = None
      
    self.Restr = RestraintDict

  def get_key_index(self,keyName,RestraintDict):
    """
      Helper method for setting up restraint calculation which finds
      key index from internal RestraintDict given string present in
      user-given yaml input.
      keyName: Key string present in the dictionary RestraintDict
      RestraintDict: An input dictionary including different values
                     related to restraint calculation, e.g.
                     Child-key: Key index indicating what type of
                                restraint function. Currently there
                                are two options - 'DOMAIN_PAIR' for 
                                domain-to-domain distance restraint
                                score indexed by integer 0, 'DOMAIN_EXCL'
                                for exlcuded volume based score indexed
                                by integer 1.
                                This child-key is the parent of the 
                                next level keys. The next-level keys
                                are the following:
                     NAME: String to include the name of restraint type,
                           either of 'DOMAIN_PAIR' or 'DOMAIN_EXCL' now
                     RESTR: An object of instance EllipsoidPairDistRestraint
                            or EllipsoidExcludedVolumeRestraint
                     EVAL: Whether to evaluate this score
                     INPUTKEY: The level of representation on which this 
                               scoring will act on.
                     SCORE: The value of the evaluated score.
       returns: The key index                                
    """
    for k in sorted(RestraintDict.keys()):
      if RestraintDict[k]['NAME'] == keyName:
        keyIndex = k;
    return keyIndex

  def generate_domain_pairs(self):
    """
      Helper method for setting up restraint calculation to generate all
      pairwise domain pairs. 
    """
    domOpt = self.opt["DOMAIN"]
    ndom = domOpt["NUMBER"]
    pairs = []
    for i in range(ndom):
      for j in range(i,ndom):
        pairs.append((i,j))
    self.opt["DOMAIN"]["DOMAIN_PAIR"] = pairs

  def evaluate_restraints(self,representation):
    """
      Given input representation object of instance Representation
      this evaluates the score from the restraint dictionary saved in class
      variable.
      representation: Input representation of type class Representation
                      Currently this representation includes domain-level and
                      /or residue-level ellipsoids
      returns: Restraint dictionary updated to include score value
    """
    Restr = copy.copy(self.Restr)
    for keyIndex in sorted(Restr.keys()):
      robj = Restr[keyIndex]
      if robj['EVAL']:
        inputkey = robj['INPUTKEY']
        if inputkey == 'ellipsoid':
          if robj['NAME'] == 'DOMAIN_PAIR' or robj['NAME'] == 'DOMAIN_EXCL':
            robj['SCORE'] = robj['RESTR'].evaluate(representation.DomEllipsoids)
          else: robj['SCORE'] = 0.0
        else:
          robj['SCORE'] = 0.0
      else: robj['SCORE'] = 0.0
      Restr[keyIndex] = robj
    return Restr

  def get_evaluated_scores(self,Restr):
    """
      Calculate sum of scores as included in restraint dictionary as input
      Restr: Input restraint dictionary
      returns: sum of scores evaluated
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
      Calculate sum of scores as included in restraint dictionary as input
      as well print the score on-screen
      Restr: Input restraint dictionary
      returns: sum of scores evaluated
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

  def setup_data(self):
    """
      Setting data dictionary to be used in MCSampler class. The data
      dictionary includes following information with the respective
      key-strings:
        ellipsoid: A dictionary with keys 0 or 1 for array of residue-
                   ellipsoids or array of domain-ellipsoids respectively
        mesh: Evaluated mesh from ellipsoids
        targetOpt: Set of options to work with target image (projected
                   image from coordinate data)
        refernOpt: Set of options to work with reference image (read input
                   image from AFM-data)
        refernMat: Reference image matrix
        restr: The dictionary of restraints
    """
    ellipsoid = {}
    if self.opt["ELLIPSOID"]["NLEVEL"] == 1:
      if self.opt["ELLIPSOID"]["LEVEL"][0] == 'domain':
        ellipsoid = {1: self.representation.DomEllipsoids}
      elif self.opt["ELLIPSOID"]["LEVEL"][0] == 'residue':
        ellipsoid = {0: self.representation.ResEllipsoids}
    elif self.opt["ELLIPSOID"]["NLEVEL"] == 2:
      ellipsoid = {1: self.representation.DomEllipsoids, 0: self.representation.ResEllipsoids}    
    #if self.opt["MODE"] == "Debug":
    #  atoms = [a for a in self.pdb.object.get_atoms()]
    #else: atoms = None
    return {"ellipsoid": ellipsoid, "mesh": self.representation.mesh,
            "targetOpt": self.targetOpt, "refernOpt": self.referenceOpt, 
            "refernMat": self.referenceMat, "restr": self.Restr}

  def reset_out_directory_runtime(self,out):
    """
      A convenient method to reset out-directory as mentioned in user-given
      input yaml data by the input string. This input string is the string 
      obtained from parsing the argument of executable (afmSampler)
      out: Output string directory
    """
    if out != self.opt['OPTIMIZE']['OUTPUT']['IMAGEPATH']:
      self.opt['OPTIMIZE']['OUTPUT']['IMAGEPATH'] = out
    if out != self.opt['OPTIMIZE']['OUTPUT']['STATPATH']:
      self.opt['OPTIMIZE']['OUTPUT']['STATPATH']  = out
    if out != self.opt['OPTIMIZE']['OUTPUT']['ELLIPPATH']:
      self.opt['OPTIMIZE']['OUTPUT']['ELLIPPATH'] = out
    if out != self.opt['OPTIMIZE']['OUTPUT']['GMDATPATH']:
      self.opt['OPTIMIZE']['OUTPUT']['GMDATPATH'] = out
    if out != self.opt['OPTIMIZE']['TRJFILE']:
      self.opt['OPTIMIZE']['TRJFILE'] = "%s.trj5" % out

  def sample(self,out="out",restartFileInput=None,restartFileOutput=".dump"):
    """
      Start sampling using class MCSampler
      out: Output string directory as obtained from parsing the argument of
           executable (afmSampler)
      restartFileInput: None for no restart. Else a string which will be
                        used to load MCSampler instance last used
      restartFileOutput: A file for dumping MCSampler instance at the 
                         completion of sampling, which by defualt is ".dump"
      Note: Restarting sampling recommends one important user work. Before
            restarting the user is requested to take back-up of STATFILE
            present within the directory given in 'out'-string argument. This
            is required if user is interested to rerun the restarted sampling.
            For continued restarted sampling (sampling for n1 steps, followed
            by n2 steps, ..., followed by nk steps,...) such backup not 
            required, although recommended.
    """
    self.reset_out_directory_runtime(out)
    if restartFileInput is not None:
      with open(restartFileInput, 'rb') as fp:
        mcs = dill.load(fp)
      mcs.restart(self.opt)
    else:
      data = self.setup_data()
      mcs = MCSampler(self.opt,data,timing=False)
      mcs.start()
    mcs.finish(restartFileOutput)          

