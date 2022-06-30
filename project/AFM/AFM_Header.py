try:
    import readline
except ImportError:
    print("Module readline not available.")
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")

import io,struct,sys,time
import random
import copy
import os
from os.path import splitext,isfile,exists
import math
from math import pi,log,sqrt,exp,ceil,floor
from operator import itemgetter
from threading import Thread
from multiprocessing import Process
from datetime import datetime
import StringIO
import re

import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.stats import multivariate_normal,chi2
from scipy.stats.stats import pearsonr
from scipy import ndimage
from scipy.spatial.distance import cdist

import Bio.PDB

# import afmtools

from PIL import Image
from skimage import filters,morphology
from skimage import measure

from markdown import markdown
import pkg_resources
import yaml
import dill

import h5py

class OptionAFM:  
  """
    Class for container of all the options
  """
  _LIMIT_TRANSLATION_ = 2.0
  _LIMIT_ROTATION_    = 10.0
  _LIMIT_SCALE_ = 0.5
  def __init__(self):
    """
      Constructor method to initialize with default options
    """
    self._LIMIT_TRANSLATION_ = OptionAFM._LIMIT_TRANSLATION_
    self._LIMIT_ROTATION_ = OptionAFM._LIMIT_ROTATION_
    self._LIMIT_SCALE_ = OptionAFM._LIMIT_SCALE_
    path = 'data.yml'
    filepath = pkg_resources.resource_filename(__name__, path)    
    with open(filepath, 'r') as stream:
      try:
        self.opt = yaml.load(stream)
      except yaml.YAMLError as exc:
        print(exc)

  def get_options(self):
    """
      return saved options
    """
    return self.opt

  def set_options(self,opt):
    """
      sets options from input options
      opt: Input option
    """
    self.opt = opt

  def reset_options(self,opt):
    """
      The method to feed in input option and reset default options
      by the input options
      opt: Input options
      returns: number of error generated if CHECK_OPTION key in default
               option is True
    """
    cur_opt = opt
    old_opt = self.opt
    c = self.compare_options(cur_opt,old_opt) 
    if c == "":
      print("Exact")
    else: 
      print("Importing options from default to the current options")      
      mod_opt = self.replace_options(cur_opt,old_opt)      
      self.set_options(mod_opt)
    self.expand_segment()
    if self.opt["CHECK_OPTION"]:
      num_err = self.check_options()
    else: num_err = 0
    return num_err      

  def compare_options(self,o1,o2,path=""):
    """
      Recursive method to compare two option dictionaries
      o1,o2: Input options to compare
    """
    err = ''; key_err = ''; val_err = '';
    old_path = path;
    o1_name = 'Current'; o2_name = "Previous";
    for k in o1.keys():
      path = old_path + "[%s]" % k
      if not o2.has_key(k):
        key_err += "Key %s%s not in %s\n" % (o1_name,path,o2_name)
      else:
        if isinstance(o1[k],dict) and isinstance(o2[k],dict):
          err += self.compare_options(o1[k],o2[k],path)
        else:
          if o1[k] != o2[k]:
            val_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                       % (o1_name,path,o1[k],o2_name,path,o2[k])
    for k in o2.keys():
      path = old_path + "[%s]" % k
      if not o1.has_key(k):
        key_err += "Key %s%s not in %s\n" % (o2_name,path,o1_name)
    return key_err + val_err + err

  def replace_options(self,curr,prev):
    """
      Recursive method to replace current dictionary by previous dictionary
    """
    for k in prev.keys():
      if not curr.has_key(k):
        curr[k] = copy.deepcopy(prev[k])
      else:
        curr_val = curr[k]; prev_val = prev[k];
        if isinstance(curr_val,dict) and isinstance(prev_val,dict):
          curr[k] = self.replace_options(curr_val,prev_val)
    return curr

  @staticmethod
  def expand_to_int_list(L):    
    """
      Expand packed string of integers to list of integers
      returns: expanded list of integers
    """
    elm = L.split(",")
    outList = []
    for e in elm:
      elist = []
      if e.find("-") != -1:
        ebeg = int(e.split("-")[0])
        eend = int(e.split("-")[1])+1
        elist = range(ebeg,eend)
        outList.extend(elist)
      else: 
        outList.append(int(e))
    return sorted(list(set(outList)))

  def expand_segment(self):
    """
      Make domain definition from packed string of integers to residue indices
    """
    if self.opt.has_key("DOMAIN"):
      domain_def = self.opt["DOMAIN"]
      if domain_def.has_key("SEGMENT"):
        segments = domain_def["SEGMENT"]
        if isinstance(segments,list):
          nseg = len(segments)
          domain_def["NUMBER"] = nseg
          for i,s in enumerate(segments):
            segments[i] = self.expand_to_int_list(s)
        else:
          print("Domain segments must be a list")
        self.opt["DOMAIN"]["SEGMENT"] = copy.deepcopy(segments)

  def check_options(self):
    """
      Looking for common options set-up
      returns: number of errors
    """
    total_err = 0
    if self.opt.has_key("PDBFILE"):
      pdbfile = self.opt["PDBFILE"]
      if pdbfile == 'unk':
        print("option %-30s is set ... ... ... [Fail]" % "PDBFILE")
        total_err += 1
      else:
        print("option %-30s is set ... ... ... [Pass]" % "PDBFILE")        
    if self.opt.has_key("IMGFILE"):
      imgfile = self.opt["IMGFILE"]
      if imgfile == 'unk':
        print("option %-30s is set ... ... ... [Fail]" % "IMGFILE")
        total_err += 1
      else:
        print("option %-30s is set ... ... ... [Pass]" % "IMGFILE")
    if self.opt.has_key("ELLIPSOID"):
      nlevel = 0
      if self.opt["ELLIPSOID"].has_key("NLEVEL"):
        nlevel = self.opt["ELLIPSOID"]["NLEVEL"]
        try:
          nlevel = int(nlevel)
          print("option %-30s is set ... ... ... [Pass]" % "ELLIPSOID/NLEVEL")
        except ValueError:
          total_err += 1
          print("option %-30s is set ... ... ... [Fail]" % "ELLIPSOID/NLEVEL")
      else: 
        print("option %-30s is set ... ... ... [Fail]" % "ELLIPSOID/NLEVEL")
        total_err += 1
      if self.opt["ELLIPSOID"].has_key("LEVEL"):
        levels = self.opt["ELLIPSOID"]["LEVEL"]
        err = 0
        if isinstance(levels,list):
          n = len(levels)
          if n == nlevel:
            for l in levels:
              if l == 'unk': err = 1
          else: err = 1
        else: err = 1
        if err > 0:
          print("option %-30s is set ... ... ... [Fail]" % "ELLIPSOID/LEVEL")
        else:
          print("option %-30s is set ... ... ... [Pass]" % "ELLIPSOID/LEVEL")
        total_err += err
      if self.opt["ELLIPSOID"].has_key("SCALESD"):
        scalesd = self.opt["ELLIPSOID"]["SCALESD"]
        err = 0
        if isinstance(scalesd,list):
          n = len(scalesd)
          if n == nlevel:
            for s in scalesd:
              try:
                s = float(s)
              except ValueError:
                err = 1
          else: err = 1
        else: err = 1
        if err > 0:
          print("option %-30s is set ... ... ... [Fail]" % "ELLIPSOID/SCALESD")
        else:
          print("option %-30s is set ... ... ... [Pass]" % "ELLIPSOID/SCALESD")
        total_err += err
      if self.opt["ELLIPSOID"].has_key("WEIGHT"):
        wt = self.opt["ELLIPSOID"]["WEIGHT"]
        err = 0
        if isinstance(wt,list):
          n = len(wt)
          if n == nlevel:
            for w in wt:
              if w != 'mass' and w != 'volume': 
                err = 1
          else: err = 1
        else: err = 1
        if err > 0:
          print("option %-30s is set ... ... ... [Fail]" % "ELLIPSOID/WEIGHT")
        else:
          print("option %-30s is set ... ... ... [Pass]" % "ELLIPSOID/WEIGHT")
        total_err += err
    if self.opt.has_key("MESH"):
      if self.opt["MESH"].has_key("LEVEL"):
        if self.opt["MESH"]["LEVEL"] == 'unk':
          print("option %-30s is set ... ... ... [Fail]" % "MESH/LEVEL")
          total_err += 1
        elif self.opt["MESH"]["LEVEL"] != 'domain' and self.opt["MESH"]["LEVEL"] != 'residue':
          print("option %-30s is set ... ... ... [Fail]" % "MESH/LEVEL")
          total_err += 1
        else: 
          print("option %-30s is set ... ... ... [Pass]" % "MESH/LEVEL")
    if self.opt.has_key("MODE"):
      mode = self.opt['MODE']
      if mode != 'Debug' and mode != 'Production':
        self.opt['MODE'] = 'Debug'
    else: self.opt['MODE'] = 'Debug'
    if self.opt.has_key("NTHR"):
      nthr = self.opt['NTHR']
      try:
        nthr = int(nthr)
      except ValueError:
        print("NTHR improper, can't cast to int, falling back to 1")
        nthr = 1
      self.opt['NTHR'] = nthr
    else: self.opt['NTHR'] = 1
    if self.opt["IMGCOMPARE"]["CORRELATION"]["Q"] == 'none':
      self.opt["IMGCOMPARE"]["CORRELATION"]["Q"] = None
    if self.opt["IMGCOMPARE"]["PLOTCORR"]["Q"] == 'none':
      self.opt["IMGCOMPARE"]["PLOTCORR"]["Q"] = None
    return total_err
