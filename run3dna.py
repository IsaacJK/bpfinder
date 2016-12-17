#!/usr/bin/env python

#get library modules
import os, sys, string

class Run3DNA:
  def __init__(self):
    self.curDir = os.getcwd() + "/"
    self.realDir = os.path.dirname(os.path.realpath(__file__)) + "/"
    self.pdbPath = ""
    self.mbpDir = self.realDir + "mbpFiles/"
    self.outpDir = self.realDir + "outpFiles/"
    self.dssrDir = self.realDir + "dssrFiles/"
    self.findPair = ""
    self.analyze = ""
    #create array to store filenames
    self.no_pairs = set([])
  
  def fileExists(self, pathToFile):
    try: 
      with open(pathToFile) as f: return(True)
    except IOError as e: return(False)    
  
  # Parse PDB files with 3DNA/DSSR
  # Purge overwrites existing file exclusion
  def genFiles(self, pdbPath, names, purge=False):
    self.pdbPath = pdbPath
    self.names = names
      #create Dir to store files
    if not os.path.exists(self.mbpDir): os.makedirs(self.mbpDir) 
    if not os.path.exists(self.outpDir): os.makedirs(self.outpDir)
    if not os.path.exists(self.dssrDir): os.makedirs(self.dssrDir)        
    self.findPair = "find_pair -p -original_coordinate %s.pdb %s.mbp" 
    self.dssr = "x3dna-dssr --long-idstr -i=%s.pdb -o=%s.dssr"
    for i in self.names:
        #Check if MBP file already exists
      temp_mbp_path = os.path.join(self.mbpDir, i + ".mbp")
      if self.fileExists(temp_mbp_path) == False or purge == True:
        os.system(self.findPair % ((self.pdbPath + i),
                                   (self.mbpDir + i)))
        os.system("analyze -c allpairs.ana")
        os.system(self.dssr % ((self.pdbPath + i),
                          (self.dssrDir + i)))
        temp_outp_path = os.path.join(i + ".outp")
        if self.fileExists(temp_outp_path):
          os.system("mv %s.outp %s" % (i, self.outpDir))

      #cleanup
    if self.fileExists(self.curDir + "allpairs.ana") == True:
      os.system("rm %s" % (self.curDir + "allpairs.ana"))
    if self.fileExists(self.curDir + "auxiliary.par") == True:
      os.system("rm %s" % (self.curDir + "auxiliary.par"))
    if self.fileExists(self.curDir + "bp_helical.par") == True:
      os.system("rm %s" % (self.curDir + "bp_helical.par"))
    if self.fileExists(self.curDir + "bp_step.par") == True:
      os.system("rm %s" % (self.curDir + "bp_step.par"))
    if self.fileExists(self.curDir + "cf_7methods.par") == True:
      os.system("rm %s" % (self.curDir + "cf_7methods.par"))
    if self.fileExists(self.curDir + "hstacking.pdb") == True:
      os.system("rm %s" % (self.curDir + "hstacking.pdb"))
    if self.fileExists(self.curDir + "mref_frames.dat") == True:
      os.system("rm %s" % (self.curDir + "mref_frames.dat"))
    if self.fileExists(self.curDir + "mulbp.inp") == True:
      os.system("rm %s" % (self.curDir + "mulbp.inp"))
    if self.fileExists(self.curDir + "multiplets.pdb") == True:
      os.system("rm %s" % (self.curDir + "multiplets.pdb"))
    if self.fileExists(self.curDir + "poc_haxis.r3d") == True:
      os.system("rm %s" % (self.curDir + "poc_haxis.r3d"))
    if self.fileExists(self.curDir + "ref_frames.dat") == True:
      os.system("rm %s" % (self.curDir + "ref_frames.dat"))
    if self.fileExists(self.curDir + "stacking.pdb") == True:
      os.system("rm %s" % (self.curDir + "stacking.pdb"))
    if self.fileExists(self.curDir + "tmp_file") == True:
      os.system("rm %s" % (self.curDir + "tmp_file"))
    if self.fileExists(self.curDir + "allpairs.pdb") == True:
      os.system("rm %s" % (self.curDir + "allpairs.pdb"))
      #DSSR files
    if self.fileExists(self.curDir + "dssr-hairpins.pdb") == True:
      os.system("rm %s" % (self.curDir + "dssr-hairpins.pdb"))
    if self.fileExists(self.curDir + "dssr-helices.pdb") == True:
      os.system("rm %s" % (self.curDir + "dssr-helices.pdb"))
    if self.fileExists(self.curDir + "dssr-multiplets.pdb") == True:
      os.system("rm %s" % (self.curDir + "dssr-multiplets.pdb"))
    if self.fileExists(self.curDir + "dssr-pairs.pdb") == True:
      os.system("rm %s" % (self.curDir + "dssr-pairs.pdb"))
    if self.fileExists(self.curDir + "dssr-stems.pdb") == True:
      os.system("rm %s" % (self.curDir + "dssr-stems.pdb"))
    if self.fileExists(self.curDir + "dssr-Bfactors.dat") == True:
      os.system("rm %s" % (self.curDir + "dssr-Bfactors.dat"))
    if self.fileExists(self.curDir + "dssr-torsions.dat") == True:
      os.system("rm %s" % (self.curDir + "dssr-torsions.dat"))
      #Return directory where mbp/outp and pair pdb's are stored
    return(self.mbpDir, self.outpDir, self.dssrDir)
