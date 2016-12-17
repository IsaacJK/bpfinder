#!/usr/bin/env python

import string, sys, os

class PDBSplitter:
  
  def __init__(self, pdbNames, pdbPath):
    #Directories
    self.curDir = os.getcwd() + "/"
    self.pdbPath = pdbPath
    self.allPDB = []
    self.split_names = []
    self.pdbNames = pdbNames
    self.pdbNames = [x.strip() for x in self.pdbNames]
    self.pdbNamesSuff = [(x + ".pdb") for x in self.pdbNames]

  def restart(self):
    sys.stdout.write('\r')
    sys.stdout.flush()
    
  def runSplitter(self):   
    for i in range(len(self.pdbNamesSuff)):
      sys.stdout.write("   %s [%d of %d]" % 
                       (self.pdbNames[i], i, len(self.pdbNamesSuff)))
      sys.stdout.flush()
      self.restart()
      self.split(self.pdbPath, self.pdbNames[i])
#        #Remove the unsplit PDB file
#      os.system("rm %s" % (self.pdbPath + self.pdbNamesSuff[i]))

  def split(self, pdbPath, pdbName):
    start, end, count = 0, 0, -1
    data, models = [], []
    start_l, end_l = [], []
    FILE = open((pdbPath + pdbName + ".pdb"), "r")
    data = [x.strip() for x in FILE]
    FILE.close    
    for i in range(len(data)):
      if data[i].find("MODEL") == 0:
        start_l.append(i + 1)
      if data[i].find("ENDMDL") == 0:
        end_l.append(i)
    for i in range(len(start_l)):
      temp = []
      for x in range(start_l[i], end_l[i]):
        temp.append(data[x] + "\n")

      models.append(list(temp))
      del temp[:len(temp)]
    for i in range(len(models)):
      self.split_names.append(pdbName + "_%d" % (i+1))
      FILE = open(pdbPath + pdbName + "_%d.pdb" % (i+1), "w")
      FILE.writelines(models[i])
      FILE.writelines("END")
      FILE.close
    

