#!/usr/bin/env python

#get library modules
import os, sys, string

curDir = os.getcwd() + "/"

argc = len(sys.argv)

def helpMenu():
  print "Use following arguments:"
  print " >run3dna -mbp [inFolder] [outFolder]"
  print " >run3dna -outp [inFolder] [outFolder]"
  print " >run3dna -dssr [inFolder] [outFolder]"


def fileExists(pathToFile):
  try: 
    with open(pathToFile) as f: return(True)
  except IOError as e: return(False)    

def removeJunk():
    #cleanup
  if fileExists(curDir + "allpairs.ana") == True:
    os.system("rm %s" % (curDir + "allpairs.ana"))
  if fileExists(curDir + "auxiliary.par") == True:
    os.system("rm %s" % (curDir + "auxiliary.par"))
  if fileExists(curDir + "bp_helical.par") == True:
    os.system("rm %s" % (curDir + "bp_helical.par"))
  if fileExists(curDir + "bp_step.par") == True:
    os.system("rm %s" % (curDir + "bp_step.par"))
  if fileExists(curDir + "cf_7methods.par") == True:
    os.system("rm %s" % (curDir + "cf_7methods.par"))
  if fileExists(curDir + "hstacking.pdb") == True:
    os.system("rm %s" % (curDir + "hstacking.pdb"))
  if fileExists(curDir + "mref_frames.dat") == True:
    os.system("rm %s" % (curDir + "mref_frames.dat"))
  if fileExists(curDir + "mulbp.inp") == True:
    os.system("rm %s" % (curDir + "mulbp.inp"))
  if fileExists(curDir + "multiplets.pdb") == True:
    os.system("rm %s" % (curDir + "multiplets.pdb"))
  if fileExists(curDir + "poc_haxis.r3d") == True:
    os.system("rm %s" % (curDir + "poc_haxis.r3d"))
  if fileExists(curDir + "ref_frames.dat") == True:
    os.system("rm %s" % (curDir + "ref_frames.dat"))
  if fileExists(curDir + "stacking.pdb") == True:
    os.system("rm %s" % (curDir + "stacking.pdb"))
  if fileExists(curDir + "tmp_file") == True:
    os.system("rm %s" % (curDir + "tmp_file"))
  if fileExists(curDir + "allpairs.pdb") == True:
    os.system("rm %s" % (curDir + "allpairs.pdb"))
    #DSSR files
  if fileExists(curDir + "dssr-hairpins.pdb") == True:
    os.system("rm %s" % (curDir + "dssr-hairpins.pdb"))
  if fileExists(curDir + "dssr-helices.pdb") == True:
    os.system("rm %s" % (curDir + "dssr-helices.pdb"))
  if fileExists(curDir + "dssr-multiplets.pdb") == True:
    os.system("rm %s" % (curDir + "dssr-multiplets.pdb"))
  if fileExists(curDir + "dssr-pairs.pdb") == True:
    os.system("rm %s" % (curDir + "dssr-pairs.pdb"))
  if fileExists(curDir + "dssr-stems.pdb") == True:
    os.system("rm %s" % (curDir + "dssr-stems.pdb"))
  if fileExists(curDir + "dssr-Bfactors.dat") == True:
    os.system("rm %s" % (curDir + "dssr-Bfactors.dat"))

if argc != 4:
  helpMen()

else:
  inDir = curDir + sys.argv[2]
  outDir = curDir + sys.argv[3]

  if not os.path.exists(inDir): os.makedirs(inDir) 
  if not os.path.exists(outDir): os.makedirs(outDir) 
  mbpCmd = "find_pair -p -original_coordinate %s.pdb %s.mbp" 
  outpCmd = "analyze -c allpairs.ana"   
  dssrCmd = "x3dna-dssr --input=%s.pdb -o=%s.dssr"
  lsCmd = "ls %s > tempLs.txt"

  pdbId = []
  noPairs = []
  hasPairs = []

  if sys.argv[1] == "-mbp":
    removeJunk()

  elif sys.argv[1] == "-outp":
    removeJunk()

  elif sys.argv[1] == "-dssr":
    os.system(lsCmd % inDir)
    FILE = open((curDir + "tempLs.txt"), "rU")
    for i in FILE:
      pdbId.append(i.strip().strip(".pdb"))
    FILE.close

    for i in pdbId:
      print "----Processing:", i, "-------"
      os.system(dssrCmd % ((inDir + i), (outDir + i)))
      print ""
      firstLine = ""
      with open((outDir + i + ".dssr"), "rU") as f:
        firstLine = f.readline()
      if "List of" in firstLine:
        hasPairs.append(i)
      else:
        print "Found no basepairs for:", i
        noPairs.append(i)

      FILE = open((curDir + "no_pairs.txt"), "wb")
      for z in noPairs:
        FILE.write(z.upper() + "\n")
      FILE.close

      FILE = open((curDir + "has_pairs.txt"), "wb")
      for z in hasPairs:
        FILE.write(z.upper() + "\n")
      FILE.close
    removeJunk()

