import os, sys, fnmatch

curDir = os.getcwd()
argc = len(sys.argv)

def fileExists(pathToFile):
  try: 
    with open(pathToFile) as f: return(True)
  except IOError as e: return(False) 

def makeFolder(pathToFolder):
  if not os.path.exists(pathToFolder): 
    os.makedirs(pathToFolder) 

def help():
  print "Usage is as follows:"
  print " >findbad.py [mbpFolder] [outpFolder]"


def findFiles(directory, keyword):
  foundFiles, foundDirs = [], []
  for root, dirs, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, keyword):
      foundFiles.append(filename)
      foundDirs.append(os.path.join(root, filename))
  return(foundFiles, foundDirs)

if argc == 3:
  badfiles = set([])
  mbpfiles = set([])
  outpfiles = set([])

  mbpDir = os.path.join(curDir, sys.argv[1])
  outpDir = os.path.join(curDir, sys.argv[2])

  for name in os.listdir(mbpDir):
    mbpfiles.add(name[:-4])

  for name in os.listdir(outpDir):
    outpfiles.add(name[:-5])

  badfiles = (mbpfiles - outpfiles)

  FILE = open("badfiles.txt", "wb")
  for name in badfiles:
    FILE.write(name + "\n")
  FILE.close()

else: help()
