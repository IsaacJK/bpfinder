#!/usr/bin/env python

#python copymove.py -[mv/cp] list.txt from/ to/

#get library modules
import os, sys, string

argc = len(sys.argv)

if argc == 2:
  if sys.argv[1].find("-h"):
    print "python copymove.py -[mv/cp/--h] list.txt from/ to/"
  else:
    print "Poor args"  
    
elif argc < 5:
  print "Too few args"
  
else:
  arg1 = sys.argv[1]
  arg2 = sys.argv[2]
  arg3 = sys.argv[3]
  arg4 = sys.argv[4]
  curDir = os.getcwd()
  inputDir = curDir + "/" + arg3
  outputDir = curDir + "/" + arg4  
  
  listNames = []
  FILE = open(arg2, "r")
  listNames = FILE.readlines()
  FILE.close

  for i in range(len(listNames)):
    listNames[i] = listNames[i].strip()
  
  if arg1.find("cp") > 0:
    for i in range(len(listNames)):
      cpCmd = "cp %s %s"
      cpCmd = cpCmd % ((inputDir + listNames[i]), outputDir)
      os.system(cpCmd)

  elif arg1.find("mv") > 0:
    for i in range(len(listNames)):
      mvCmd = "mv %s %s"
      mvCmd = mvCmd % ((inputDir + listNames[i]), outputDir)
      os.system(mvCmd)
    
    
