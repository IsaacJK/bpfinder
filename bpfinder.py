'''
Program: Basepair Finder
Author: Isaac J Kimsey
'''

# StdLib modules
import os, sys, string
from time import clock, time
import datetime
# Custom modules
import control, updateBPF, menu
import cPickle as pickle
import argparse

# Start clock to keep track of run-time
startTime = clock()
argc = len(sys.argv)

  #Define the path of where the program is
  #This will be used to store and read data files in that local
  #directory.
realPath = os.path.dirname(os.path.realpath(__file__)) + "/"
  #Folder that stores all the .dat files that are needed by this program
datDir = os.path.join(realPath, "dat/")

err_logp = os.path.join(datDir, 'ErrorLog-{:%Y-%m-%d %Hh%Mm%Ss}.txt'.format(datetime.datetime.now()))
  #Current working directory
curDir = os.getcwd()

  #Define directories of the various dat files
namepath = os.path.join(datDir, "names.dat")
pairspath = os.path.join(datDir, "pairs.dat")
pdbdatapath = os.path.join(datDir, "pdbData.dat")
hairpinpath = os.path.join(datDir, "hairpins.dat")

def printHelp():
  print "  Requires curl, matplotlib and numpy"
  print "  python bpfinder.py -h"
    #PDB input file should be optional, if input PDB not given just
    #default to base-search by non-geometric parameters like
    #bonding atoms, pucker, synanti conf, pair motifs, etc.
    #These should be pickled compactly as well so that they can be
    #accessed rapidly
  print "  python bpfinder.py -run"
  print "  python bpfinder.py -up names.dat"  

def checkDependencies():
  #Depends on external modules:
  #BPFinder based:
    #updateBPF
    #menu
    #control
    #pairstatistics
    #parse
    #pymolintegration
    #run3dna
    #strucs
    #updateBPF
    #clustStruc
  #Others:
    #numpy
    #matplotlib
    #mpl_toolkits : Axes3D
  depCheck = 1

  try: import numpy
  except ImportError:  
    print "Failed to import: numpy"
    depCheck = 0

  try: import matplotlib.pyplot
  except ImportError:  
    print "Failed to import: matplotlib.pyplot"
    depCheck = 0

  try: import matplotlib.mlab
  except ImportError:  
    print "Failed to import: matplotlib.mlab"
    depCheck = 0

  try: import mpl_toolkits.mplot3d
  except ImportError:  
    print "Failed to import: mpl_toolkits.mplot3d"
    depCheck = 0

  try: import matplotlib.backends.backend_pdf
  except ImportError:  
    print "Failed to import: matplotlib.backends.backend_pdf "
    depCheck = 0

  try: import updateBPF
  except ImportError:  
    print "Failed to import: updateBPF"
    depCheck = 0

  try: import menu
  except ImportError:  
    print "Failed to import: menu"
    depCheck = 0

  try: import control
  except ImportError:  
    print "Failed to import: control"
    depCheck = 0

  try: import pairstatistics
  except ImportError:  
    print "Failed to import: pairstatistics"
    depCheck = 0

  return(depCheck)

  #Utility function that, when called, checks if a given file exists
  #and returns a corresponding boolean value
def fileExists(pathToFile):
  try: 
    with open(pathToFile) as f: return(True)
  except IOError as e: return(False)    

if checkDependencies() == 0:
  print "Failed dependency check"

else:
  if argc <= 1:
    print "Too few input parameters."
    printHelp()
  else:
    ctrl = control.Control()

    if sys.argv[1] == "-h" or sys.argv[1] == "--h": printHelp() 

      #Builds up the names.dat file based on user defined params
      #Basically the same as the -up function, but only exports names
    elif sys.argv[1] == "-names":
      opt1, opt2 = 0, 0
      lores, hires, tlo = 0.0, 3.0, 0.0
      upXray = "Y"
      upNmr = "Y"
        #Can update xray or NMR structures
      while True:
        try:
          print "\nWhat PDB names you like to add to the names.dat file?"
          print " 1. [%s] X-Ray Structures ([Y]es/[N]o)"  % upXray
          print " 2. [%s - %s] X-Ray resolution" %  (lores, hires)
          print " 3. [%s] NMR Structures ([Y]es/[N]o)" % upNmr
          print " 9. Reset Parameters"
          print " 0. Save and Return"
          opt1 = int(raw_input(" > "))
          if opt1 == 1:
            print ""
            while True:
              print "Include X-Ray Structures?"
              print " 1. Yes"
              print " 2. No"
              print " 3. Return to previous menu"
              try: 
                opt2 = int(raw_input(" > "))
                if opt2 == 3:
                  break
                if opt2 == 1:
                  upXray = "Y"
                  break
                if opt2 == 2:
                  upXray = "N"
                  break
                else: print "\nSelect a value between 1 and 3.\n"
              except ValueError: print "Please select an integer value.\n"            
 
          if opt1 == 2:
            while True:
              print "Select a resolution range."
              try:
                print " Minimum resolution"
                lores = float(raw_input(" > "))
                print " Maximum resolution"
                hires = float(raw_input(" > "))

                lores = round(lores, 2)
                hires = round(hires, 2)

                if lores < 0 or hires < 0:
                  print "Please enter a non-negative real number\n"

                elif lores > hires:
                  tlo = lores
                  lores = hires
                  hires = tlo
                  break
                else: break
              except ValueError: print "Please select a real number.\n"

          if opt1 == 3:
            print ""
            while True:
              print "Include NMR Structures?"
              print " 1. Yes"
              print " 2. No"
              print " 3. Return to previous menu"
              try: 
                opt2 = int(raw_input(" > "))
                if opt2 == 3:
                  break
                if opt2 == 1:
                  upNmr = "Y"
                  break
                if opt2 == 2:
                  upNmr = "N"
                  break
                else: print "\nSelect a value between 1 and 3.\n"
              except ValueError: print "Please select an integer value.\n"    

          if opt1 == 9: 
            lores, hires, tlo = 0.0, 3.0, 0.0
            upXray = "Y"
            upNmr = "Y"
                     
          if opt1 == 0: 
              #Just set resolution to hi-res value
            resolution = hires
            if fileExists(namepath) == False:
              print "names.dat does not exist."
            else:
              if upXray == "Y" and upNmr == "N":
                ctrl.buildNames(namepath, resolution, "XRAY")
              elif upXray == "N" and upNmr == "Y":
                ctrl.buildNames(namepath, resolution, "NMR")
              elif upXray == "Y" and upNmr == "Y":
                ctrl.buildNames(namepath, resolution, "BOTH")           
            break
        except ValueError: print "Please select an integer value."

      #Menu for updating pair data
    elif sys.argv[1] == "-up":
      opt1, opt2 = 0, 0
      lores, hires, tlo = 0.0, 3.0, 0.0
      upXray = "Y"
      upNmr = "Y"
      forceup = "N"
      purge = False # Purge flag
        #Can update xray or NMR structures
      while True:
        try:
          print "\nWhat would you like to update?"
          print " 1. [%s] X-Ray Structures ([Y]es/[N]o)"  % upXray
          print " 2. [%s - %s] X-Ray resolution" %  (lores, hires)
          print " 3. [%s] NMR Structures ([Y]es/[N]o)" % upNmr
          print " 4. [%s] Purge and force rebuild ([Y]es/[N]o)" % forceup
          print " 9. Reset Parameters"
          print " 0. Save and Return"
          opt1 = int(raw_input(" > "))
          if opt1 == 1:
            print ""
            while True:
              print "Include X-Ray Structures?"
              print " 1. Yes"
              print " 2. No"
              print " 3. Return to previous menu"
              try: 
                opt2 = int(raw_input(" > "))
                if opt2 == 3:
                  break
                if opt2 == 1:
                  upXray = "Y"
                  break
                if opt2 == 2:
                  upXray = "N"
                  break
                else: print "\nSelect a value between 1 and 3.\n"
              except ValueError: print "Please select an integer value.\n"            
 
          if opt1 == 2:
            while True:
              print "Select a resolution range."
              try:
                print " Minimum resolution"
                lores = float(raw_input(" > "))
                print " Maximum resolution"
                hires = float(raw_input(" > "))

                lores = round(lores, 2)
                hires = round(hires, 2)

                if lores < 0 or hires < 0:
                  print "Please enter a non-negative real number\n"

                elif lores > hires:
                  tlo = lores
                  lores = hires
                  hires = tlo
                  break
                else: break
              except ValueError: print "Please select a real number.\n"

          if opt1 == 3:
            print ""
            while True:
              print "Include NMR Structures?"
              print " 1. Yes"
              print " 2. No"
              print " 3. Return to previous menu"
              try: 
                opt2 = int(raw_input(" > "))
                if opt2 == 3:
                  break
                if opt2 == 1:
                  upNmr = "Y"
                  break
                if opt2 == 2:
                  upNmr = "N"
                  break
                else: print "\nSelect a value between 1 and 3.\n"
              except ValueError: print "Please select an integer value.\n"  

          # Force purge of all old database items and redownload
          if opt1 == 4:
              # Blank out names and pairs / etc
              ctrl = control.Control()
              ctrl.buildBlank(pairspath, namepath, hairpinpath)            
              # Change purge flag
              purge = True

          if opt1 == 9: 
            lores, hires, tlo = 0.0, 3.0, 0.0
            upXray = "Y"
            upNmr = "Y"
                     
          if opt1 == 0: 
              #Just set resolution to hi-res value
            resolution = hires
            if fileExists(namepath) == False or fileExists(pairspath) == False:
              print "names.dat or pairs.dat does not exist."
            elif fileExists(pdbdatapath) == False:
              print "pdbData.dat does not exist."
            elif fileExists(hairpinpath) == False:
              print "hairpins.dat does not exist."
            else:
              if upXray == "Y" and upNmr == "N":
                ctrl.fullUp(namepath, pairspath, resolution, pdbdatapath,
                                  hairpinpath, "XRAY", datDir, purge=purge)
              elif upXray == "N" and upNmr == "Y":
                ctrl.fullUp(namepath, pairspath, resolution, pdbdatapath, 
                                  hairpinpath, "NMR", datDir, purge=purge)
              elif upXray == "Y" and upNmr == "Y":
                ctrl.fullUp(namepath, pairspath, resolution, pdbdatapath,
                                  hairpinpath, "BOTH", datDir, purge=purge)                
            break
        except ValueError: print "Please select an integer value."

    elif sys.argv[1] == "-list":
      if argc == 3:
        listpath = os.path.join(curDir, sys.argv[2])
      
      ctrl = control.Control()
      ctrl.fullUp(namepath, pairspath, 3, pdbdatapath, hairpinpath, 
                        "LIST", listpath)

        #Use for dev mode.
        #WARNING erases current DAT files
        #Uses small number of PDB files to work with, to speed up loading
        #and processing time
    elif sys.argv[1] == "-dev":
      # Create an error log file for this update
      if fileExists(hairpinpath) == False:
        print "hairpins.dat does not exist."
      purge = True
      ctrl = control.Control()
      ctrl.buildBlank(pairspath, namepath, hairpinpath)
        #Tells control to tell the updater to only update with HYB XRAY strucs
        #less than 2.0 angstrom, <50 strucs.
      ctrl.fullUp(namepath, pairspath, 2.0, pdbdatapath, hairpinpath, 
                        "DEV", datDir, purge=purge, err_logp = err_logp)    

    elif sys.argv[1] == "-pdbdata":
      if len(sys.argv) == 2:     
        if fileExists(pdbdatapath) == True and fileExists(namepath) == True:
          ctrl = control.Control()
          FILE = open(namepath, "rb")
          pdbNames = pickle.load(FILE)
          FILE.close()
          print len(pdbNames)
          ctrl.buildPDBData(pdbNames, pdbdatapath)

        else: print "pdbData.dat or names.dat does not exist. Create them."

    elif sys.argv[1] == "-run":
      if len(sys.argv) == 2:
          #List of base pair information
        load_pairs = []
          #List of pdb data (IE temp, solvent, etc)
        load_pdbData = []
          #List of PDB structural information like hairpins, etc
        load_pdbStrucs = []
        if fileExists(pairspath) == True:
          FILE = open(pairspath, "rb")        
          load_pairs = pickle.load(FILE)
          FILE.close()
          if fileExists(pdbdatapath) == True:
            FILE = open(pdbdatapath, "rb")
            load_pdbData = pickle.load(FILE)
            FILE.close()
            menuob = menu.MenuC()
            menuob.mainMenu(menuob, load_pairs, 0, load_pdbData)
          else: print "pdbData.dat does not exist."
        else: print "pairs.dat does not exist."

      elif len(sys.argv) == 3:
        if ".CSV" in sys.argv[2].upper():
          load_pairs = []
          load_pdbData = []
          if (fileExists(pairspath) == True 
              and fileExists(pdbdatapath) == True):
            FILE = open(pairspath, "rb")
            load_pairs = pickle.load(FILE)
            FILE.close()
            FILE = open(pdbdatapath, "rb")
            load_pdbData = pickle.load(FILE)
            FILE.close()
            menuob = menu.MenuC()
            menuob.mainMenu(menuob, load_pairs, sys.argv[2], load_pdbData)
          else: print "pairs.dat or pdbData.dat does not exist"
    elif sys.argv[1] == "-blank":
      ctrl = control.Control()
      ctrl.buildBlank(pairspath, namepath, hairpinpath)
    else:
      print "Poor command arguments given."
      printHelp()
  # Print out run-time
  elapsedTime = (clock() - startTime)
  print elapsedTime
