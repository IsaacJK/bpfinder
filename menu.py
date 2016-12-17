import updateBPF, menu, datetime, math, os, sys, pairstatistics, csv
import pymolintegration as bpfPym
import cPickle as pickle
import numpy as np
import Tkinter, tkFileDialog

class MenuC:
  def __init__(self):
      #Structure List Variables
    self.selStr = ["?", "?", "?", "?", "N"]
    #Sets containing PDB names pulled from the PDB
      #DNA, RNA, Protein, resolution, etc list
    self.strucsPDBlist = set([])
      #List of PDB names from user defined search of temps
    self.tempPairs = set([])
      #List of PDB names for user defined search of pH
    self.phPairs = set([])
      #List of PDB names for sequence search
    self.seqPairs = set([])
      #Parent name return list
    self.pSName = set([])
      #Ligand name return list
    self.ligPDBlist = set([])
      #General search list
    self.genSearchPDB = set([])
      #Mod resi return list
    self.modResiPDB = set([])
      #Ribosome list of strucs
    self.ribosomePDB = set([])
      #Search by chemical ID
    self.chemID = set([])
      #Structure molecular weight names
    self.searchWtPDB = set([])
    self.readInPairNames = set([])
    self.upClass = updateBPF.UpdateBPF()
    self.resLo, self.resHi, self.tempLoRes = 0.0, 0.0, 0.0 
    self.allPairs = {}
    self.pdbData = {}
    self.refPairs1, self.refPairs2, self.refPairs3 = set([]), set([]), set([])
    self.writeOutName, self.writeOutFolder = "", ""
    self.statsC = pairstatistics.pairStats()
    self.singleC = []
    self.areshiftspresent = 0
    self.shiftData = []
    self.curDir = os.getcwd() + "/"
    self.realDir = os.path.dirname(os.path.realpath(__file__)) + "/"
    self.datDir = self.realDir + "dat/"

    ########### Functional Definitions Start ###########
    #This function merges '|' the sets inside two dictionaries if they overlap
    # IE. if dict1 = {'a':set([1,2,3]), 'b':set([5,6])}
    #    and dict2 = {'b':set([1,2])}
    # Then mergeDictSets(dict1, dict2) would return:
    # returndict = {'a:set([1,2,3]), 'b':set([1,2,5,6])}
  def mergeDictSets(self, dict1, dict2):
    tdict1 = dict1
    tdict2 = dict2

    for x in dict2:
      if x in dict1:
        dict1[x] = dict1[x] | dict2[x]
      else:
        dict1[x] = dict2[x]

    return(dict1)

    ########### Functional Definitions End ###########

    ########### Main Menu Start ###########
  def mainMenu(self, mainClass, inc_allPairs, database, inc_pdbData):
    self.allPairs = inc_allPairs
    self.pdbData = inc_pdbData
    
    # for i in self.allPairs:
    #   print self.allPairs[i].uniqueID
    #   print " ", self.allPairs[i].multType

    mainClass.resetList()
    opt1, opt2 = 0, 0
    self.areshiftspresent = database
    while True:
      try:
        print "\nYou are currently searching", len(self.refPairs1), "of", len(self.allPairs), "base pairs."
        print " 1. Define Experimental Method & Molecule Type"
        print " 2. Search by PDB Query"
        print " 3. Search by Structural and Geometric Base Pair Parameters"
        print " 4. Base Pair Presets"
        print " 5. Read in pairs"
        print " 100. Undo"
        print " 101. Reset refinement list"
        print " 102. Write raw results"
        print " 103. Write raw results & Pymol files"
        print " 104. Write out raw results & statistics"
        print " 105. Write out raw results & chemical shifts"
        print " 106. Write out statistics & clusters (req. R & MClust package)"
        print " 107. Save Search Results"
        print " 108. Load Search Results"
        print " 0. Quit Program"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          mainClass.writeOutPDBNames()
          break
        if opt1 == 1:
          self.refPairs3 = self.refPairs1
          self.strucsPDBlist = set([])
          mainClass.strucList()
          mainClass.pruneListPDB(self.strucsPDBlist, "AND")
        if opt1 == 2:
          mainClass.pdbRefMenu(mainClass) 
        if opt1 == 3:
          mainClass.matchGeomM(mainClass)
        if opt1 == 4:
          mainClass.matchBPTypes(mainClass)
        if opt1 == 5:
          self.refPairs3 = self.refPairs1
          mainClass.readInPairs(mainClass)
          #Read in external pairs menu goes here
        if opt1 == 100:
          self.refPairs1 = self.refPairs3
          print ""
        if opt1 == 101:
          print ""
          mainClass.resetList()
        if opt1 == 102:
          print ""
          mainClass.writeOut()
        if opt1 == 103:
          print ""
          mainClass.writePymol(mainClass)
        if opt1 == 104:
          print ""
          mainClass.writeStats(mainClass)
        if opt1 == 105:
          self.shiftData = []
          if self.areshiftspresent == 0:
            print "\nNo chemical shift database file found\n"
          else:
            print ""
            mainClass.readShiftDB(mainClass, database)
            mainClass.writeShifts(self.shiftData)      
            self.statsC.plotStats((self.writeOutFolder + self.writeOutName), 1)
        if opt1 == 106:
          print ""
          mainClass.writeClusts(mainClass)
        if opt1 == 107:
          print ""
          mainClass.saveSearch(mainClass)
        if opt1 == 108:
          print ""
          mainClass.loadSearch(mainClass)
      except ValueError: print "Please select a valid menu option.\n"
    ########### Main Menu End ###########
  def readInPairs(self, mainC):
    opt1, opt2, opt3, cntr = 0, 0, 0, 0
    while True:
      try:
        print "You are currently searching", len(self.refPairs1), "of", len(self.allPairs), "base pairs."
        print "\nYou have read in", len(self.readInPairNames), "base pairs."
        print " 1. Read-in external pair data"
        print " 2. Include only pairs read-in"
        print " 3. Exclude pairs read-in"
        print " 4. Read-in old pair format"
        print " 8. Undo read-in pairs data"
        print " 9. Undo changes to main search pairs"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        if opt1 == 8:
          self.readInPairNames = set([])
          print ""
        if opt1 == 9:
          self.refPairs1 = self.refPairs3
          print ""
        if opt1 == 2:
          self.refPairs1 = self.refPairs1 & self.readInPairNames

        if opt1 == 3:
          self.refPairs1 = self.refPairs1 - self.readInPairNames

          #Old format
          # 1A36_B (DT 22) -- (DA 101) C
          #New format
          # 1A36_AT_B_DT_22_C_DA_101
        if opt1 == 4:
          tempNames = []
          tempstr = ""
          self.readInPairNames = set([])
          print "Input local file name"
          opt2 = str(raw_input(" > "))
          if mainC.fileExists(self.curDir + opt2) == True:
            while True:
              try:
                print " Define column containing unique pair IDs"
                opt3 = int(raw_input(" > "))
                infile = open((self.curDir + opt2), "rU")
                reader = csv.reader(infile, delimiter=",")
                for i in reader:
                  cntr += 1
                  if cntr != 1:
                    if opt3 <= len(i):
                      tempNames.append(i[opt3 - 1].strip())
                break
              except ValueError: print "Please select a valid menu option.\n"
          else: print "  File does not exist.\n"         
          for xxx in tempNames:
            temp1 = xxx.replace("_", ",")
            temp2 = temp1.replace(") -- (", ",")
            temp3 = temp2.replace(" (", ",")
            temp4 = temp3.replace(") ", ",")
            temp5 = temp4.replace(" ", ",")
            temp5 = temp5.split(",")
            tempstr = (temp5[1] + "_" + temp5[2] + "_" + temp5[3] + "_" + temp5[6]
                              + "_" + temp5[4] + "_" + temp5[5])
            temppdb = temp5[0]

            for yyy in self.refPairs1:
              if tempstr in yyy and temppdb in yyy:
                self.readInPairNames.add(yyy)

        if opt1 == 1:
          self.readInPairNames = set([])
          print "Input local file name"
          root = Tkinter.Tk()
          root.attributes('-topmost', 1)
          root.attributes('-topmost', 0)
          inLocFile = tkFileDialog.askopenfilename()
          root.destroy()
          # opt2 = str(raw_input(" > "))
          # if mainC.fileExists(self.curDir + opt2) == True:
          if mainC.fileExists(inLocFile) == True:
            while True:
              try:
                print " Define column containing unique pair IDs"
                opt3 = int(raw_input(" > "))
                infile = open(inLocFile, "rU")
                reader = csv.reader(infile, delimiter=",")
                for i in reader:
                  cntr += 1
                  if cntr != 1:
                    if opt3 <= len(i):
                      self.readInPairNames.add(i[opt3 - 1].strip())
                break
              except ValueError: print "Please select a valid menu option.\n"
          else: print "  File does not exist.\n"

      except ValueError: print "Please select a valid menu option.\n"
      
    ########### Parent Struc Type&Res Menu Start ###########
  def strucList(self):      
      #Menu option integers
    opt1, opt2 = 0, 0
    while True:
      try:
        print "\nSelect the following data for structure list refinement."
        print " 1. [%s] DNA ([Y]es/[N]o/[?]Ignore)"  % self.selStr[0]
        print " 2. [%s] RNA ([Y]es/[N]o/[?]Ignore)" % self.selStr[1]
        print " 3. [%s] Protein ([Y]es/[N]o/[?]Ignore)" % self.selStr[2]
        print " 4. [%s] Hybrid ([Y]es/[N]o/[?]Ignore)" % self.selStr[3]
        print " 5. [%s - %s] X-Ray resolution" %  (str(self.resLo), str(self.resHi))
        print " 6. [%s] NMR Structures ([Y]es/[N]o)" % self.selStr[4]
        print " 9. Reset Parameters"
        print " 0. Save and Return"
        opt1 = int(raw_input(" > "))
        if opt1 == 1:
          print ""
          while True:
            print "Include DNA in structure list?"
            print " 1. Yes"
            print " 2. No"
            print " 3. Ignore"
            print " 4. Return to previous menu"
            try: 
              opt2 = int(raw_input(" > "))
              if opt2 == 4:
                break
              if opt2 == 1:
                self.selStr[0] = "Y"
                break
              if opt2 == 2:
                self.selStr[0] = "N"
                break
              if opt2 == 3:
                self.selStr[0] = "?"
                break
              else: print "\nSelect a value between 1 and 4.\n"
            except ValueError: print "Please select an integer value.\n"
     
        if opt1 == 2:
          print ""
          while True:
            print "Include RNA in structure list?"
            print " 1. Yes"
            print " 2. No"
            print " 3. Ignore"
            print " 4. Return to previous menu"
            try: 
              opt2 = int(raw_input(" > "))
              if opt2 == 4:
                break
              if opt2 == 1:
                self.selStr[1] = "Y"
                break
              if opt2 == 2:
                self.selStr[1] = "N"
                break
              if opt2 == 3:
                self.selStr[1] = "?"
                break
              else: print "\nSelect a value between 1 and 4.\n"
            except ValueError: print "Please select an integer value.\n"

        if opt1 == 3:
          print ""
          while True:
            print "Include Protein in structure list?"
            print " 1. Yes"
            print " 2. No"
            print " 3. Ignore"
            print " 4. Return to previous menu"
            try: 
              opt2 = int(raw_input(" > "))
              if opt2 == 4:
                break
              if opt2 == 1:
                self.selStr[2] = "Y"
                break
              if opt2 == 2:
                self.selStr[2] = "N"
                break
              if opt2 == 3:
                self.selStr[2] = "?"
                break
              else: print "\nSelect a value between 1 and 4.\n"
            except ValueError: print "Please select an integer value.\n"

        if opt1 == 4:
          print ""
          while True:
            print "Include Hybrid DNA/RNA in structure list?"
            print " 1. Yes"
            print " 2. No"
            print " 3. Ignore"
            print " 4. Return to previous menu"
            try: 
              opt2 = int(raw_input(" > "))
              if opt2 == 4:
                break
              if opt2 == 1:
                self.selStr[3] = "Y"
                break
              if opt2 == 2:
                self.selStr[3] = "N"
                break
              if opt2 == 3:
                self.selStr[3] = "?"
                break
              else: print "\nSelect a value between 1 and 4.\n"
            except ValueError: print "Please select an integer value.\n"

        if opt1 == 5:
          while True:
            print "Select a resolution range."
            try:
              print " Minimum resolution"
              self.resLo = float(raw_input(" > "))
              print " Maximum resolution"
              self.resHi = float(raw_input(" > "))

              self.resLo = round(self.resLo, 2)
              self.resHi = round(self.resHi, 2)

              if self.resLo < 0 or self.resHi < 0:
                print "Please enter a non-negative real number\n"

              elif self.resLo > self.resHi:
                self.tempLoRes = self.resLo
                self.resLo = self.resHi
                self.resHi = self.tempLoRes
                break
              else: break
            except ValueError: print "Please select a real number.\n"

        if opt1 == 6:
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
                self.selStr[4] = "Y"
                break
              if opt2 == 2:
                self.selStr[4] = "N"
                break
              else: print "\nSelect a value between 1 and 3.\n"
            except ValueError: print "Please select an integer value.\n"

        if opt1 == 9: None
          # print ""
          # self.strucsPDBlist = self.strucList.clear()
        if opt1 == 0: 
          self.strucsPDBlist = self.upClass.returnList(self.selStr, self.resLo,
                                            self.resHi)
          break
      except ValueError: print "Please select an integer value."
    ########### Parent Struc Type&Res Menu End ###########

    ########### PDB Menu Start ###########
  def pdbRefMenu(self, mainClass):
    opt1, opt2 = 0, 0
    while True:
      try:
        print "\nYou are currently searching", len(self.refPairs1), "of", len(self.allPairs), "base pairs."
        print " 1. Search by PDB Name"
        print " 2. General Search"
        print " 3. Parent Structure Name"
        print " 4. Search by Chemical ID"
        print " 5. Search by Ligands"
        print " 6. Search by Modified residues"
        print " 7. Search by Temperature"
        print " 8. Search by pH"
        print " 9. Search by Sequence"
        print " 10. Search by Molecular Weight"
        print " 11. Add/Remove Large Ribosomal Subunits"
        print " 100. Undo previous changes"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        if opt1 == 100:
          self.refPairs1 = self.refPairs3
        if opt1 == 1:
          self.refPairs3 = self.refPairs1
          while True:
            try:
              print "\nChoose an option:"
              print " 1. Type in PDB names"
              print " 2. Read-in PDB list of names"
              opt2 = int(raw_input(" > "))
              if opt2 == 1:
                mainClass.mRetPDBPairs()
                break
              elif opt2 == 2:
                mainClass.readPDBlist()
                break
              else:
                print "Invalid menu option."
            except ValueError: print "Please select a valid menu option.\n"
        if opt1 == 2:
          self.refPairs3 = self.refPairs1
          self.genSearchPDB = set([])
          mainClass.genSearch()
          mainClass.pruneListPDB(self.genSearchPDB, "AND")      
        if opt1 == 3:
          self.refPairs3 = self.refPairs1
          self.pSName = set([])
          mainClass.pStrucName(mainClass)
        if opt1 == 4:
          self.refPairs3 = self.refPairs1
          self.chemID = set([])
          mainClass.findChemID(mainClass)
        if opt1 == 5:
          self.refPairs3 = self.refPairs1
          self.ligPDBlist = set([])
          mainClass.findLigands(mainClass)
        if opt1 == 6:
          self.refPairs3 = self.refPairs1
          self.modResiPDB = set([])
          mainClass.findModResi(mainClass)
        if opt1 == 7:
          self.refPairs3 = self.refPairs1
          self.tempPairs = set([])
          mainClass.tempSearchMenu(mainClass)
        if opt1 == 8:
          self.refPairs3 = self.refPairs1
          self.phPairs = set([])
          mainClass.phSearchMenu(mainClass)
        if opt1 == 9:
          self.refPairs3 = self.refPairs1
          self.seqPairs = set([])
          mainClass.searchSeq(mainClass)
        if opt1 == 10:
          self.refPairs3 = self.refPairs1
          self.searchWtPDB = set([])
          mainClass.strucWtSearch(mainClass)
        if opt1 == 11:
          self.refPairs3 = self.refPairs1
          self.ribosomePDB = set([])
          mainClass.refRibosome(mainClass)
      except ValueError: print "Please select a valid menu option.\n"
    ########### PDB Menu End ###########

    ########### Geometry Menu Start ###########
  def matchGeomM(self, mainC):
    # self.singleC = mainC.findPairValsPDB()
    opt1, opt2 = 0, 0
    while True:
      try:
        print "\nSelect which geometric parameters to search for"
        print " You are currently searching", len(self.refPairs1), "of", len(self.allPairs), "base pairs."
        print " 1. Search by Local Base Pair Parameters"
        print " 2. Match base pair motif"
        print " 3. Match H-Bond Atoms and Distances"
        print " 4. Match H-Bonding Atoms, Bases and Distances"
        print " 5. Sugar Pucker"
        print " 6. Match Specific Chi Torsions from a PDB"
        print " 7. Match Chi Torsion Range and Base Type"
        print " 8. Match by Atom Chemical Shifts"
        print " 9. Match C1'-C1' Range"
        print " 10. Match Multiplet Order"
        print " 11. Match SMILES base pairs"
        print " 12. Match Anti & Syn Parameters"
        print " 100. Undo previous changes"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        if opt1 == 100:
          self.refPairs1 = self.refPairs3
        if opt1 == 1:

          self.refPairs3 = self.refPairs1
          mainC.findLocalGeom(mainC)
        if opt1 == 2:
          self.refPairs3 = self.refPairs1
          mainC.matchMotif()
        if opt1 == 3:
          self.refPairs3 = self.refPairs1
          mainC.searchHBonds()
        if opt1 == 4:
          self.refPairs3 = self.refPairs1
          mainC.searchHBondsFine()
        if opt1 == 5:
          self.refPairs3 = self.refPairs1
          mainC.mPucker(mainC)
        if opt1 == 6:
          self.refPairs3 = self.refPairs1
          mainC.mChiTors(mainC)
        if opt1 == 7:
          self.refPairs3 = self.refPairs1
          mainC.chiRange(mainC)
        if opt1 == 8:
          self.refPairs3 = self.refPairs1
          self.shiftData = []
          if self.areshiftspresent == 0:
            print "\nNo chemical shift database file found\n"
          else:
            print ""
            mainC.readShiftDB(mainC, self.areshiftspresent)
            mainC.matchShifts(mainC)
        if opt1 == 9:
          self.refPairs3 = self.refPairs1
          mainC.mC1pC1p(mainC)
        if opt1 == 10:
          self.refPairs3 = self.refPairs1
          mainC.mMultiplet(mainC)
        if opt1 == 11:
          self.refPairs3 = self.refPairs1
          mainC.mSMILES(mainC)
        if opt1 == 12:
          self.refPairs3 = self.refPairs1
          mainC.mAntiSyn(mainC)
      except ValueError: print "Please select a valid menu option\n"
    ########### Geometry Menu End ###########
  def mAntiSyn(self, mainC):
    menuopt = 0
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    matchPairs = set([])

    while True:
      try:
        print "\nSelect a syn/anti combination"
        print " 1. Anti / Anti"
        print " 2. Anti / Syn"
        print " 3. Syn / Syn"
        print " 9. Return to previous menu without saving."
        print " 0. Return to previous menu and update pairs."
        menuopt = int(raw_input(" > "))
        if menuopt == 1:
          for i in tempList:
            if i.bbChiAconf == "Anti" and i.bbChiBconf == "Anti":
              matchPairs.add(i.uniqueID)
        elif menuopt == 2:
          for i in tempList:
            if i.bbChiAconf == "Anti" and i.bbChiBconf == "Syn":
              matchPairs.add(i.uniqueID)
            if i.bbChiAconf == "Syn" and i.bbChiBconf == "Anti":
              matchPairs.add(i.uniqueID)
        elif menuopt == 3:
          for i in tempList:
            if i.bbChiAconf == "Syn" and i.bbChiBconf == "Syn":
              matchPairs.add(i.uniqueID)
        elif menuopt == 9:
          break
        elif menuopt == 0:
          self.refPairs1 = self.refPairs1 & matchPairs
          matchPairs, tempList = set([]), []
          break
        else:
          print "Please select a valid menu option\n"
      except ValueError: print "Please select a valid menu option\n"    

  def mSMILES(self, mainC):
      #Temporary list of pairs
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    matchPairs = set()
    sMenuOpt = None
      #Numbers of found smiles, ligands and base pairs
    fSMILES, fLIGANDS, fPAIRS = 0,0,0
      #Name of CSV file with SMILES names
    tSmileFile = None
      #Individual SMILES string
    smilesstr = None
      #Column number the SMILES strings are stored in in the csv file
      #Also the number of columns in the CSV file being read in
    smilecol, numcols = 0, 0
      #List for storing the SMILES strings
    smilesdata, tempsmiles = [], []
    smilesnames, ligands = set([]), set([])
    tsmiles, tligs = set([]), set([])
      #Smile to ligand correlation
    smilelig, tsmilelig = {}, {}
      #Failed SMILES
    failsmile = set([])

    while True:
      fSMILES, fLIGANDS, fPAIRS = len(smilesnames), len(ligands), len(matchPairs)
      f404SMILE = len(failsmile)
      print ("\nFound (%s) SMILES strings corresponding to (%s) ligands across (%s) base pairs." 
             % (fSMILES, fLIGANDS, fPAIRS))
      print " Of (%s) SMILES, (%s) were not found in the PDB" % (fSMILES, f404SMILE)
      print "Select a menu option:"
      print " 1. Input SMILES string directly."
      print " 2. Input a .csv file containing multiple SMILES strings."
      print " 3. Correlate SMILES strings to base pairs."
      print " 4. Write-out SMILES corresponding Ligand ID's"
      print " 9. Return to previous menu without making any changes."
      print " 0. Save changes and return to previous menu."
      try:
        sMenuOpt = int(raw_input(" > "))
          #Search for ligands with an indv SMILES string
        if sMenuOpt == 1:
          while True:
            try:
              print "  Input SMILES string."
              smilesstr = raw_input(" > ")
                #Add the SMILES name to the list
              smilesnames.add(smilesstr)

                #Fetch the ligands corresponding to this string
              tligs, tsmilelig = self.upClass.matchSmiles(smilesstr)
              if tligs != None and tsmilelig != None:
                  #Merge the ligands found here with the previous set
                ligands = ligands | tligs
                  #Merge the sets within the dictionary of smiles for these two dicts
                smilelig = self.mergeDictSets(smilelig, tsmilelig)
              else:
                failsmile.add(smilesstr)
              break
            except ValueError: print "  Invalid input."

          #This menu option reads in a SMILES .csv file specified by the user.
        elif sMenuOpt == 2:
            #Grab name of csv file with smiles columns
          while True:
            try:
              print "\n  Input .csv file name containing SMILES."
              tSmileFile = raw_input(" > ")
              if not tSmileFile[-4:].lower() == ".csv":
                tSmileFile = tSmileFile + ".csv"
              if mainC.fileExists(os.path.join(self.curDir, tSmileFile)):
                break
              else: print "  Invalid filename, try again."
            except ValueError: print "  Invalid input."

          FILE = open(tSmileFile, "rU")
          tempsmiles = [x.strip().split(",") for x in FILE]
          FILE.close

            #Grab the smiles names and combine them with the previous list
          smilesdata = smilesdata + tempsmiles

            #Now port them over to the set containing this data to remove
            # any duplicates
          for name in smilesdata:
            smilesnames.add(name[0])

            #Grab the number of columns in the .csv file
            # this is to make sure the user does not enter
            # a column # that doesn't exist in the file
          numcols = len(smilesdata[0])

          while True:
            try:
              print "  Input column number containing SMILES IDs (ex. 1)"
              smilecol = int(raw_input(" > "))
              if smilecol <= numcols:
                  #Normalize for base zero vector ops
                smilecol = smilecol - 1
                break
              else:
                print "  Invalid column number. Try again."

            except ValueError: print "  Invalid number. Please specify an integer value."

            #Now we need to go back through all these SMILES and pull
            # the corresponding ligands
          for name in smilesnames:
              #Fetch the ligands corresponding to this string
            tligs, tsmilelig = self.upClass.matchSmiles(name)
            if tligs != None and tsmilelig != None:
                #Merge the ligands found here with the previous set
              ligands = ligands | tligs
                #Merge the sets within the dictionary of smiles for these two dicts
              smilelig = self.mergeDictSets(smilelig, tsmilelig)
            else:
              failsmile.add(name)

          #Correlate ligand ID's with base pairs
        elif sMenuOpt == 3:
          for i in tempList:
            if i.bAsubtype in ligands:
              matchPairs.add(i.uniqueID)
            if i.bBsubtype in ligands:
              matchPairs.add(i.uniqueID)

          #Write out the Ligand ID's corresponding to each SMILE
        elif sMenuOpt == 4:
          print ""
          if len(smilelig) != 0:
            now = datetime.datetime.now()
            outname = ("SMILES_Ligands_" + str(now.year) + str(now.month) 
                       + str(now.day) + str(now.hour) + str(now.minute) 
                       + str(now.second) + ".csv")
            print "Writing out SMILES and Ligand ID's to %s" % outname
            outsmile = os.path.join(self.curDir, outname)

            FILE = open(outsmile, "wb")
            for line in smilelig:
              FILE.write(line + ",")
              for lig in smilelig[line]:
                FILE.write(lig + ",")
              FILE.write("\n")
            FILE.close()
          if len(failsmile) != 0:
            now = datetime.datetime.now()
            outfail = ("Failed_SMILES_" + str(now.year) + str(now.month) 
                       + str(now.day) + str(now.hour) + str(now.minute) 
                       + str(now.second) + ".csv")    
            print "Writing out irretrievable SMILES to %s" % outfail
            FILE = open(os.path.join(self.curDir, outfail), "wb")
            for line in failsmile:
              FILE.write(line + ",\n")
            FILE.close()
          else:
            print "Nothing to write out."

          #This menu option returns the user to the previous menu without
          # making any changes to their search pairs.
        elif sMenuOpt == 9:
          break

        elif sMenuOpt == 0:
          self.refPairs1 = self.refPairs1 & matchPairs
          matchPairs = set([])
          tempList = []
          break

        else: print "Invalid menu option."

        #Exception for menu option not being a number
      except ValueError: print "  Invalid number. Please specify an integer value."
    
    # self.upClass.matchSmiles(smilesstr)


  def mMultiplet(self, mainC):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    multnum = 0
    mhigh = ""

    while True:
      try:
        print "\nEnter multiplet order to search for"
        print " (IE. '2' for a doublet)"
        multnum = int(raw_input(" > "))
        break
      except ValueError: print "Please enter a valid option.\n"    

    while True:
      try:
        print "\nSearch for", multnum, "and higher multiplets?"
        print " (Yes or No)"
        mhigh = str(raw_input(" > "))
        if mhigh.upper() == "YES" or mhigh.upper() == "Y":
          [self.refPairs2.add(x.uniqueID) for x in tempList if 
            (multnum <= x.multNum)]
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
          tempList = []
          break

        elif mhigh.upper() == "NO" or mhigh.upper() == "N":
          [self.refPairs2.add(x.uniqueID) for x in tempList if 
            (multnum == x.multNum)]
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
          tempList = []
          break
        else:
          print "Invalid option given, please specify Yes or No."
      except ValueError: print "Please enter a valid option.\n"        

  def chiRange(self, mainC):
    lo1, lo2, hi1, hi2 = 0.0, 0.0, 0.0, 0.0
    bA, bB = "", ""
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    tempList = [x for x in tempList if (x.bbChiA != "---" and x.bbChiB != "---")]
    self.refPairs2 = set([])

    while True:
      try:
        print "\nDefine Base A abbreviation (IE. A, DA, CBR, etc)"
        print " (Use * at end of residue to desinate a wild-card base)"
        print " (IE. A* will match all adenosine variants)"
        bA = str(raw_input(" Base A: ")).upper()
        break
      except ValueError: print "Invalid base type. Try again."   

    while True:
      print "\nDefine Base A - Chi Angle Lower Bound"
      print "  (Vertical bars around the chi define an absolute. Ex |90|)"
        #Define lower bounds of base A
      lo1 = raw_input(" Chi A lower bound (degrees): ")
      if lo1[:1] == "|" and lo1[-1:] == "|":
        try:
          lo1 = math.fabs(float(lo1[1:-1]))
          break
        except ValueError: print "Invalid input. Try again."   
      elif lo1[:1] != "|" and lo1[-1:] != "|":
        try:
          lo1 = float(lo1)
          break
        except ValueError: print "Invalid input. Try again."   
      else: print "Invalid input. Try again."

    while True:
      print "\nDefine Base A - Chi Angle Upper Bound"
      print "  (Vertical bars around the chi define an absolute. Ex |180|)"
        #Define upper bounds of base A
      hi1 = raw_input(" Chi A upper bound (degrees): ")
      if hi1[:1] == "|" and hi1[-1:] == "|":
        try:
          hi1 = math.fabs(float(hi1[1:-1]))
          break
        except ValueError: print "Invalid input. Try again."   
      elif hi1[:1] != "|" and hi1[-1:] != "|":
        try:
          hi1 = float(hi1)
          break
        except ValueError: print "Invalid input. Try again."   
      else: print "Invalid input. Try again."  

    while True:
      try:
        print "\nDefine Base B abbreviation"
        print " (Use * at end of residue to desinate a wild-card base)"
        print " (Use ** only to search for any Base B)"
        bB = str(raw_input(" Base B: ")).upper()
        break
      except ValueError: print "Invalid base type. Try again." 

    while True:
      print "\nDefine Base B - Chi Angle Lower Bound"
      print "  (Vertical bars around the chi define an absolute. Ex |90|)"
        #Define lower bounds of base B
      lo2 = raw_input(" Chi B lower bound (degrees): ")
      if lo2[:1] == "|" and lo2[-1:] == "|":
        try:
          lo2 = math.fabs(float(lo2[1:-1]))
          break
        except ValueError: print "Invalid input. Try again."   
      elif lo2[:1] != "|" and lo2[-1:] != "|":
        try:
          lo2 = float(lo2)
          break
        except ValueError: print "Invalid input. Try again."   
      else: print "Invalid input. Try again."

    while True:
      print "\nDefine Base B - Chi Angle Upper Bound"
      print "  (Vertical bars around the chi define an absolute. Ex |180|)"
        #Define upper bounds of base B
      hi2 = raw_input(" Chi B upper bound (degrees): ")
      if hi2[:1] == "|" and hi2[-1:] == "|":
        try:
          hi2 = math.fabs(float(hi2[1:-1]))
          break
        except ValueError: print "Invalid input. Try again."   
      elif hi2[:1] != "|" and hi2[-1:] != "|":
        try:
          hi2 = float(hi2)
          break
        except ValueError: print "Invalid input. Try again."   
      else: print "Invalid input. Try again."  

    lo1, hi1 = mainC.sortLoHi(lo1, hi1, "FLOAT")
    lo2, hi2 = mainC.sortLoHi(lo2, hi2, "FLOAT")

      #If no wildcards are used
    if bA[-1:] != "*" and bB[-1:] != "*":
      print "If A and B"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA == x.bAsubtype and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (bB == x.bBsubtype and (lo2 <= float(x.bbChiB) <= hi2))
        )
       or
        #2nd statement
        ( 
           (bA == x.bBsubtype and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (bB == x.bAsubtype and (lo2 <= float(x.bbChiA) <= hi2))
        )      
      ]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

      #if A* and B*
    elif bA[-1:] == "*" and bB[-1:] == "*" and bB != "**":
      print "If A* and B*"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA[:-1] == x.bAtype.upper() and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (bB[:-1] == x.bBtype.upper() and (lo2 <= float(x.bbChiB) <= hi2))
        )
       or
        #2nd statement
        ( 
           (bA[:-1] == x.bBtype.upper() and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (bB[:-1] == x.bAtype.upper() and (lo2 <= float(x.bbChiA) <= hi2))
        )      
      ]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

      #If A* and B
    elif bA[-1:] == "*" and bB[-1:] != "*":
      print "If A* and B"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA[:-1] == x.bAtype.upper() and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (bB == x.bBsubtype.upper() and (lo2 <= float(x.bbChiB) <= hi2))
        )
       or
        #2nd statement
        ( 
           (bA[:-1] == x.bBtype.upper() and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (bB == x.bAsubtype.upper() and (lo2 <= float(x.bbChiA) <= hi2))
        )      
      ]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

      #If A and B*
    elif bA[-1:] != "*" and bB[-1:] == "*" and bB != "**":
      print "If A and B*"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA == x.bAsubtype.upper() and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (bB[:-1] == x.bBtype.upper() and (lo2 <= float(x.bbChiB) <= hi2))
        )
       or
        #2nd statement
        ( 
           (bA == x.bBsubtype.upper() and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (bB[:-1] == x.bAtype.upper() and (lo2 <= float(x.bbChiA) <= hi2))
        )      
      ]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

      #If A and **
    elif bA[-1:] != "*" and bB == "**":
      print "If A and **"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA == x.bAsubtype and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (lo2 <= float(x.bbChiB) <= hi2)
        )
       or
        #2nd statement
        ( 
           (bA == x.bBsubtype and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (lo2 <= float(x.bbChiA) <= hi2)
        )      
      ]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

      #if A* and **
    elif bA[-1:] == "*" and bB == "**":
      print "If A* and **"
      [self.refPairs2.add(x.uniqueID) for x in tempList if
        #1st statement
        ( 
           (bA[:-1] == x.bAtype.upper() and (lo1 <= float(x.bbChiA) <= hi1))
           and
           (lo2 <= math.fabs(float(x.bbChiB)) <= hi2)
        )
       or
        #2nd statement
        ( 
           (bA[:-1] == x.bBsubtype.upper() and (lo1 <= float(x.bbChiB) <= hi1))
           and
           (lo2 <= float(x.bbChiA) <= hi2)
        )      
      ]

      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

    else: print "Invalid combination used."

  def mC1pC1p(self, mainC):
    c1lo, c1hi, tc1p = 0.0, 0.0, 0.0
    while True:
      try:
        print "\nEnter a chemical shift range"
        c1lo = float(raw_input(" C1' lower bound (angstroms): "))
        break
      except ValueError: print "Invalid distance"   
    while True:
      try:
        c1hi = float(raw_input(" C1' upper bound (angstroms): "))
        break
      except ValueError: print "Invalid distance"       

    if c1lo > c1hi:
      tc1p = c1lo
      c1lo = c1hi
      c1hi = tc1p

    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}

    for i in tempList:
      if i.c1c1 != "---":
        if (float(i.c1c1) >= c1lo and float(i.c1c1) <= c1hi):
          self.refPairs2.add(i.uniqueID)

    # self.refPairs2 = set([])
    # [self.refPairs2.add(x.uniqueID) for x in tempList if (( float(x.c1c1) >= c1lo ) and ( float(x.c1c1) <= c1hi))]
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

    ########### Base Pair Presets Menu Start ###########
  def matchBPTypes(self, mainC):
    opt1 = 0

    while True:
      try:
        print "\nSearch for base pair type."
        print " You are currently searching", len(self.refPairs1), "of", len(self.allPairs), "base pairs."
        print " 1. Canonical Pairs"
        print " 2. Non-Canonical Pairs"
        print " 3. Protonated"
        print " 4. Deprotonated"
        print " 5. Rare-Tautomeric"
        # print " 6. Hoogsteen"
        # print " 7. Reverse Hoogsteen"
        print " 9. Undo"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        if opt1 == 9:
          self.refPairs1 = self.refPairs3
        
        if opt1 == 1:
          self.refPairs3 = self.refPairs1
          mainC.findCanonical(True)
          print ""

        if opt1 == 2:
          self.refPairs3 = self.refPairs1
          mainC.findCanonical(False)
          print ""

        if opt1 == 6:
          None #Hoogsteen

        if opt1 == 7:
          None #Reverse Hoogsteen

          #Find protonated
        if opt1 == 3: 
          self.refPairs3 = self.refPairs1
          opt2 = 0
          print ""
          while True:
            try:
              print " 1. Strict search criteria"
              print " 2. Loose search criteria"
              opt2 = int(raw_input(" > "))
              if opt2 == 1:
                A_a = ["N1", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                A_i = ["N1", "N7"]
                G_a = ["O6", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                G_i = ["N7"]
                U_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                U_i = [""]
                T_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                T_i = [""]
                C_a = ["O2", "N3", "O3'", "O4'", "O1P", "O2P"]
                C_i = ["N3"]
                print ""
                mainC.findIonized(2, 3.4, A_a, A_i, G_a, G_i, U_a, U_i, T_a, T_i, C_a, C_i)
                break
                #Loose search criterion
              elif opt2 == 2:
                A_a = ["N1", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                A_i = ["N1", "N7"]
                G_a = ["O6", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                G_i = ["N7"]
                U_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                U_i = [""]
                T_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                T_i = [""]
                C_a = ["O2", "N3", "O3'", "O4'", "O1P", "O2P"]
                C_i = ["N3"]
                # A_a = ["N1", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                # A_i = ["N1", "N7"]
                # G_a = ["O6", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                # G_i = ["N7", "O6"]
                # U_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                # U_i = ["O4"]
                # T_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                # T_i = ["O4"]
                # C_a = ["O2", "N3", "O3'", "O4'", "O1P", "O2P"]
                # C_i = ["N3"]
                print ""
                mainC.findIonized(1, 3.4, A_a, A_i, G_a, G_i, U_a, U_i, T_a, T_i, C_a, C_i)
                break
            except ValueError: print "Invalid menu option"    
        if opt1 == 4: 
          self.refPairs3 = self.refPairs1
          opt2 = 0
          print ""
          A_a = ["N6"]
          A_i = []
          G_a = ["N1", "N2"]
          G_i = ["N1"]
          U_a = ["N3"]
          U_i = ["N3"]
          T_a = ["N3"]
          T_i = ["N3"]
          C_a = ["N4"]
          C_i = []
          while True:
            try:
              print " 1. Strict search criteria"
              print " 2. Loose search criteria"
              opt2 = int(raw_input(" > "))
              if opt2 == 1:
                print ""
                mainC.findIonized(2, 3.4, A_a, A_i, G_a, G_i, U_a, U_i, T_a, T_i, C_a, C_i)
                break
              elif opt2 == 2:
                print ""
                mainC.findIonized(1, 3.4, A_a, A_i, G_a, G_i, U_a, U_i, T_a, T_i, C_a, C_i)
                break
            except ValueError: print "Invalid menu option"    
          #Find tautomers
        if opt1 == 5: 
          opt2 = 0
          print ""
          self.refPairs3 = self.refPairs1

          while True:
            try:
              print " 1. Strict search criteria"
              print " 2. Loose search criteria"
              opt2 = int(raw_input(" > "))
              if opt2 == 1:
                print ""
                A_a = ["N1", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                G_a = ["O6", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                U_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                T_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                C_a = ["O2", "N3", "O3'", "O4'", "O1P", "O2P"]
                mainC.findIonized(2, 3.4, A_a, A_a, G_a, G_a, U_a, U_a, T_a, T_a, C_a, C_a)
                A_a = ["N6"]
                G_a = ["N1", "N2"]
                U_a = ["N3"]
                T_a = ["N3"]
                C_a = ["N4"]
                mainC.findIonized(2, 3.4, A_a, A_a, G_a, G_a, U_a, U_a, T_a, T_a, C_a, C_a)
                break
              elif opt2 == 2:
                print ""
                A_a = ["N1", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                G_a = ["O6", "O3'", "N3", "N7", "O4'", "O1P", "O2P"]
                U_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                T_a = ["O2", "O4", "O3'", "O4'", "O1P", "O2P"]
                C_a = ["O2", "N3", "O3'", "O4'", "O1P", "O2P"]
                mainC.findIonized(1, 3.6, A_a, A_a, G_a, G_a, U_a, U_a, T_a, T_a, C_a, C_a)
                A_a = ["N6"]
                G_a = ["N1", "N2"]
                U_a = ["N3"]
                T_a = ["N3"]
                C_a = ["N4"]
                mainC.findIonized(1, 3.6, A_a, A_a, G_a, G_a, U_a, U_a, T_a, T_a, C_a, C_a)
                break
            except ValueError: print "Invalid menu option"    
      except ValueError: print "Invalid menu option"    
    ########### Base Pair Presets Menu End ###########

    #Takes two input values: float or integer, depending on specified 'type'
    #Returns correct lo / hi values in specified 'type'
    #Returns: LoValue, HiValue
  def sortLoHi(self, inLo, inHi, typeD):
    if typeD == "INT":
      if (inLo <= inHi):
        return(int(inLo), int(inHi))
      else:
        return(int(inHi), int(inLo))

    elif typeD == "FLOAT":
      if (inLo <= inHi):
        return(float(inLo), float(inHi))
      else:
        return(float(inHi), float(inLo))

  def strucWtSearch(self, mainC):
    loWt, hiWt = 0.0, 0.0
    while True:
      try:
        loWt = float(raw_input(" Min Weight (Daltons): "))
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        hiWt = float(raw_input(" Max Weight (Daltons): "))
        break
      except ValueError: print "Invalid number"

      #Ensure min to max values
    if loWt > hiWt:
      tpW = loWt
      loWt = hiWt
      hiWt = tpW

    self.searchWtPDB = self.upClass.matchWeight(loWt, hiWt)
    mainC.pruneListPDB(self.searchWtPDB, "AND")

  def refRibosome(self, mainC):
    opt1 = 0
    while True:
      try:
        print "\nInclude or Exclude Large Ribosomal Subunit Structures"
        print " 1. Include only large ribosomal subunits"
        print " 2. Exclude all large ribosomal subunits"
        print " 0. Return"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        elif opt1 == 1:
          self.ribosomePDB = self.upClass.pruneRibosome()
          mainC.pruneListPDB(self.ribosomePDB, "AND") 
          break
        elif opt1 == 2:
          self.ribosomePDB = self.upClass.pruneRibosome()
          mainC.pruneListPDB(self.ribosomePDB, "NOT") 
          break    
      except ValueError: print "Please select a valid menu option.\n"    

  def saveSearch(self, mainC):
    opt1 = 0
    while True:
      try:
        print "\nSave current search results to Slot"
        print " 1. Slot 1"
        print " 2. Slot 2"
        print " 3. Slot 3"
        print " 4. Slot 4"
        print " 5. Slot 5"
        print " 6. Slot 6"
        print " 7. Slot 7"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break     
        if opt1 == 1:
          FILE = open(self.datDir + "Save1.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close

        if opt1 == 2:
          FILE = open(self.datDir + "Save2.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close     
        if opt1 == 3:
          FILE = open(self.datDir + "Save3.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close    
        if opt1 == 4:
          FILE = open(self.datDir + "Save4.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close                
        if opt1 == 5:
          FILE = open(self.datDir + "Save5.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close     
        if opt1 == 6:
          FILE = open(self.datDir + "Save6.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close    
        if opt1 == 7:
          FILE = open(self.datDir + "Save7.dat", "wb")
          pickle.dump(self.refPairs1, FILE, protocol=2)
          FILE.close    
      except ValueError: print "Please select a valid menu option\n"        

  def loadSearch(self, mainC):
    opt1 = 0
    while True:
      try:
        print "Load results from a previous save"
        print " 1. Slot 1"
        print " 2. Slot 2"
        print " 3. Slot 3"
        print " 4. Slot 4"
        print " 5. Slot 5"
        print " 6. Slot 6"
        print " 7. Slot 7"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break     
        if opt1 == 1:
          if mainC.fileExists(self.datDir + "Save1.dat") == True:
            FILE = open(self.datDir + "Save1.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave1.dat does not exist."

        if opt1 == 2:
          if mainC.fileExists(self.datDir + "Save2.dat") == True:
            FILE = open(self.datDir + "Save2.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave2.dat does not exist."

        if opt1 == 3:
          if mainC.fileExists(self.datDir + "Save3.dat") == True:
            FILE = open(self.datDir + "Save3.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave3.dat does not exist."

        if opt1 == 4:
          if mainC.fileExists(self.datDir + "Save4.dat") == True:
            FILE = open(self.datDir + "Save4.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave4.dat does not exist."

        if opt1 == 5:
          if mainC.fileExists(self.datDir + "Save5.dat") == True:
            FILE = open(self.datDir + "Save5.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave5.dat does not exist."

        if opt1 == 6:
          if mainC.fileExists(self.datDir + "Save6.dat") == True:
            FILE = open(self.datDir + "Save6.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave6.dat does not exist."

        if opt1 == 7:
          if mainC.fileExists(self.datDir + "Save7.dat") == True:
            FILE = open(self.datDir + "Save7.dat", "rb")
            self.refPairs1 = pickle.load(FILE)
            FILE.close
          else: print "\nSave7.dat does not exist."

      except ValueError: print "Please select a valid menu option\n"        

  def findCanonical(self, wcOrnot):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
      #If true, return all canonical WC pairs
    if wcOrnot == True:
      [self.refPairs2.add(x.uniqueID) for x in tempList if x.isWC == wcOrnot]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []
    else:
      [self.refPairs2.add(x.uniqueID) for x in tempList if x.isWC == wcOrnot]
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
      tempList = []

  def findIonized(self, numB, bDist, A_a, A_i, G_a, G_i, U_a, U_i, T_a, T_i, C_a, C_i):   
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1} 
    for i in tempList:
      if int(i.numBonds) >= numB and int(i.numBonds) <= 4 and math.fabs(i.stagger) < 2.0:
        if i.pairMotif.upper() == "AA":
          for y in range(i.numBonds):
            if (((i.bondAtLHS[y] in A_i and i.bondAtRHS[y] in A_a) 
                   or (i.bondAtLHS[y] in A_a and i.bondAtRHS[y] in A_i)) 
                and i.bondDists[y] <= bDist):
              self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "AC":
          if i.bAtype.upper() + i.bBtype.upper() == "AC":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in A_i and i.bondAtRHS[y] in C_a 
                  or i.bondAtLHS[y] in A_a and i.bondAtRHS[y] in C_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in C_i and i.bondAtRHS[y] in A_a 
                  or i.bondAtLHS[y] in C_a and i.bondAtRHS[y] in A_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "AG":
          if i.bAtype.upper() + i.bBtype.upper() == "AG":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in A_i and i.bondAtRHS[y] in G_a 
                  or i.bondAtLHS[y] in A_a and i.bondAtRHS[y] in G_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in G_i and i.bondAtRHS[y] in A_a 
                  or i.bondAtLHS[y] in G_a and i.bondAtRHS[y] in A_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "AT":
          if i.bAtype.upper() + i.bBtype.upper() == "AT":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in A_i and i.bondAtRHS[y] in T_a 
                  or i.bondAtLHS[y] in A_a and i.bondAtRHS[y] in T_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in T_i and i.bondAtRHS[y] in A_a 
                  or i.bondAtLHS[y] in T_a and i.bondAtRHS[y] in A_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "AU":
          if i.bAtype.upper() + i.bBtype.upper() == "AU":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in A_i and i.bondAtRHS[y] in U_a 
                  or i.bondAtLHS[y] in A_a and i.bondAtRHS[y] in U_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in U_i and i.bondAtRHS[y] in A_a 
                  or i.bondAtLHS[y] in U_a and i.bondAtRHS[y] in A_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "CC":
          if i.bAtype.upper() + i.bBtype.upper() == "CC":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in C_i and i.bondAtRHS[y] in C_a 
                  or i.bondAtLHS[y] in C_a and i.bondAtRHS[y] in C_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "CG":
          if i.bAtype.upper() + i.bBtype.upper() == "CG":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in C_i and i.bondAtRHS[y] in G_a 
                  or i.bondAtLHS[y] in C_a and i.bondAtRHS[y] in G_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in G_i and i.bondAtRHS[y] in C_a 
                  or i.bondAtLHS[y] in G_a and i.bondAtRHS[y] in C_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "CT":
          if i.bAtype.upper() + i.bBtype.upper() == "CT":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in C_i and i.bondAtRHS[y] in T_a 
                  or i.bondAtLHS[y] in C_a and i.bondAtRHS[y] in T_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in T_i and i.bondAtRHS[y] in C_a 
                  or i.bondAtLHS[y] in T_a and i.bondAtRHS[y] in C_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "CU":
          if i.bAtype.upper() + i.bBtype.upper() == "CU":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in C_i and i.bondAtRHS[y] in U_a 
                  or i.bondAtLHS[y] in C_a and i.bondAtRHS[y] in U_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in U_i and i.bondAtRHS[y] in C_a 
                  or i.bondAtLHS[y] in U_a and i.bondAtRHS[y] in C_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "GG":
          if i.bAtype.upper() + i.bBtype.upper() == "GG":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in G_i and i.bondAtRHS[y] in G_a 
                  or i.bondAtLHS[y] in G_a and i.bondAtRHS[y] in G_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "GT":
          if i.bAtype.upper() + i.bBtype.upper() == "GT":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in G_i and i.bondAtRHS[y] in T_a 
                  or i.bondAtLHS[y] in G_a and i.bondAtRHS[y] in T_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in T_i and i.bondAtRHS[y] in G_a 
                  or i.bondAtLHS[y] in T_a and i.bondAtRHS[y] in G_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "GU":
          if i.bAtype.upper() + i.bBtype.upper() == "GU":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in G_i and i.bondAtRHS[y] in U_a 
                  or i.bondAtLHS[y] in G_a and i.bondAtRHS[y] in U_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in U_i and i.bondAtRHS[y] in G_a 
                  or i.bondAtLHS[y] in U_a and i.bondAtRHS[y] in G_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "TT":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in T_i and i.bondAtRHS[y] in T_a 
                  or i.bondAtLHS[y] in T_a and i.bondAtRHS[y] in T_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "TU":
          if i.bAtype.upper() + i.bBtype.upper() == "TU":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in T_i and i.bondAtRHS[y] in U_a 
                  or i.bondAtLHS[y] in T_a and i.bondAtRHS[y] in U_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)
          else:
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in U_i and i.bondAtRHS[y] in T_a 
                  or i.bondAtLHS[y] in U_a and i.bondAtRHS[y] in T_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

        elif i.pairMotif.upper() == "UU":
            for y in range(i.numBonds):
              if ((i.bondAtLHS[y] in U_i and i.bondAtRHS[y] in U_a 
                  or i.bondAtLHS[y] in U_a and i.bondAtRHS[y] in U_i) 
                and i.bondDists[y] <= bDist):
                self.refPairs2.add(i.uniqueID)

    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []     

  def matchShifts(self, mainC):
    baseType = ""
    atomStr = ""
    shiftLo, shiftHi = 0.0, 0.0
    tempShifts = []

    while True:
      try:
        print "\nEnter an atom type (ex. C8)"
        atomStr = str(raw_input(" Atom: "))
        atomStr = atomStr.upper()
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        print "\nEnter a chemical shift range"
        shiftLo = float(raw_input(" Shift lower bound (ppm): "))
        break
      except ValueError: print "Invalid shift"   
    while True:
      try:
        shiftHi = float(raw_input(" Shift upper bound (ppm): "))
        break
      except ValueError: print "Invalid shift"   
    while True:
      try:
        print "\nEnter a base pair type to match to shift (ex. G)"
        print " (Input '*' to skip matching base pair type)"
        baseType = str(raw_input(" Base type: "))
        atomStr = atomStr.upper()
        break
      except ValueError: print "Invalid number"

    shiftLo, shiftHi = mainC.sortLoHi(shiftLo, shiftHi, "FLOAT")

            #0 PDB Name
            #1 Pair Motif Sorted Name
            #2 Base A
            #3 Base B
            #4 Temp
            #5 pH
            #6 Atom type
            #7 Chemical Shift
            #8 Shear, 9 stretch
            #10 Stagger, 11 buckle
            #12 Propeller, 13 opening
            #14 Sugar pucker
            #15 Chi Torsion
    
    if baseType == "*":
      tempShifts = [x for x in self.shiftData if x[6] == atomStr and (shiftLo <= x[7] <= shiftHi)]
    else:
      tempShifts = [x for x in self.shiftData if x[6] == atomStr and (shiftLo <= x[7] <= shiftHi) 
                              and (baseType in x[2])]

    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    ([self.refPairs2.add(x.uniqueID) for x in tempList for y in tempShifts 
      if (x.bbChiAconf == y[15] or x.bbChiBconf == y[15]) and
      (x.puckerA == y[14] or x.puckerB == y[14]) and
      (mainC.exceedsZ6(1.5, x.shear, x.stretch, x.stagger,
                                            x.buckle, x.propeller, x.opening,
                                            y[8], y[9], y[10],
                                            y[11], y[12], y[13],
                                            1.1,1.75,1.75,5,5,5) == True)])
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []


  def readShiftDB(self, menuC, fileStr):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    motifs = set([])
    counter = 1
    cnt = 0
    tLo, tHi, phLo, phHi = 0, 0, 0.0, 0.0

    while True:
      try:
        print "\nDefine a temperature range for the shifts"
        tLo = int(raw_input(" Low temp(K): "))
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        tHi = int(raw_input(" High temp(K): "))
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        phLo = float(raw_input(" Low pH: "))
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        phHi = float(raw_input(" High pH: "))
        break
      except ValueError: print "Invalid number"

      #Ensure min to max values
    tLo, tHi = menuC.sortLoHi(tLo, tHi, "FLOAT")
    phLo, phHi = menuC.sortLoHi(phLo, phHi, "FLOAT")

    infile = open(fileStr, "rU")
    reader = csv.reader(infile, delimiter=",")
    for i in reader:
      if counter != 1:
        if (i[10] != "N/A" and i[35] != "N/A" and i[69] != "N/A"
            and (tLo <= int(i[18]) <= tHi)
            and (phLo <= float(i[19]) <= phHi)):
          self.shiftData.append([])
            #0 PDB Name
          self.shiftData[cnt].append(i[31])
            #1 Pair Motif Sorted Name
          if i[10] < i[5]:
            self.shiftData[cnt].append(i[10] + i[5])
          else:
            self.shiftData[cnt].append(i[5] + i[10])
            #2 Base A
          self.shiftData[cnt].append(i[5])
            #3 Base B
          self.shiftData[cnt].append(i[10])
            #4 Temp
          self.shiftData[cnt].append(int(i[18]))
            #5 pH
          self.shiftData[cnt].append(float(i[19]))
            #6 Atom type
          self.shiftData[cnt].append(i[4])
            #7 Chemical Shift
          self.shiftData[cnt].append(float(i[1]))
            #8 Shear, 9 stretch
          self.shiftData[cnt].append(float(i[35])), self.shiftData[cnt].append(float(i[34]))
            #10 Stagger, 11 buckle
          self.shiftData[cnt].append(float(i[33])), self.shiftData[cnt].append(float(i[36]))
            #12 Propeller, 13 opening
          self.shiftData[cnt].append(float(i[37])), self.shiftData[cnt].append(float(i[38]))
            #14 Sugar pucker
          self.shiftData[cnt].append(i[32])
            #15 Chi Torsion
          if (math.fabs(float(i[69])) >= 90):
            self.shiftData[cnt].append("Anti")
          else: self.shiftData[cnt].append("Syn")
          cnt += 1
      counter += 1
    infile.close

  def checkShifts(self, inT, loT, hiT, inph, lopH, hipH, 
                             v1, v2, v3, v4, v5, v6, m1, m2, m3, m4, m5, m6):
    if ((self.exceedsZ6(1.25, v1, v2, v3, v4, v5, v6, m1, m2, m3, m4, m5, m6, 2, 2, 2, 5, 5, 5) == True)
        and ((loT <= inT <= hiT) == True) and ((lopH <= inph <= hipH) == True)):
      return(True)
    else: return(False)

  def mChiTors(self, mainC):
    chiAS, chiBS = 0.0, 0.0
    print "Input a base pair PDB ID to search for"
    self.singleC = mainC.findPairValsPDB()
    print self.singleC[0].uniqueID   
    chiA = float(self.singleC[0].bbChiA)
    chiB = float(self.singleC[0].bbChiB)
    while True:
      try:   
        print ("Input standard deviation to match ChiA: %s(%s%s, %s)" 
                   % (self.singleC[0].bAstrand, self.singleC[0].bAtype,
                       str(self.singleC[0].bAnum), str(chiA)))
        chiAS = float(raw_input(" > "))
        break
      except ValueError: print "Please use a real number\n"

    while True:
      try:   
        print ("Input standard deviation to match ChiV: %s(%s%s, %s)" 
                   % (self.singleC[0].bBstrand, self.singleC[0].bBtype,
                       str(self.singleC[0].bBnum), str(chiB)))
        chiBS = float(raw_input(" > "))
        break
      except ValueError: print "Please use a real number\n"

    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1 and 
                        (self.allPairs[x].bbChiA != "---" and self.allPairs[x].bbChiB != "---")}
    ([self.refPairs2.add(x.uniqueID) for x in tempList if 
      ((mainC.exceedsZ2(1, float(x.bbChiA), float(x.bbChiB),
                                       chiA, chiB, chiAS, chiBS) == True) or
      (mainC.exceedsZ2(1, float(x.bbChiB), float(x.bbChiA),
                                       chiA, chiB, chiAS, chiBS) == True))])
  # def exceedsZ2(self, zThresh, inV1, inV2, mV1, mV2, vS1, vS2):
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

  def mPucker(self, mainC):
    print "Input a base pair PDB ID to search for"
    self.singleC = mainC.findPairValsPDB()
    print self.singleC[0].uniqueID
    pA, pB = "", ""
    pA = self.singleC[0].puckerA
    pB = self.singleC[0].puckerB
    print "", pA, pB
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    ([self.refPairs2.add(x.uniqueID) for x in tempList if 
      ((x.puckerA == pA and x.puckerB == pB) or (x.puckerA == pB and x.puckerB == pA))])
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

    #Can search by range of any of the local pair geometry variables (ie shear)
  def findLocalGeom(self, mainC):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    sheV, strV, staV, bucV, proV, opeV = 0., 0., 0., 0., 0., 0.
    opt1, opt2 = None, None
    vLo, vHi = None, None
    [self.refPairs2.add(x.uniqueID) for x in tempList]

    while True:
      try:
        print "\nWhich variable to do wish to refine by?"
        print (" You are currently searching %d of %d base pairs" % 
              (len(self.refPairs2), len(self.refPairs1)))
        print " 1)  Shear"
        print " 2)  Stretch"
        print " 3)  Stagger"
        print " 4)  Buckle"
        print " 5)  Propeller"
        print " 6)  Opening"
        print " 0)  Return to previous menu."
        opt1 = int(raw_input(" > "))

    #Takes two input values: float or integer, depending on specified 'type'
    #Returns correct lo / hi values in specified 'type'
    #Returns: LoValue, HiValue
  # def sortLoHi(self, inLo, inHi, typeD):

        if opt1 == 0:
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
          tempList = []  
          break
        elif 1 <= opt1 <= 6:
          while True:
            try:
              print "\n  1) Include range"
              print "  2) Exclude range"
              opt2 = int(raw_input(" > "))
              if 1 <= opt2 <= 2:
                break
              else: print "Invalid choice, try again."
            except ValueError: print "Please provide a valid menu option\n"
            #Get shear lower value
          while True:
            try:
              vLo = float(raw_input(" Lower Bounds > "))
              break
            except ValueError: print "Please provide a valid number\n"
            #Get shear upper value
          while True:
            try:
              vHi = float(raw_input(" Upper Bounds > "))
              break
            except ValueError: print "Please provide a valid number\n"
          vLo, vHi = mainC.sortLoHi(vLo, vHi, "FLOAT")
            #Include these
          if opt2 == 1:
            if opt1 == 1:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.shear <= vHi]
            elif opt1 == 2:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.stretch <= vHi]
            elif opt1 == 3:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.stagger <= vHi]
            elif opt1 == 4:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.buckle <= vHi]
            elif opt1 == 5:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.propeller <= vHi]
            elif opt1 == 6:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if not vLo <= x.opening <= vHi]
            #Exclude these
          elif opt2 == 2:
            if opt1 == 1:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.shear <= vHi]
            elif opt1 == 2:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.stretch <= vHi]
            elif opt1 == 3:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.stagger <= vHi]
            elif opt1 == 4:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.buckle <= vHi]
            elif opt1 == 5:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.propeller <= vHi]
            elif opt1 == 6:
              [self.refPairs2.discard(x.uniqueID) for x in tempList 
               if vLo <= x.opening <= vHi]

      except ValueError: print "Please input a integer\n"

  def mLocPairPars(self, mainC):
    opt1, fVal = "", 0.0
    shearS, stretchS, staggerS, buckleS, propellerS, openingS = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    print "Input a base pair PDB ID to search for"
    self.singleC = mainC.findPairValsPDB()
    print self.singleC[0].uniqueID
    while True:
      try:
        print "\nRefine search parameters (define standard deviations)"
        print " (Use format 'She1.5' to set a shear standard deviation to 1.5, etc)"
        print (" [SHE] %s - [STR] %s - [STA] %s - [BUC] %s - [PRO] %s - [OPE] %s" 
                    % (str(shearS), str(stretchS), str(staggerS),
                       str(buckleS), str(propellerS), str(openingS)))
        print " Type 'Q' to complete changes"
        opt1 = str(raw_input(" > "))
        if (((opt1[:3].upper() == "SHE") or (opt1[:3].upper() == "STR")
            or (opt1[:3].upper() == "STA") or (opt1[:3].upper() == "BUC") 
            or (opt1[:3].upper() == "PRO") or  (opt1[:3].upper() == "OPE")
            or (opt1[:3].upper() == "ALL"))
            and (type(float(opt1[3:]) == type(fVal)))):
          if opt1[:3].upper() == "SHE":
            shearS = float(opt1[3:])
          if opt1[:3].upper() == "STR":
            stretchS = float(opt1[3:])
          if opt1[:3].upper() == "STA":
            staggerS = float(opt1[3:])
          if opt1[:3].upper() == "BUC":
            buckleS = float(opt1[3:])    
          if opt1[:3].upper() == "PRO":
            propellerS = float(opt1[3:])
          if opt1[:3].upper() == "OPE":
            openingS = float(opt1[3:])
          if opt1[:3].upper() == "ALL":
            shearS = float(opt1[3:])        
            stretchS = float(opt1[3:])
            staggerS = float(opt1[3:])
            buckleS = float(opt1[3:])
            propellerS = float(opt1[3:])
            openingS = float(opt1[3:])         
        elif (opt1.upper() == "Q"):
          mainC.rangePairLoc(mainC, shearS, stretchS, staggerS, buckleS, propellerS, openingS,
                                              self.singleC[0].shear, self.singleC[0].stretch, 
                                              self.singleC[0].stagger, self.singleC[0].buckle, 
                                              self.singleC[0].propeller, self.singleC[0].opening)
          break
        else: 
          print "\nPlease input a string that begins with either:"
          print " [She]ar"
          print " [Stre]tch"
          print " [Sta]gger"
          print " [Buc]kle"
          print " [Pro]peller"
          print " [Ope]ning"
      except ValueError: print "Please input a valid string\n"

  def readPDBlist(self):
    pdbListName = None
    listPath = None
    tPDBlist = []
    tPDBset = set([])

    while True:
      try:
        print "\nInput name of local .txt file containing list of PDB names."
        pdbListName = str(raw_input(" > "))

        if self.fileExists(os.path.join(self.curDir, pdbListName)) == True:
          FILE = open(os.path.join(self.curDir, pdbListName), "rU")
          tPDBlist = [x.strip() for x in FILE]
          FILE.close()

          for xname in tPDBlist:
            if xname[-4:].lower() == ".pdb":
              tPDBset.add(xname[:-4].upper())
            else:
              tPDBset.add(xname.upper())
          {self.refPairs2.add(j) for j in self.refPairs1 if j[:4] in tPDBset}
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
          break
        else:
          print "Invalid file. Try again."
          pdbListName = None
      except ValueError: print "Incorrect name\n" 

  def mRetPDBPairs(self):
    listNames = []
    tempName = ""
    while True:
      try:
        print "\nEnter one PDB name per line. Type 'Quit to stop entering names."
        tempName = str(raw_input(" > "))
        tempName = tempName.upper()
        if len(tempName) > 4:
          print "PDB name", tempName, "is too long. Try again."
        elif len(tempName) < 4:
          print "PDB name", tempName, "is too short. Try again."
        elif len(tempName) == 4 and tempName != "QUIT":
          listNames.append(tempName)
        else:
          print ""
          break

      except ValueError: print "Incorrect name\n"      
    {self.refPairs2.add(j) for j in self.refPairs1 if j[:4] in listNames}
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
  
  def rangePairLoc(self, mainC, sheS, strS, staS, bucS, proS, opeS,
                                sheV, strV, staV, bucV, proV, opeV):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}

    ([self.refPairs2.add(x.uniqueID) for x in tempList 
      if mainC.exceedsZ6(1.5, x.shear, x.stretch, x.stagger, x.buckle,
                                                  x.propeller, x.opening, sheV, strV, staV,
                                                  bucV, proV, opeV, sheS, strS, staS, bucS,
                                                  proS, opeS) == True])  

    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

    #Returns a true/false boolean based on if the calc'ed Z-Score
    #exceeds given Z-score threshold. Comparative values are given
    #by incVal, meanVal is the comparators value, valStdev is the defined
    #standard deviation
  def exceedsZ6(self, zThresh, inV1, inV2, inV3, inV4, inV5, inV6, 
                            mV1, mV2, mV3, mV4, mV5, mV6,
                            vS1, vS2, vS3, vS4, vS5, vS6):
    Z1T, Z2T, Z3T, Z4T, Z5T, Z6T = False, False, False, False, False, False
    calcZ1 = ((inV1 - mV1) / vS1)
    calcZ2 = ((inV2 - mV2) / vS2)
    calcZ3 = ((inV3 - mV3) / vS3)
    calcZ4 = ((inV4 - mV4) / vS4)
    calcZ5 = ((inV5 - mV5) / vS5)
    calcZ6 = ((inV6 - mV6) / vS6)

    if math.fabs(calcZ1) <= zThresh:
      Z1T = True
    if math.fabs(calcZ2) <= zThresh:
      Z2T = True
    if math.fabs(calcZ3) <= zThresh:
      Z3T = True
    if math.fabs(calcZ4) <= zThresh:
      Z4T = True
    if math.fabs(calcZ5) <= zThresh:
      Z5T = True
    if math.fabs(calcZ6) <= zThresh:
      Z6T = True

    if ((Z1T == True) and (Z2T == True) and (Z3T == True) and 
        (Z4T == True) and (Z5T == True) and (Z6T == True)):
      return (True)

  def exceedsZ2(self, zThresh, inV1, inV2, mV1, mV2, vS1, vS2):
    Z1T, Z2T = False, False
    
      #Need to account for the fact that chi rotations -160 to -180 are nearly the same
      # as chi torsions from +160 to +180, and same for 0 to 20.
    if math.fabs(inV1) >= 160 or math.fabs(inV1) <= 20:
      calcZ1 = ((math.fabs(inV1) - math.fabs(mV1)) / vS1)
    else:
      calcZ1 = ((inV1 - mV1) / vS1)

    if math.fabs(inV2) >= 160 or math.fabs(inV2) <= 20:
      calcZ2 = ((math.fabs(inV2) - math.fabs(mV2)) / vS1)
    else:
      calcZ2 = ((inV2 - mV2) / vS2)

    if math.fabs(calcZ1) <= zThresh:
      Z1T = True
    if math.fabs(calcZ2) <= zThresh:
      Z2T = True

    if ((Z1T == True) and (Z2T == True)):
      return (True)

  def findPairValsPDB(self):
    pdbName, strandA, strandB, resiAnum, resiBnum = "", "", "", "", ""
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    while True:
      pdbName = str(raw_input("PDB: ")).upper()
      strandA = str(raw_input("StrandA: "))
      strandB = str(raw_input("StrandB: "))
      resiAnum = str(raw_input("ResiA: "))
      resiBnum = str(raw_input("ResiB: "))

      sPairVals = ([x for x in tempList if (x.pName == pdbName) and (((x.bAstrand == strandA)
                                                                and (x.bBstrand == strandB) 
                                                                and (str(x.bAnum) == resiAnum)
                                                                and (str(x.bBnum) == resiBnum))
                                                                or ((x.bBstrand == strandA)
                                                                and (x.bAstrand == strandB) 
                                                                and (str(x.bBnum) == resiAnum)
                                                                and (str(x.bAnum) == resiBnum)))])
      if len(sPairVals) == 0:
        print "\nNo basepair found for those parameters, try again."
      if len(sPairVals) > 1:
        print "\nFound multiple results for this search. Using first value."
        return(sPairVals)
        tempList = []
      if len(sPairVals) == 1:
        return(sPairVals)
        tempList = []    

  def searchHBondsFine(self):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}

    atom1, atom2, base1, base2, baseCombo, atomStr = "", "", "", "", "", ""
    hbDistLo, hbDistHi = 0.,0.

    while True:
      try:
        print "\nInput Hydrogen Bonding Base #1"
        print "  (IE. G, DA, C, U, etc. Wildcard * matches all base subtypes)"
        base1 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Hydrogen Bonding Atom #1"
        print "  (IE. N1, N7, O2')"
        atom1 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Hydrogen Bonding Base #2"
        print "  (IE. G, DA, C, U, etc. Wildcard * matches all base subtypes)"
        base2 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Hydrogen Bonding Atom #2"
        print "  (IE. N1, N7, O2')"
        atom2 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Lower Limit Hydrogen Bonding Distance (angstroms)"
        hbDistLo = float(raw_input(" > "))
        break
      except ValueError: print "Please input a real number.\n"
    while True:
      try:
        print "\nInput Upper Limit Hydrogen Bonding Distance (angstroms)"
        hbDistHi = float(raw_input(" > "))
        break
      except ValueError: print "Please input a real number.\n"
    hbDistLo, hbDistHi = self.sortLoHi(hbDistLo, hbDistHi, "FLOAT")
    if atom1 < atom2:
      atomStr = atom1 + atom2
    else:
      atomStr = atom2 + atom1

    # Map bases to uniform types
    base_map = {"DG": "G", "RG": "G", "G": "G",
                "DT": "T", "RT": "T", "T": "T",
                "DC": "C", "RC": "C", "C": "C",
                "DA": "A", "RA": "A", "A": "A",
                "DU": "U", "RU": "U", "U": "U"}
    # Define types of bases (base type or subtype)
    # to match against user defined base
    # Default is subtype
    bAt = "bAsubtype"
    bBt = "bBsubtype"

    # Check for base subtypes
    if "*" in base1:
      base1 = base1.replace("*", "").upper()
      # Make sure uniform base is correctly mapped
      if base1 in base_map.keys():
        base1 = base_map[base1]
        bAt = "bAtype"
      else:
        print "Bad base type given."
    # Check base2 types
    if "*" in base2:
      base2 = base2.replace("*", "").upper()
      # Make sure uniform base is correctly mapped
      if base2 in base_map.keys():
        base2 = base_map[base2]
        bBt = "bBtype"
      else:
        print "Bad base type given."

    #List comprehension for finding match for both bases and bonding atoms and distances
    #Have to index to each of 4 possible hbonds
    for x in tempList:
      # Assign the matching base A type or subtype
      if bAt == "bAsubtype":
        bAt = x.bAsubtype
      else:
        bAt = x.bAtype
      # Assign matching base B type/subtype
      if bBt == "bBsubtype":
        bBt = x.bBsubtype
      else:
        bBt = x.bBtype

      # Add set #1
      if ((
          (
            (x.bondAtLHS[0] == atom1 and bAt.upper() == base1) 
            and (x.bondAtRHS[0] == atom2 and bBt.upper() == base2)
          ) 
          or 
          (
            (x.bondAtRHS[0] == atom1 and bBt.upper() == base1) 
            and (x.bondAtLHS[0] == atom2 and bAt.upper() == base2)
          )
         )
         and (hbDistLo <= x.bondDists[0] <= hbDistHi)):
        self.refPairs2.add(x.uniqueID)

      # add set #2
      if ((
          (
            (x.bondAtLHS[1] == atom1 and bAt.upper() == base1) 
            and (x.bondAtRHS[1] == atom2 and bBt.upper() == base2)
          ) 
          or 
          (
            (x.bondAtRHS[1] == atom1 and bBt.upper() == base1) 
            and (x.bondAtLHS[1] == atom2 and bAt.upper() == base2)
          )
         )
         and (hbDistLo <= x.bondDists[1] <= hbDistHi)):
        self.refPairs2.add(x.uniqueID)

      # Add set #3
      if ((
          (
            (x.bondAtLHS[2] == atom1 and bAt.upper() == base1) 
            and (x.bondAtRHS[2] == atom2 and bBt.upper() == base2)
          ) 
          or 
          (
            (x.bondAtRHS[2] == atom1 and bBt.upper() == base1) 
            and (x.bondAtLHS[2] == atom2 and bAt.upper() == base2)
          )
         )
         and (hbDistLo <= x.bondDists[2] <= hbDistHi)):
        self.refPairs2.add(x.uniqueID)

      # Add set #4
      if ((
          (
            (x.bondAtLHS[3] == atom1 and bAt.upper() == base1) 
            and (x.bondAtRHS[3] == atom2 and bBt.upper() == base2)
          ) 
          or 
          (
            (x.bondAtRHS[3] == atom1 and bBt.upper() == base1) 
            and (x.bondAtLHS[3] == atom2 and bAt.upper() == base2)
          )
         )
         and (hbDistLo <= x.bondDists[3] <= hbDistHi)):
        self.refPairs2.add(x.uniqueID)

    # ([self.refPairs2.add(x.uniqueID) for x in tempList 
    #   if (
    #       (
    #         (x.bondAtLHS[0] == atom1 and x.bAsubtype == base1) 
    #         and (x.bondAtRHS[0] == atom2 and x.bBsubtype == base2)
    #       ) 
    #       or 
    #       (
    #         (x.bondAtRHS[0] == atom1 and x.bBsubtype == base1) 
    #         and (x.bondAtLHS[0] == atom2 and x.bAsubtype == base2)
    #       )
    #      )
    #      and (hbDistLo <= x.bondDists[0] <= hbDistHi)
    #   ])

    # ([self.refPairs2.add(x.uniqueID) for x in tempList 
    #   if (
    #       (
    #         (x.bondAtLHS[1] == atom1 and x.bAsubtype == base1) 
    #         and (x.bondAtRHS[1] == atom2 and x.bBsubtype == base2)
    #       ) 
    #       or 
    #       (
    #         (x.bondAtRHS[1] == atom1 and x.bBsubtype == base1) 
    #         and (x.bondAtLHS[1] == atom2 and x.bAsubtype == base2)
    #       )
    #      )
    #      and (hbDistLo <= x.bondDists[1] <= hbDistHi)
    #   ])

    # ([self.refPairs2.add(x.uniqueID) for x in tempList 
    #   if (
    #       (
    #         (x.bondAtLHS[2] == atom1 and x.bAsubtype == base1) 
    #         and (x.bondAtRHS[2] == atom2 and x.bBsubtype == base2)
    #       ) 
    #       or 
    #       (
    #         (x.bondAtRHS[2] == atom1 and x.bBsubtype == base1) 
    #         and (x.bondAtLHS[2] == atom2 and x.bAsubtype == base2)
    #       )
    #      )
    #      and (hbDistLo <= x.bondDists[2] <= hbDistHi)
    #   ])
   
    # ([self.refPairs2.add(x.uniqueID) for x in tempList 
    #   if (
    #       (
    #         (x.bondAtLHS[3] == atom1 and x.bAsubtype == base1) 
    #         and (x.bondAtRHS[3] == atom2 and x.bBsubtype == base2)
    #       ) 
    #       or 
    #       (
    #         (x.bondAtRHS[3] == atom1 and x.bBsubtype == base1) 
    #         and (x.bondAtLHS[3] == atom2 and x.bAsubtype == base2)
    #       )
    #      )
    #      and (hbDistLo <= x.bondDists[3] <= hbDistHi)
    #   ])

    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

  def searchHBonds(self):
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    atom1, atom2, atomStr = "", "", ""
    hbDist = 0.0
    while True:
      try:
        print "\nInput Hydrogen Bonding Atom #1"
        print "  (IE. N1, N7, O2')"
        atom1 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Hydrogen Bonding Atom #2"
        atom2 = str(raw_input(" > "))
        break
      except ValueError: print "Please input a correct string.\n"
    while True:
      try:
        print "\nInput Hydrogen Bonding Distance (angstroms)"
        hbDist = float(raw_input(" > "))
        break
      except ValueError: print "Please input a real number.\n"

    if atom1 < atom2:
      atomStr = atom1 + atom2
    else:
      atomStr = atom2 + atom1
    [self.refPairs2.add(x.uniqueID) for x in tempList if x.bondAtCombo[0] 
    == atomStr and x.bondDists[0] <= hbDist]
    [self.refPairs2.add(x.uniqueID) for x in tempList if x.bondAtCombo[1] 
    == atomStr and x.bondDists[1] <= hbDist]
    [self.refPairs2.add(x.uniqueID) for x in tempList if x.bondAtCombo[2] 
    == atomStr and x.bondDists[2] <= hbDist]
    [self.refPairs2.add(x.uniqueID) for x in tempList if x.bondAtCombo[3] 
    == atomStr and x.bondDists[3] <= hbDist]
    self.refPairs1 = self.refPairs1 & self.refPairs2
    self.refPairs2 = set([])
    tempList = []

  def writeStats(self, mClass):
    mClass.writeOut()
    self.statsC.plotStats((self.writeOutFolder + self.writeOutName), 0)

  def writeClusts(self, mClass):
    numClusts = 0
    while True:
      try:
        print "\nInput maximum number of clusters for convergence"
        print " (warning: >100 clusters can become extremely slow)"
        numClusts = int(raw_input(" > "))
        mClass.writeOut()
        self.statsC.plotStats((self.writeOutFolder + self.writeOutName), 0)       
        if not os.path.exists(self.writeOutFolder + "clusters/"):
          os.makedirs(self.writeOutFolder  + "clusters/")
        bpfPym.writeClusts((self.writeOutFolder + self.writeOutName).strip(".csv"), 
                                          (self.writeOutFolder + "clusters/" + self.writeOutName.strip(".csv")),
                                              (self.writeOutFolder + "clusters/"), 
                                              (self.curDir), numClusts)

        os.system("R --arch x86_64 --silent < tempClust.r --no-save > clustOut.txt")
        if mClass.fileExists("tempClust.r") == True: pass
          # os.system("rm tempClust.r")        
        if mClass.fileExists("clustOut.txt") == True:
          os.system("rm clustOut.txt")    
        break
      except ValueError: print "Please select a valid number.\n"

  def writePymol(self, mClass):
    tempFiletype = "PSE"
    opt1, opt2 = 0, 0
    tempDist = 0.0

    while True:
      try:
        print "\nSelect the following option for creating Pymol files"
        print " 1. [%s] Filetype (PSE or PDB)"  % tempFiletype
        print " 2. [%s] Inclusion radius around base pair (in angstroms)" % str(tempDist)
        print " 8. Quit without saving"
        print " 9. Reset Parameters"
        print " 0. Save and Return"
        opt1 = int(raw_input(" > "))
        if opt1 == 8:
          print ""
          break
        if opt1 == 9:
          tempDist = 0.0
          tempFiletype = "PSE"
        if opt1 == 1:
          while True:
            try:
              print "\nExport as PSE or PDB file format"
              print " 1. PSE"
              print " 2. PDB"
              opt2 = int(raw_input(" > "))
              if opt2 == 1:
                tempFiletype = "PSE"
                print ""
                break
              elif opt2 == 2:
                tempFiletype = "PDB"
                print ""
                break
              else:
                print "Invalid option."
            except ValueError: print "Please select a correct menu option.\n"
        if opt1 == 2:
          while True:
            try:
              print "\nInput an inclusion radius around base pair to add to %s files" % tempFiletype
              tempDist = float(raw_input(" > "))
              tempDist = round(tempDist, 2)
              break
            except ValueError: print "Please input a real number.\n"
        if opt1 == 0:
          mClass.writeOut()
          self.statsC.plotStats((self.writeOutFolder + self.writeOutName), 0)
          if not os.path.exists(self.writeOutFolder + "PDB/"):
            os.makedirs(self.writeOutFolder  + "PDB/")
          bpfPym.writeProgram((self.writeOutFolder + self.writeOutName), 
                                                (self.writeOutFolder + "PDB/"), tempDist, tempFiletype,
                                                (self.realDir + "tempPDB/"), self.curDir)
            #Try to decide whether Linux or OS X is running, and launch pymol accordingly
          platformStr = sys.platform
          if "darwin" in platformStr:
            os.system("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL -qcr tempBpfPym.py")
          elif "linux" in platformStr:
            os.system("/home/usr/share/pymol/pymol -qcr tempBpfPym.py")
          if mClass.fileExists("tempBpfPym.py") == True:
            os.system("rm tempBpfPym.py")
          print ""
          break
      except ValueError: print "Please select a correct menu option.\n"   

  def fileExists(self, pathToFile):
    try: 
      with open(pathToFile) as f: return(True)
    except IOError as e: return(False)    

  def matchMotif(self): 
    str1, str2, combo = "", "", ""
    tempList = {self.allPairs[x] for x in self.allPairs if x in self.refPairs1}
    while True:
      try:
        print "\nEnter 1st base pair type (ie. G, DA, PSU, etc)"
        print "Asterisk (*) after base pair designation denotes a wildcard residue (IE. C* will match DC and C)"
        str1 = str(raw_input(" > "))
        print "\nEnter 2nd base pair type (ie. G, DA, PSU, etc)"
        print " (Enter '**' to search for any base)"
        str2 = str(raw_input(" > "))
        print ""
        
        if ((str1[-1:] == "*" and str1[:1] == "*") or (str2[-1:] == "*" and str2[:1] == "*")):
          if (str1[-1:] == "*" and str1[:1] == "*"):
            combo = str2.upper()

          else:
            combo = str1.upper()

          [self.refPairs2.add(x.uniqueID) for x in tempList if ((x.bAsubtype == combo)
                                                                                               or (x.bBsubtype == combo))]
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])

        elif str1[-1:] == "*" and str2[-1:] == "*":
          if str1.upper() < str2.upper():
            combo = str1[:1] + str2[:1]
          else:
            combo = str2[:1] + str1[:1]
          [self.refPairs2.add(x.uniqueID) for x in tempList if x.pairMotif.upper() == combo]
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
        
        elif str1[-1:] == "*" and str2[-1:] != "*":
          if str1.upper() < str2.upper():
            combo = str1[:1] + str2
          else:
            combo = str2 + str1[:1]
          #if bAsubtype + bBtype | or bAtype + bBsubtype | bBtype + bAsubtype | bBsubtype + bAtype
          ([self.refPairs2.add(x.uniqueID) for x in tempList 
            if (x.bAsubtype + x.bBtype == combo) or
               (x.bAtype + x.bBsubtype == combo) or
               (x.bBsubtype + x.bAtype == combo) or
               (x.bBtype + x.bAsubtype == combo)])
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])

        elif str1[-1:] != "*" and str2[-1:] == "*":
          if str1.upper() < str2.upper():
            combo = str1 + str2[:1]
          else:
            combo = str2[:1] + str1
          ([self.refPairs2.add(x.uniqueID) for x in tempList 
            if (x.bAsubtype + x.bBtype == combo) or
               (x.bAtype + x.bBsubtype == combo) or
               (x.bBsubtype + x.bAtype == combo) or
               (x.bBtype + x.bAsubtype == combo)])
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])

        else:
          if str1.upper() < str2.upper():
            combo = str1 + str2
          else:
            combo = str2 + str1
          [self.refPairs2.add(x.uniqueID) for x in tempList if x.subPMotif == combo]
          self.refPairs1 = self.refPairs1 & self.refPairs2
          self.refPairs2 = set([])
        tempList = []
        break
      except ValueError: print "Please search for a name.\n"   

      #Prune inpPairs (a list of pairs class objects) by PDB names in compNames
  def pruneListPDB(self, compNames, setType):
    if setType == "AND":
      {self.refPairs2.add(j) for j in self.refPairs1 if j[:4] in compNames}
      self.refPairs1 = self.refPairs1 & self.refPairs2
      self.refPairs2 = set([])
    elif setType == "NOT":
      {self.refPairs2.add(j) for j in self.refPairs1 if j[:4] in compNames}
      listnames = set([])
      self.refPairs1 = self.refPairs1 - self.refPairs2
      self.refPairs2 = set([])

    #Resets search list back to original values
  def resetList(self):
    {self.refPairs1.add(x) for x in self.allPairs}
    self.refPairs2 = set([])

  def tempSearchMenu(self, mainC):
    tLo, tHi = 0.0, 0.0
    opt1 = 0
    while True:
      try:
        print "\nSelect to search by crystal growing temp or crystal temp at beamline."
        print " 1. Temperature (Growing Conditions)"
        print " 2. Temperature (At Beamline)"
        print " 0. Return to previous menu"
        opt1 = int(raw_input(" > "))
        if opt1 == 0:
          print ""
          break
        elif opt1 == 1:
          while True:
            try:
              print "\nDefine a temperature range"
              tLo = int(raw_input(" Low temp(K): "))
              break
            except ValueError: print "Invalid number"
          while True:
            try:
              tHi = int(raw_input(" High temp(K): "))
              break
            except ValueError: print "Invalid number"
            #Ensure min to max values
          if tLo > tHi:
            tTemp = tLo
            tLo = tHi
            tHi = tTemp
          self.tempPairs = self.upClass.tempCrystal(tLo, tHi)
          mainC.pruneListPDB(self.tempPairs, "AND")    
          break
        elif opt1 == 2:
          while True:
            try:
              print "\nDefine a temperature range"
              tLo = int(raw_input(" Low temp(K): "))
              break
            except ValueError: print "Invalid number"
          while True:
            try:
              tHi = int(raw_input(" High temp(K): "))
              break
            except ValueError: print "Invalid number"
            #Ensure min to max values
          if tLo > tHi:
            tTemp = tLo
            tLo = tHi
            tHi = tTemp
          self.tempPairs = self.upClass.tempBeamline(tLo, tHi)
          mainC.pruneListPDB(self.tempPairs, "AND")  
          break
      except ValueError: print "Please select a valid menu option\n"      

  def phSearchMenu(self, mainC):
    phLo, phHi = 0.0, 0.0
    while True:
      try:
        phLo = float(raw_input(" Low pH: "))
        break
      except ValueError: print "Invalid number"
    while True:
      try:
        phHi = float(raw_input(" High pH: "))
        break
      except ValueError: print "Invalid number"

      #Ensure min to max values
    phLo, phHi = mainC.sortLoHi(phLo, phHi, "FLOAT")
    self.phPairs = self.upClass.phCrystal(phLo, phHi)
    mainC.pruneListPDB(self.phPairs, "AND")

  def searchSeq(self, mainC):
    str1 = ""
    while True:
      try:
        print "\nSearch for parent PDB containing a user-defined sequence"
        print " (Ex. 'XUUCGX' searches for UUCG with wild-card flanking residues)"
        print " (Regex can be used. Ex. 'CG{4}AT' searches for 'CGXXXXAT')"
        print " ( '-' before string defines NOT including)"
        str1 = str(raw_input(" > "))
        print ""
        if str1[:1] != "-":
          self.seqPairs = self.upClass.seqSearch(str1)
          mainC.pruneListPDB(self.seqPairs, "AND") 
        elif str1[:1] == "-":
          self.seqPairs = self.upClass.seqSearch(str1[1:])
          mainC.pruneListPDB(self.seqPairs, "NOT") 
        break
      except ValueError: print "Please search for a name.\n"

  def writeOutPDBNames(self):
    tempSet = set([])
    FILE = open("AllPDBNames.txt", "wb")  
    for i in self.refPairs1:
      tempSet.add(i[:4])
    for i in tempSet:
      FILE.write(i + "\n")
    FILE.close

  def genSearch(self):
    str1 = ""
    while True:
      try:
        print "\nDefine a general search parameter"
        print " (Booleans are supported. Ex 'Hoogsteen AND DNA')"
        str1 = str(raw_input(" > "))
        print ""
        self.genSearchPDB = self.upClass.generalSearch(str1)
        break
      except ValueError: print "Please search for a name.\n"   

  def findModResi(self, mainC):
    str1 = ""
    while True:
      try:
        print "\nInput modified residue name (ex. pseudouridine)"
        str1 = str(raw_input(" > "))
        print ""
        if str1[:1] != "-":
          self.modResiPDB = self.upClass.searchModResi(str1)
          mainC.pruneListPDB(self.modResiPDB, "AND") 
        elif str1[:1] == "-":
          self.modResiPDB = self.upClass.searchModResi(str1[1:])
          mainC.pruneListPDB(self.modResiPDB, "NOT") 
        break
      except ValueError: print "Please search for a name.\n"

  def findLigands(self, mainC):
    str1 = ""
    while True:
      try:
        print "\nInput ligand name (ex. Paromomycin)"
        print " ( '-' before string defines NOT including)"
        str1 = str(raw_input(" > "))
        print ""
        if str1[:1] != "-":
          self.ligPDBlist = self.upClass.searchLigands(str1)
          mainC.pruneListPDB(self.ligPDBlist, "AND") 
        elif str1[:1] == "-":
          self.ligPDBlist = self.upClass.searchLigands(str1[1:])
          mainC.pruneListPDB(self.ligPDBlist, "NOT") 
        break
      except ValueError: print "Please search for a name.\n"

  def findChemID(self, mainC):
    str1 = ""
    while True:
      try:
        print "\nInput 1-3 character chemical ID (IE. MG for magnesium)"
        print " ( '-' before string defines NOT including)"
        str1 = str(raw_input(" > "))
        print ""
        if str1[:1] != "-":
          self.chemID = self.upClass.searchChemID(str1)
          mainC.pruneListPDB(self.chemID, "AND") 
        elif str1[:1] == "-":
          self.chemID = self.upClass.searchChemID(str1[1:])
          mainC.pruneListPDB(self.chemID, "NOT") 
        break
      except ValueError: print "Please search for a name.\n"

  def pStrucName(self, mainC):
    str1 = ""
    while True:
      try:
        print "\nInput parent structure type (ex. Ribosome)"
        print " ( '-' before string defines NOT including)"
        str1 = str(raw_input(" > "))
        print ""
        if str1[:1] != "-":
          self.pSName = self.upClass.searchStrucName(str1)
          mainC.pruneListPDB(self.pSName, "AND") 
        elif str1[:1] == "-":
          self.pSName = self.upClass.searchStrucName(str1[1:])
          mainC.pruneListPDB(self.pSName, "NOT") 
        break
      except ValueError: print "Please search for a name.\n"

  def writeShifts(self, shiftData):
    tempList = []
    allMotifs = set([])
    {tempList.append(self.allPairs[x]) for x in self.allPairs if x in self.refPairs1}

    now = datetime.datetime.now()
    outStr = (str(now.year) + str(now.month) + str(now.day) + 
                    str(now.hour) + str(now.minute) + str(now.second))

    self.writeOutFolder = self.curDir + outStr + "/"
    if not os.path.exists(self.writeOutFolder): os.makedirs(self.writeOutFolder) 
    FILE = open(((self.writeOutFolder + "%s.csv") % outStr), "w")
    self.writeOutName = ("%s.csv" % outStr)

      #Header line
    FILE.write("PDB ID" + "," + "SubPair Motif" + "," + "Pair Motif" + "," + "Pair ID" + "," 
                     + "Base A Strand" + "," + "Base B Strand" + "," + "Base A Type" + "," 
                     + "Base B Type" + "," + "Base A Subtype" + "," + "Base B Subtype" + "," 
                     + "Base A Num"  + "," + "Base B Num" + ","
                     + "Shear" + "," + "Stretch" + "," + "Stagger" + "," 
                     + "Buckle" + "," + "Propeller" + "," + "Opening" + "," 
                     + "Pucker A" +"," + "Pucker B" + "," 
                     + "V0 A" + "," + "V1 A" + "," + "V2 A"  + "," + "V3 A" + "," +  "V4 A" + "," 
                     + "V0 B" + ","  + "V1 B" + "," + "V2 B" + "," + "V3 B" + "," + "V4 B" + "," 
                     +  "Chi A"  + "," + "Chi B" + "," + "ChiConf A" +"," + "ChiConf B" + "," 
                     + "bbAlpha A" + ","  + "bbBeta A" + "," + "bbGamma A" + "," 
                     + "bbDelta A" + "," + "bbEpsilon A" + "," + "bbZeta A" + "," 
                     + "bbAlpha B" + ","  + "bbBeta B" + "," + "bbGamma B" + "," 
                     + "bbDelta B" + "," + "bbEpsilon B" + "," + "bbZeta B" + ","
                     + "C1'-C1'" + "," + "Num HBonds" + "," 
                     + "HBondAt1 A" + "," + "HBond1 Type" + "," + "HBondAt1 B" + "," 
                     + "HBond1 Dist"  + ","
                     + "HBondAt2 A" + "," + "HBond2 Type" + "," + "HBondAt2 B" + "," 
                     + "HBond2 Dist"  + ","
                     + "HBondAt3 A" + "," + "HBond3 Type" + "," + "HBondAt3 B" + "," 
                     + "HBond3 Dist" + ","
                     + "HBondAt4 A" + "," + "HBond4 Type" + "," + "HBondAt4 B" + "," 
                     + "HBond4 Dist" + ","
                     # pdbData starts here
                     + "Resolution" + "," + "pH" + "," + "Temp" + "," + "Expt Conditions" + ","
                     + "Pubmed ID" + ","
                     #Shifts start here
                     + "A_C1'" + "," + "A_C2" + "," + "A_C2'" + "," + "A_C3'" + "," + "A_C4" + "," 
                     + "A_C4'" + ","
                     + "A_C5" + "," + "A_C5'" + "," + "A_C6" + "," + "A_C8" + "," + "A_H1" + ","
                     + "A_H1'" + "," + "A_H2" + "," + "A_H2'" + "," + "A_H21" + "," + "A_H22" + ","
                     + "A_H3" + "," + "A_H3'" + "," + "A_H4'" + "," + "A_H41" + ","
                     + "A_H42" + "," + "A_H5" + "," + "A_H5'" + "," + "A_H5''" + "," + "A_H6" + ","
                     + "A_H61" + "," + "A_H62" + "," + "A_H8" + "," + "A_HO2''" + "," + "A_N1" + ","
                     + "A_N2" + "," + "A_N3" + "," + "A_N4" + "," + "A_N6'" + "," + "A_N7" + ","
                     + "A_N9" + "," + "A_P" + ","
                     + "A_stdC1'" + "," + "A_stdC2" + "," + "A_stdC2'" + "," + "A_stdC3'" + "," 
                     + "A_stdC4" + "," + "A_stdC4'" + "," + "A_stdC5" + "," + "A_stdC5'" + "," 
                     + "A_stdC6" + "," + "A_stdC8" + "," + "A_stdH1" + "," + "A_stdH1'" + "," 
                     + "A_stdH2" + "," + "A_stdH2'" + "," + "A_stdH21" + "," + "A_stdH22" + ","
                     + "A_stdH3" + "," + "A_stdH3'" + "," + "A_stdH4'" + "," + "A_stdH41" + ","
                     + "A_stdH42" + "," + "A_stdH5" + "," + "A_stdH5'" + "," + "A_stdH5''" + "," 
                     + "A_stdH6" + "," + "A_stdH61" + "," + "A_stdH62" + "," + "A_stdH8" + "," 
                     + "A_stdHO2''" + "," + "A_stdN1" + "," + "A_stdN2" + "," + "A_stdN3" + "," 
                     + "A_stdN4" + "," + "A_stdN6'" + "," + "A_stdN7" + "," + "A_stdN9" + "," 
                     + "A_stdP" + ","
                     + "B_C1'" + "," + "B_C2" + "," + "B_C2'" + "," + "B_C3'" + "," + "B_C4" + "," 
                     + "B_C4'" + ","
                     + "B_C5" + "," + "B_C5'" + "," + "B_C6" + "," + "B_C8" + "," + "B_H1" + ","
                     + "B_H1'" + "," + "B_H2" + "," + "B_H2'" + "," + "B_H21" + "," + "B_H22" + ","
                     + "B_H3" + "," + "B_H3'" + "," + "B_H4'" + "," + "B_H41" + ","
                     + "B_H42" + "," + "B_H5" + "," + "B_H5'" + "," + "B_H5''" + "," + "B_H6" + ","
                     + "B_H61" + "," + "B_H62" + "," + "B_H8" + "," + "B_HO2''" + "," + "B_N1" + ","
                     + "B_N2" + "," + "B_N3" + "," + "B_N4" + "," + "B_N6'" + "," + "B_N7" + ","
                     + "B_N9" + "," + "B_P" + ","
                     + "B_stdC1'" + "," + "B_stdC2" + "," + "B_stdC2'" + "," + "B_stdC3'" + "," 
                     + "B_stdC4" + "," + "B_stdC4'" + "," + "B_stdC5" + "," + "B_stdC5'" + "," 
                     + "B_stdC6" + "," + "B_stdC8" + "," + "B_stdH1" + "," + "B_stdH1'" + "," 
                     + "B_stdH2" + "," + "B_stdH2'" + "," + "B_stdH21" + "," + "B_stdH22" + ","
                     + "B_stdH3" + "," + "B_stdH3'" + "," + "B_stdH4'" + "," + "B_stdH41" + ","
                     + "B_stdH42" + "," + "B_stdH5" + "," + "B_stdH5'" + "," + "B_stdH5''" + "," 
                     + "B_stdH6" + "," + "B_stdH61" + "," + "B_stdH62" + "," + "B_stdH8" + "," 
                     + "B_stdHO2''" + "," + "B_stdN1" + "," + "B_stdN2" + "," + "B_stdN3" + "," 
                     + "B_stdN4" + "," + "B_stdN6'" + "," + "B_stdN7" + "," + "B_stdN9" + "," 
                     + "B_stdP" + ","
                     + "\n")
    
    for i in tempList:        
          #PDB name
      FILE.write(i.pName + ',')

          #Assign subpair motifs
      if i.bAsubtype < i.bBsubtype:
        i.subPMotif = (i.bAsubtype + i.bBsubtype)
        FILE.write(i.bAsubtype + i.bBsubtype + ",")
      else:
        i.subPMotif = (i.bBsubtype + i.bAsubtype)
        FILE.write(i.bBsubtype + i.bAsubtype + ",")

        #Write out base pairing motif
        #AA AC AG AP AT AU CC CG CT CU GG GP GU II IU PU UU
      if i.bAtype.upper() < i.bBtype.upper():       
        FILE.write(i.bAtype.upper() + i.bBtype.upper() + ",")
        allMotifs.add(i.bAtype.upper() + i.bBtype.upper())
      else:
        FILE.write(i.bBtype.upper() + i.bAtype.upper() + ",")
        allMotifs.add(i.bBtype.upper() + i.bAtype.upper())
        #Write readable format for pair ID
      # tempStr = (i.bAstrand + " (" + i.bAsubtype + " " + str(i.bAnum) + i.bAnumX + ") -- ("
      #                  + i.bBsubtype + " " + str(i.bBnum) + i.bBnumX + ") " + i.bBstrand)
      # FILE.write(i.pName + "_" + tempStr + ",")

      FILE.write(i.uniqueID + ",")

      FILE.write(i.bAstrand + "," + i.bBstrand + "," + i.bAtype + "," + i.bBtype + "," + 
                        i.bAsubtype + "," + i.bBsubtype + "," + str(i.bAnum) + i.bAnumX + "," + 
                        str(i.bBnum) + i.bBnumX + ",") 
        #Write local base pair parameters
      FILE.write(str(i.shear) + "," + str(i.stretch) + "," + str(i.stagger) + "," + str(i.buckle) + ","
                        + str(i.propeller) + "," + str(i.opening) + ",")
        #Write sugar parameters
      FILE.write(i.puckerA + "," + i.puckerB + "," + str(i.vee0a) + "," + str(i.vee1a) + "," 
                      + str(i.vee2a) + "," + str(i.vee3a) + "," + str(i.vee4a) + "," 
                      + str(i.vee0b) + "," + str(i.vee1b) + "," + str(i.vee2b) + "," 
                      + i.vee3b + "," + i.vee4b + "," )

        #Write backbone torsion angles
      FILE.write(str(i.bbChiA) + "," + str(i.bbChiB) + ",")
      if i.bbChiA != "---":
        if (math.fabs(float(i.bbChiA)) >= 90):
          FILE.write("Anti" + ",")
        else:
          FILE.write("Syn" + ",")
      else: FILE.write("---" + ",")
      if i.bbChiB != "---":
        if (math.fabs(float(i.bbChiB)) >= 90):
          FILE.write("Anti" + ",")
        else:
          FILE.write("Syn" + ",")
        
      else: FILE.write("---" + ",")
      FILE.write(i.bbAlphaA + "," + i.bbBetaA + "," + i.bbGammaA + "," + i.bbDeltaA + "," 
                      + i.bbEpsilonA + "," + i.bbZetaA + "," + i.bbAlphaB + "," + i.bbBetaB + ","
                      + i.bbGammaB + "," + i.bbDeltaB + "," + i.bbEpsilonB + "," + i.bbZetaB + ",")

        #Write out Virtual distances and angles
      FILE.write(str(i.c1c1) + ",")

        #Write out H-Bond parameters and values
      FILE.write(str(i.numBonds) + ",")
      if i.numBonds == 1:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
      if i.numBonds == 2:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")                
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
      if i.numBonds == 3:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")
        FILE.write(str(i.bondAtLHS[2]) + "," + str(i.bondType[2]) + "," + str(i.bondAtRHS[2])
                         + "," + str(i.bondDists[2]) + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
      if i.numBonds == 4:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")
        FILE.write(str(i.bondAtLHS[2]) + "," + str(i.bondType[2]) + "," + str(i.bondAtRHS[2])
                         + "," + str(i.bondDists[2]) + ",")
        FILE.write(str(i.bondAtLHS[3]) + "," + str(i.bondType[3]) + "," + str(i.bondAtRHS[3])
                         + "," + str(i.bondDists[3]) + ",")

      #Write structure data
      if self.pdbData[i.pName] == None:
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
      else:
        if self.pdbData[i.pName].resolution == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].resolution + ",")
        if self.pdbData[i.pName].phVal == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].phVal + ",")
        if self.pdbData[i.pName].temp == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].temp + ",")
        if self.pdbData[i.pName].conditions == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].conditions + ",")
        if self.pdbData[i.pName].pubmedId == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].pubmedId + ",")

            #0 PDB Name
            #1 Pair Motif Sorted Name
            #2 Base A
            #3 Base B
            #4 Temp
            #5 pH
            #6 Atom type
            #7 Chemical Shift
            #8 Shear, 9 stretch
            #10 Stagger, 11 buckle
            #12 Propeller, 13 opening
            #14 Sugar pucker
            #15 Chi Torsion

          #Process and write shift data
      aC1p, aC2, aC2p, aC3p = [], [], [], []
      aC4, aC4p, aC5, aC5p = [], [], [], []
      aC6, aC8, aH1, aH1p = [], [], [], []
      aH2, aH2p, aH21, aH22 = [], [], [], []
      aH3, aH3p, aH4p, aH41 = [], [], [], []
      aH42, aH5, aH5p, aH5pp = [], [], [], []
      aH6, aH61, aH62, aH8 = [], [], [], []
      aHO2p, aN1, aN2, aN3, aN4 = [], [], [], [], []
      aN6, aN7, aN9, aP = [], [], [], []
      bC1p, bC2, bC2p, bC3p = [], [], [], []
      bC4, bC4p, bC5, bC5p = [], [], [], []
      bC6, bC8, bH1, bH1p = [], [], [], []
      bH2, bH2p, bH21, bH22 = [], [], [], []
      bH3, bH3p, bH4p, bH41 = [], [], [], []
      bH42, bH5, bH5p, bH5pp = [], [], [], []
      bH6, bH61, bH62, bH8 = [], [], [], []
      bHO2p, bN1, bN2, bN3, bN4 = [], [], [], [], []
      bN6, bN7, bN9, bP = [], [], [], []

      for j in shiftData:
        if j[1] == i.pairMotif:
          if (self.exceedsZ6(1.5, i.shear, i.stretch, i.stagger,
                                            i.buckle, i.propeller, i.opening,
                                            j[8], j[9], j[10],
                                            j[11], j[12], j[13],
                                            1.1,1.75,1.75,5,5,5) == True):
            # if j[2] == i.bAsubtype:
            #   if (i.puckerA == j[14] and i.bbChiAconf == j[15]): 
            if ((j[2] == i.bAsubtype)
               and (i.puckerA == j[14] and i.bbChiAconf == j[15])): 
                if (j[6] == "C1'"): aC1p.append(j[7])
                if (j[6] == "C2"): aC2.append(j[7])
                if (j[6] == "C2'"): aC2p.append(j[7])
                if (j[6] == "C3'"): aC3p.append(j[7])
                if (j[6] == "C4"): aC4.append(j[7])
                if (j[6] == "C4'"): aC4p.append(j[7])
                if (j[6] == "C5"): aC5.append(j[7])         
                if (j[6] == "C5'"): aC5p.append(j[7]) 
                if (j[6] == "C6"): aC6.append(j[7])
                if (j[6] == "C8"): aC8.append(j[7])
                if (j[6] == "H1"): aH1.append(j[7])
                if (j[6] == "H1'"): aH1p.append(j[7])
                if (j[6] == "H2"): aH2.append(j[7])
                if (j[6] == "H2'"): aH2p.append(j[7])
                if (j[6] == "H21"): aH21.append(j[7])
                if (j[6] == "H22"): aH22.append(j[7])
                if (j[6] == "H3"): aH3.append(j[7])
                if (j[6] == "H3'"): aH3p.append(j[7])
                if (j[6] == "H4'"): aH4p.append(j[7])                  
                if (j[6] == "H41"): aH41.append(j[7])
                if (j[6] == "H5"): aH5.append(j[7])
                if (j[6] == "H5'"): aH5p.append(j[7])
                if (j[6] == "H5''"): aH5pp.append(j[7])
                if (j[6] == "H6"): aH6.append(j[7])
                if (j[6] == "H61"): aH61.append(j[7])
                if (j[6] == "H62"): aH62.append(j[7])
                if (j[6] == "H8"): aH8.append(j[7])
                if (j[6] == "HO2'"): aHO2p.append(j[7])
                if (j[6] == "N1"): aN1.append(j[7])
                if (j[6] == "N2"): aN2.append(j[7])               
                if (j[6] == "N3"): aN3.append(j[7]) 
                if (j[6] == "N4"): aN4.append(j[7]) 
                if (j[6] == "N6"): aN6.append(j[7]) 
                if (j[6] == "N7"): aN7.append(j[7]) 
                if (j[6] == "N9"): aN9.append(j[7]) 
                if (j[6] == "P"): aP.append(j[7]) 
            elif ((j[2] == i.bBsubtype) 
                    and (i.puckerB == j[14] and i.bbChiBconf == j[15])):
                if (j[6] == "C1'"): bC1p.append(j[7])
                if (j[6] == "C2"): bC2.append(j[7])
                if (j[6] == "C2'"): bC2p.append(j[7])
                if (j[6] == "C3'"): bC3p.append(j[7])
                if (j[6] == "C4"): bC4.append(j[7])
                if (j[6] == "C4'"): bC4p.append(j[7])
                if (j[6] == "C5"): bC5.append(j[7])         
                if (j[6] == "C5'"): bC5p.append(j[7]) 
                if (j[6] == "C6"): bC6.append(j[7])
                if (j[6] == "C8"): bC8.append(j[7])
                if (j[6] == "H1"): bH1.append(j[7])
                if (j[6] == "H1'"): bH1p.append(j[7])
                if (j[6] == "H2"): bH2.append(j[7])
                if (j[6] == "H2'"): bH2p.append(j[7])
                if (j[6] == "H21"): bH21.append(j[7])
                if (j[6] == "H22"): bH22.append(j[7])
                if (j[6] == "H3"): bH3.append(j[7])
                if (j[6] == "H3'"): bH3p.append(j[7])
                if (j[6] == "H4'"): bH4p.append(j[7])                  
                if (j[6] == "H41"): bH41.append(j[7])
                if (j[6] == "H5"): bH5.append(j[7])
                if (j[6] == "H5'"): bH5p.append(j[7])
                if (j[6] == "H5''"): bH5pp.append(j[7])
                if (j[6] == "H6"): bH6.append(j[7])
                if (j[6] == "H61"): bH61.append(j[7])
                if (j[6] == "H62"): bH62.append(j[7])
                if (j[6] == "H8"): bH8.append(j[7])
                if (j[6] == "HO2'"): bHO2p.append(j[7])
                if (j[6] == "N1"): bN1.append(j[7])
                if (j[6] == "N2"): bN2.append(j[7])               
                if (j[6] == "N3"): bN3.append(j[7]) 
                if (j[6] == "N4"): bN4.append(j[7]) 
                if (j[6] == "N6"): bN6.append(j[7]) 
                if (j[6] == "N7"): bN7.append(j[7]) 
                if (j[6] == "N9"): bN9.append(j[7]) 
                if (j[6] == "P"): bP.append(j[7]) 
      
      # Write out shifts here
      if len(aC1p) != 0:
        FILE.write(str(np.mean(aC1p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC2) != 0:
        FILE.write(str(np.mean(aC2)) + ",")
      else: FILE.write("---" + ",")

      if len(aC2p) != 0:
        FILE.write(str(np.mean(aC2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC3p) != 0:
        FILE.write(str(np.mean(aC3p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC4) != 0:
        FILE.write(str(np.mean(aC4)) + ",")
      else: FILE.write("---" + ",")

      if len(aC4p) != 0:
        FILE.write(str(np.mean(aC4p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC5) != 0:
        FILE.write(str(np.mean(aC5)) + ",")
      else: FILE.write("---" + ",")

      if len(aC5p) != 0:
        FILE.write(str(np.mean(aC5p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC6) != 0:
        FILE.write(str(np.mean(aC6)) + ",")
      else: FILE.write("---" + ",")
      
      if len(aC8) != 0:
        FILE.write(str(np.mean(aC8)) + ",")
      else: FILE.write("---" + ",")

      if len(aH1) != 0:
        FILE.write(str(np.mean(aH1)) + ",")
      else: FILE.write("---" + ",")

      if len(aH1p) != 0:
        FILE.write(str(np.mean(aH1p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH2) != 0:
        FILE.write(str(np.mean(aH2)) + ",")
      else: FILE.write("---" + ",")

      if len(aH2p) != 0:
        FILE.write(str(np.mean(aH2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH21) != 0:
        FILE.write(str(np.mean(aH21)) + ",")
      else: FILE.write("---" + ",")

      if len(aH22) != 0:
        FILE.write(str(np.mean(aH22)) + ",")
      else: FILE.write("---" + ",")

      if len(aH3) != 0:
        FILE.write(str(np.mean(aH3)) + ",")
      else: FILE.write("---" + ",")

      if len(aH3p) != 0:
        FILE.write(str(np.mean(aH3p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH4p) != 0:
        FILE.write(str(np.mean(aH4p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH41) != 0:
        FILE.write(str(np.mean(aH41)) + ",")
      else: FILE.write("---" + ",")

      if len(aH42) != 0:
        FILE.write(str(np.mean(aH42)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5) != 0:
        FILE.write(str(np.mean(aH5)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5p) != 0:
        FILE.write(str(np.mean(aH5p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5pp) != 0:
        FILE.write(str(np.mean(aH5pp)) + ",")
      else: FILE.write("---" + ",")

      if len(aH6) != 0:
        FILE.write(str(np.mean(aH6)) + ",")
      else: FILE.write("---" + ",")

      if len(aH61) != 0:
        FILE.write(str(np.mean(aH61)) + ",")
      else: FILE.write("---" + ",")

      if len(aH62) != 0:
        FILE.write(str(np.mean(aH62)) + ",")
      else: FILE.write("---" + ",")

      if len(aH8) != 0:
        FILE.write(str(np.mean(aH8)) + ",")
      else: FILE.write("---" + ",")

      if len(aHO2p) != 0:
        FILE.write(str(np.mean(aHO2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aN1) != 0:
        FILE.write(str(np.mean(aN1)) + ",")
      else: FILE.write("---" + ",")

      if len(aN2) != 0:
        FILE.write(str(np.mean(aN2)) + ",")
      else: FILE.write("---" + ",")

      if len(aN3) != 0:
        FILE.write(str(np.mean(aN3)) + ",")
      else: FILE.write("---" + ",")

      if len(aN4) != 0:
        FILE.write(str(np.mean(aN4)) + ",")
      else: FILE.write("---" + ",")

      if len(aN6) != 0:
        FILE.write(str(np.mean(aN6)) + ",")
      else: FILE.write("---" + ",")

      if len(aN7) != 0:
        FILE.write(str(np.mean(aN7)) + ",")
      else: FILE.write("---" + ",")

      if len(aN9) != 0:
        FILE.write(str(np.mean(aN9)) + ",")
      else: FILE.write("---" + ",")

      if len(aP) != 0:
        FILE.write(str(np.mean(aP)) + ",")
      else: FILE.write("---" + ",")

      if len(aC1p) != 0:
        FILE.write(str(np.std(aC1p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC2) != 0:
        FILE.write(str(np.std(aC2)) + ",")
      else: FILE.write("---" + ",")

      if len(aC2p) != 0:
        FILE.write(str(np.std(aC2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC3p) != 0:
        FILE.write(str(np.std(aC3p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC4) != 0:
        FILE.write(str(np.std(aC4)) + ",")
      else: FILE.write("---" + ",")

      if len(aC4p) != 0:
        FILE.write(str(np.std(aC4p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC5) != 0:
        FILE.write(str(np.std(aC5)) + ",")
      else: FILE.write("---" + ",")

      if len(aC5p) != 0:
        FILE.write(str(np.std(aC5p)) + ",")
      else: FILE.write("---" + ",")

      if len(aC6) != 0:
        FILE.write(str(np.std(aC6)) + ",")
      else: FILE.write("---" + ",")
      
      if len(aC8) != 0:
        FILE.write(str(np.std(aC8)) + ",")
      else: FILE.write("---" + ",")

      if len(aH1) != 0:
        FILE.write(str(np.std(aH1)) + ",")
      else: FILE.write("---" + ",")

      if len(aH1p) != 0:
        FILE.write(str(np.std(aH1p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH2) != 0:
        FILE.write(str(np.std(aH2)) + ",")
      else: FILE.write("---" + ",")

      if len(aH2p) != 0:
        FILE.write(str(np.std(aH2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH21) != 0:
        FILE.write(str(np.std(aH21)) + ",")
      else: FILE.write("---" + ",")

      if len(aH22) != 0:
        FILE.write(str(np.std(aH22)) + ",")
      else: FILE.write("---" + ",")

      if len(aH3) != 0:
        FILE.write(str(np.std(aH3)) + ",")
      else: FILE.write("---" + ",")

      if len(aH3p) != 0:
        FILE.write(str(np.std(aH3p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH4p) != 0:
        FILE.write(str(np.std(aH4p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH41) != 0:
        FILE.write(str(np.std(aH41)) + ",")
      else: FILE.write("---" + ",")

      if len(aH42) != 0:
        FILE.write(str(np.std(aH42)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5) != 0:
        FILE.write(str(np.std(aH5)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5p) != 0:
        FILE.write(str(np.std(aH5p)) + ",")
      else: FILE.write("---" + ",")

      if len(aH5pp) != 0:
        FILE.write(str(np.std(aH5pp)) + ",")
      else: FILE.write("---" + ",")

      if len(aH6) != 0:
        FILE.write(str(np.std(aH6)) + ",")
      else: FILE.write("---" + ",")

      if len(aH61) != 0:
        FILE.write(str(np.std(aH61)) + ",")
      else: FILE.write("---" + ",")

      if len(aH62) != 0:
        FILE.write(str(np.std(aH62)) + ",")
      else: FILE.write("---" + ",")

      if len(aH8) != 0:
        FILE.write(str(np.std(aH8)) + ",")
      else: FILE.write("---" + ",")

      if len(aHO2p) != 0:
        FILE.write(str(np.std(aHO2p)) + ",")
      else: FILE.write("---" + ",")

      if len(aN1) != 0:
        FILE.write(str(np.std(aN1)) + ",")
      else: FILE.write("---" + ",")

      if len(aN2) != 0:
        FILE.write(str(np.std(aN2)) + ",")
      else: FILE.write("---" + ",")

      if len(aN3) != 0:
        FILE.write(str(np.std(aN3)) + ",")
      else: FILE.write("---" + ",")

      if len(aN4) != 0:
        FILE.write(str(np.std(aN4)) + ",")
      else: FILE.write("---" + ",")

      if len(aN6) != 0:
        FILE.write(str(np.std(aN6)) + ",")
      else: FILE.write("---" + ",")

      if len(aN7) != 0:
        FILE.write(str(np.std(aN7)) + ",")
      else: FILE.write("---" + ",")

      if len(aN9) != 0:
        FILE.write(str(np.std(aN9)) + ",")
      else: FILE.write("---" + ",")

      if len(aP) != 0:
        FILE.write(str(np.std(aP)) + ",")
      else: FILE.write("---" + ",")

      if len(bC1p) != 0:
        FILE.write(str(np.mean(bC1p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC2) != 0:
        FILE.write(str(np.mean(bC2)) + ",")
      else: FILE.write("---" + ",")

      if len(bC2p) != 0:
        FILE.write(str(np.mean(bC2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC3p) != 0:
        FILE.write(str(np.mean(bC3p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC4) != 0:
        FILE.write(str(np.mean(bC4)) + ",")
      else: FILE.write("---" + ",")

      if len(bC4p) != 0:
        FILE.write(str(np.mean(bC4p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC5) != 0:
        FILE.write(str(np.mean(bC5)) + ",")
      else: FILE.write("---" + ",")

      if len(bC5p) != 0:
        FILE.write(str(np.mean(bC5p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC6) != 0:
        FILE.write(str(np.mean(bC6)) + ",")
      else: FILE.write("---" + ",")
      
      if len(bC8) != 0:
        FILE.write(str(np.mean(bC8)) + ",")
      else: FILE.write("---" + ",")

      if len(bH1) != 0:
        FILE.write(str(np.mean(bH1)) + ",")
      else: FILE.write("---" + ",")

      if len(bH1p) != 0:
        FILE.write(str(np.mean(bH1p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH2) != 0:
        FILE.write(str(np.mean(bH2)) + ",")
      else: FILE.write("---" + ",")

      if len(bH2p) != 0:
        FILE.write(str(np.mean(bH2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH21) != 0:
        FILE.write(str(np.mean(bH21)) + ",")
      else: FILE.write("---" + ",")

      if len(bH22) != 0:
        FILE.write(str(np.mean(bH22)) + ",")
      else: FILE.write("---" + ",")

      if len(bH3) != 0:
        FILE.write(str(np.mean(bH3)) + ",")
      else: FILE.write("---" + ",")

      if len(bH3p) != 0:
        FILE.write(str(np.mean(bH3p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH4p) != 0:
        FILE.write(str(np.mean(bH4p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH41) != 0:
        FILE.write(str(np.mean(bH41)) + ",")
      else: FILE.write("---" + ",")

      if len(bH42) != 0:
        FILE.write(str(np.mean(bH42)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5) != 0:
        FILE.write(str(np.mean(bH5)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5p) != 0:
        FILE.write(str(np.mean(bH5p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5pp) != 0:
        FILE.write(str(np.mean(bH5pp)) + ",")
      else: FILE.write("---" + ",")

      if len(bH6) != 0:
        FILE.write(str(np.mean(bH6)) + ",")
      else: FILE.write("---" + ",")

      if len(bH61) != 0:
        FILE.write(str(np.mean(bH61)) + ",")
      else: FILE.write("---" + ",")

      if len(bH62) != 0:
        FILE.write(str(np.mean(bH62)) + ",")
      else: FILE.write("---" + ",")

      if len(bH8) != 0:
        FILE.write(str(np.mean(bH8)) + ",")
      else: FILE.write("---" + ",")

      if len(bHO2p) != 0:
        FILE.write(str(np.mean(bHO2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bN1) != 0:
        FILE.write(str(np.mean(bN1)) + ",")
      else: FILE.write("---" + ",")

      if len(bN2) != 0:
        FILE.write(str(np.mean(bN2)) + ",")
      else: FILE.write("---" + ",")

      if len(bN3) != 0:
        FILE.write(str(np.mean(bN3)) + ",")
      else: FILE.write("---" + ",")

      if len(bN4) != 0:
        FILE.write(str(np.mean(bN4)) + ",")
      else: FILE.write("---" + ",")

      if len(bN6) != 0:
        FILE.write(str(np.mean(bN6)) + ",")
      else: FILE.write("---" + ",")

      if len(bN7) != 0:
        FILE.write(str(np.mean(bN7)) + ",")
      else: FILE.write("---" + ",")

      if len(bN9) != 0:
        FILE.write(str(np.mean(bN9)) + ",")
      else: FILE.write("---" + ",")

      if len(bP) != 0:
        FILE.write(str(np.mean(bP)) + ",")
      else: FILE.write("---" + ",")

      if len(bC1p) != 0:
        FILE.write(str(np.std(bC1p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC2) != 0:
        FILE.write(str(np.std(bC2)) + ",")
      else: FILE.write("---" + ",")

      if len(bC2p) != 0:
        FILE.write(str(np.std(bC2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC3p) != 0:
        FILE.write(str(np.std(bC3p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC4) != 0:
        FILE.write(str(np.std(bC4)) + ",")
      else: FILE.write("---" + ",")

      if len(bC4p) != 0:
        FILE.write(str(np.std(bC4p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC5) != 0:
        FILE.write(str(np.std(bC5)) + ",")
      else: FILE.write("---" + ",")

      if len(bC5p) != 0:
        FILE.write(str(np.std(bC5p)) + ",")
      else: FILE.write("---" + ",")

      if len(bC6) != 0:
        FILE.write(str(np.std(bC6)) + ",")
      else: FILE.write("---" + ",")
      
      if len(bC8) != 0:
        FILE.write(str(np.std(bC8)) + ",")
      else: FILE.write("---" + ",")

      if len(bH1) != 0:
        FILE.write(str(np.std(bH1)) + ",")
      else: FILE.write("---" + ",")

      if len(bH1p) != 0:
        FILE.write(str(np.std(bH1p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH2) != 0:
        FILE.write(str(np.std(bH2)) + ",")
      else: FILE.write("---" + ",")

      if len(bH2p) != 0:
        FILE.write(str(np.std(bH2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH21) != 0:
        FILE.write(str(np.std(bH21)) + ",")
      else: FILE.write("---" + ",")

      if len(bH22) != 0:
        FILE.write(str(np.std(bH22)) + ",")
      else: FILE.write("---" + ",")

      if len(bH3) != 0:
        FILE.write(str(np.std(bH3)) + ",")
      else: FILE.write("---" + ",")

      if len(bH3p) != 0:
        FILE.write(str(np.std(bH3p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH4p) != 0:
        FILE.write(str(np.std(bH4p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH41) != 0:
        FILE.write(str(np.std(bH41)) + ",")
      else: FILE.write("---" + ",")

      if len(bH42) != 0:
        FILE.write(str(np.std(bH42)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5) != 0:
        FILE.write(str(np.std(bH5)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5p) != 0:
        FILE.write(str(np.std(bH5p)) + ",")
      else: FILE.write("---" + ",")

      if len(bH5pp) != 0:
        FILE.write(str(np.std(bH5pp)) + ",")
      else: FILE.write("---" + ",")

      if len(bH6) != 0:
        FILE.write(str(np.std(bH6)) + ",")
      else: FILE.write("---" + ",")

      if len(bH61) != 0:
        FILE.write(str(np.std(bH61)) + ",")
      else: FILE.write("---" + ",")

      if len(bH62) != 0:
        FILE.write(str(np.std(bH62)) + ",")
      else: FILE.write("---" + ",")

      if len(bH8) != 0:
        FILE.write(str(np.std(bH8)) + ",")
      else: FILE.write("---" + ",")

      if len(bHO2p) != 0:
        FILE.write(str(np.std(bHO2p)) + ",")
      else: FILE.write("---" + ",")

      if len(bN1) != 0:
        FILE.write(str(np.std(bN1)) + ",")
      else: FILE.write("---" + ",")

      if len(bN2) != 0:
        FILE.write(str(np.std(bN2)) + ",")
      else: FILE.write("---" + ",")

      if len(bN3) != 0:
        FILE.write(str(np.std(bN3)) + ",")
      else: FILE.write("---" + ",")

      if len(bN4) != 0:
        FILE.write(str(np.std(bN4)) + ",")
      else: FILE.write("---" + ",")

      if len(bN6) != 0:
        FILE.write(str(np.std(bN6)) + ",")
      else: FILE.write("---" + ",")

      if len(bN7) != 0:
        FILE.write(str(np.std(bN7)) + ",")
      else: FILE.write("---" + ",")

      if len(bN9) != 0:
        FILE.write(str(np.std(bN9)) + ",")
      else: FILE.write("---" + ",")

      if len(bP) != 0:
        FILE.write(str(np.std(bP)) + ",")
      else: FILE.write("---" + ",")

      aC1p, aC2, aC2p, aC3p = [], [], [], []
      aC4, aC4p, aC5, aC5p = [], [], [], []
      aC6, aC8, aH1, aH1p = [], [], [], []
      aH2, aH2p, aH21, aH22 = [], [], [], []
      aH3, aH3p, aH4p, aH41 = [], [], [], []
      aH42, aH5, aH5p, aH5pp = [], [], [], []
      aH6, aH61, aH62, aH8 = [], [], [], []
      aHO2p, aN1, aN2, aN3, aN4 = [], [], [], [], []
      aN6, aN7, aN9, aP = [], [], [], []
      bC1p, bC2, bC2p, bC3p = [], [], [], []
      bC4, bC4p, bC5, bC5p = [], [], [], []
      bC6, bC8, bH1, bH1p = [], [], [], []
      bH2, bH2p, bH21, bH22 = [], [], [], []
      bH3, bH3p, bH4p, bH41 = [], [], [], []
      bH42, bH5, bH5p, bH5pp = [], [], [], []
      bH6, bH61, bH62, bH8 = [], [], [], []
      bHO2p, bN1, bN2, bN3, bN4 = [], [], [], [], []
      bN6, bN7, bN9, bP = [], [], [], []

          #Write terminal line delimiter
      FILE.write('\n')    
    FILE.close

  def writeOut(self):
    tempList = []
    allMotifs = set([])
    {tempList.append(self.allPairs[x]) for x in self.allPairs if x in self.refPairs1}
    now = datetime.datetime.now()
    outStr = (str(now.year) + str(now.month) + str(now.day) + 
                    str(now.hour) + str(now.minute) + str(now.second))

    self.writeOutFolder = self.curDir + outStr + "/"
    if not os.path.exists(self.writeOutFolder): os.makedirs(self.writeOutFolder) 
    FILE = open(((self.writeOutFolder + "%s.csv") % outStr), "w")
    self.writeOutName = ("%s.csv" % outStr)

      #Header line
    FILE.write("PDB ID" + "," + "SubPair Motif" + "," + "Pair Motif" + "," + "Pair ID" + "," + "Base A Strand" + ","
                     + "Base B Strand" + "," + "Base A Type" + "," + "Base B Type" + "," 
                     + "Base A Subtype" + "," + "Base B Subtype" + "," + "Base A Num"  + "," 
                     + "Base B Num" + ","
                     + "Shear" + "," + "Stretch" + "," + "Stagger" + "," 
                     + "Buckle" + "," + "Propeller" + "," + "Opening" + "," 
                     + "Pucker A" +"," + "Pucker B" + "," 
                     + "V0 A" + "," + "V1 A" + "," + "V2 A"  + "," + "V3 A" + "," +  "V4 A" + "," 
                     + "V0 B" + ","  + "V1 B" + "," + "V2 B" + "," + "V3 B" + "," + "V4 B" + "," 
                     +  "Chi A"  + "," + "Chi B" + "," + "ChiConf A" +"," + "ChiConf B" + "," 
                     + "bbAlpha A" + ","  + "bbBeta A" + "," + "bbGamma A" + "," 
                     + "bbDelta A" + "," + "bbEpsilon A" + "," + "bbZeta A" + "," 
                     + "bbAlpha B" + ","  + "bbBeta B" + "," + "bbGamma B" + "," 
                     + "bbDelta B" + "," + "bbEpsilon B" + "," + "bbZeta B" + ","
                     + "C1'-C1'" + "," + "Num HBonds" + "," 
                     + "HBondAt1 A" + "," + "HBond1 Type" + "," + "HBondAt1 B" + "," 
                     + "HBond1 Dist"  + ","
                     + "HBondAt2 A" + "," + "HBond2 Type" + "," + "HBondAt2 B" + "," 
                     + "HBond2 Dist"  + ","
                     + "HBondAt3 A" + "," + "HBond3 Type" + "," + "HBondAt3 B" + "," 
                     + "HBond3 Dist" + ","
                     + "HBondAt4 A" + "," + "HBond4 Type" + "," + "HBondAt4 B" + "," 
                     + "HBond4 Dist" + ","
                     # pdbData starts here
                     + "Resolution" + "," + "pH" + "," + "Temp" + "," + "Expt Conditions" + ","
                     + "Pubmed ID" + ","
                     + "\n")
    
    for i in tempList:        
          #PDB name
      FILE.write(i.pName + ',')

          #Assign subpair motifs
      if i.bAsubtype < i.bBsubtype:
        i.subPMotif = (i.bAsubtype + i.bBsubtype)
        FILE.write(i.bAsubtype + i.bBsubtype + ",")
      else:
        i.subPMotif = (i.bBsubtype + i.bAsubtype)
        FILE.write(i.bBsubtype + i.bAsubtype + ",")

        #Write out base pairing motif
        #AA AC AG AP AT AU CC CG CT CU GG GP GU II IU PU UU
      if i.bAtype.upper() < i.bBtype.upper():       
        FILE.write(i.bAtype.upper() + i.bBtype.upper() + ",")
        allMotifs.add(i.bAtype.upper() + i.bBtype.upper())
      else:
        FILE.write(i.bBtype.upper() + i.bAtype.upper() + ",")
        allMotifs.add(i.bBtype.upper() + i.bAtype.upper())
        #Write readable format for pair ID
      # tempStr = (i.bAstrand + " (" + i.bAsubtype + " " + str(i.bAnum) + i.bAnumX + ") -- ("
      #                  + i.bBsubtype + " " + str(i.bBnum) + i.bBnumX + ") " + i.bBstrand)
      # FILE.write(i.pName + "_" + tempStr + ",")
      FILE.write(i.uniqueID + ",")

      FILE.write(i.bAstrand + "," + i.bBstrand + "," + i.bAtype + "," + i.bBtype + "," + 
                        i.bAsubtype + "," + i.bBsubtype + "," + str(i.bAnum) + i.bAnumX + "," + 
                        str(i.bBnum) + i.bBnumX + ",") 
        #Write local base pair parameters
      FILE.write(str(i.shear) + "," + str(i.stretch) + "," + str(i.stagger) + "," + str(i.buckle) + ","
                        + str(i.propeller) + "," + str(i.opening) + ",")
        #Write sugar parameters
      FILE.write(i.puckerA + "," + i.puckerB + "," + str(i.vee0a) + "," + str(i.vee1a) + "," 
                      + str(i.vee2a) + "," + str(i.vee3a) + "," + str(i.vee4a) + "," 
                      + str(i.vee0b) + "," + str(i.vee1b) + "," + str(i.vee2b) + "," 
                      + i.vee3b + "," + i.vee4b + "," )

        #Write backbone torsion angles
      FILE.write(str(i.bbChiA) + "," + str(i.bbChiB) + ",")
      if i.bbChiA != "---":
        if (math.fabs(float(i.bbChiA)) >= 90):
          FILE.write("Anti" + ",")
        else:
          FILE.write("Syn" + ",")
      else: FILE.write("---" + ",")
      if i.bbChiB != "---":
        if (math.fabs(float(i.bbChiB)) >= 90):
          FILE.write("Anti" + ",")
        else:
          FILE.write("Syn" + ",")
        
      else: FILE.write("---" + ",")
      FILE.write(i.bbAlphaA + "," + i.bbBetaA + "," + i.bbGammaA + "," + i.bbDeltaA + "," 
                      + i.bbEpsilonA + "," + i.bbZetaA + "," + i.bbAlphaB + "," + i.bbBetaB + ","
                      + i.bbGammaB + "," + i.bbDeltaB + "," + i.bbEpsilonB + "," + i.bbZetaB + ",")

        #Write out Virtual distances and angles
      FILE.write(str(i.c1c1) + ",")

        #Write out H-Bond parameters and values
      FILE.write(str(i.numBonds) + ",")
      if i.numBonds == 1:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
      if i.numBonds == 2:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")                
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        # if i.pName == "1FUF":
          # print "----", i.uniqueID, "----"
          # print ("1: ", str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
          #                + "," + str(i.bondDists[0]) + ",")
          # print ("2: ", str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
          #                + "," + str(i.bondDists[1]) + ",")    
      if i.numBonds == 3:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")
        FILE.write(str(i.bondAtLHS[2]) + "," + str(i.bondType[2]) + "," + str(i.bondAtRHS[2])
                         + "," + str(i.bondDists[2]) + ",")
        FILE.write("---" + "," + "---" + "," + "---" + "," + "---" + ",")
        # if i.pName == "1FUF":
          # print "----", i.uniqueID, "----"
          # print ("1: ", str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
          #                + "," + str(i.bondDists[0]) + ",")
          # print ("2: ", str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
          #                + "," + str(i.bondDists[1]) + ",")    
          # print ("3: ", str(i.bondAtLHS[2]) + "," + str(i.bondType[2]) + "," + str(i.bondAtRHS[2])
          #                + "," + str(i.bondDists[2]) + ",")
      if i.numBonds == 4:
        FILE.write(str(i.bondAtLHS[0]) + "," + str(i.bondType[0]) + "," + str(i.bondAtRHS[0])
                         + "," + str(i.bondDists[0]) + ",")
        FILE.write(str(i.bondAtLHS[1]) + "," + str(i.bondType[1]) + "," + str(i.bondAtRHS[1])
                         + "," + str(i.bondDists[1]) + ",")
        FILE.write(str(i.bondAtLHS[2]) + "," + str(i.bondType[2]) + "," + str(i.bondAtRHS[2])
                         + "," + str(i.bondDists[2]) + ",")
        FILE.write(str(i.bondAtLHS[3]) + "," + str(i.bondType[3]) + "," + str(i.bondAtRHS[3])
                         + "," + str(i.bondDists[3]) + ",")

      #Write structure data
      if self.pdbData[i.pName] == None:
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
        FILE.write("null" + ",")
      else:
        if self.pdbData[i.pName].resolution == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].resolution + ",")
        if self.pdbData[i.pName].phVal == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].phVal + ",")
        if self.pdbData[i.pName].temp == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].temp + ",")
        if self.pdbData[i.pName].conditions == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].conditions + ",")
        if self.pdbData[i.pName].pubmedId == None: FILE.write("null" + ",")
        else: FILE.write(self.pdbData[i.pName].pubmedId + ",")
          #Write terminal line delimiter
      FILE.write('\n')    
    FILE.close