import os, sys, string, marshal, math
import cPickle as pickle
import strucs, parse, updateBPF, run3dna
import datetime

class Control:
  def __init__(self):
      #all base pairs
    self.allpairs = []
      #New pairs processed from an incremental update
    self.newpairs = []
      #PDB data structure objects list, incl temp, pH ,etc
    self.pdbDataObs = {}
      #Store new PDB data (temp, pH, etc)
    self.newpdbdata = {}
      #Set of names of PDB structures with exp info like temp, pH, etc
    self.pdbExpNames = set([])
      #Stores a list of structures with no pairs that need not be re-downloaded
    self.nopairs = set([])
      #Just a temporary set of PDB names, may be deleted on the fly
    self.tempnames = set([])
      #Dictionary of hairpin classes
    self.hairpins = {}
      #List of new hairpin classes
    self.newHairpins = []
    self.rawData = []
    self.allMotifs = set([])
      #Set of PDB names from names.dat
    self.nameDatPairs = set([])
      #Set of PDB names pulled from the PDB, fresh
    self.newPDBNames = set([])
      #Difference between names.dat and newPDBNames
    self.pdbDiffName = set([])
    self.curDir = os.getcwd()
    self.realDir = os.path.dirname(os.path.realpath(__file__))
    self.datDir = os.path.join(self.realDir, "dat/")

    self.tempPDB = os.path.join(self.realDir, "tempPDB/")
    # self.tempPDB = "/Users/kimseyij/Del/PDB_NA_XRAY_3A/"
    self.mbpDir, self.outpDir, self.dssrDir = None, None, None
      #Define Run3DNA class
    self.r3dna = run3dna.Run3DNA()
      #Define parser class
    self.parse = parse.ParserC()
      #Define an update class
    self.up = updateBPF.UpdateBPF()

    #This function, when called, blanks out pairs.dat and names.dat files
    #by writing empty to lists to each
  def buildBlank(self, pairspath, namepath, hairpinpath):
    blankNames = set([])
    blankPairs = {}
    blankPDBStruc = {}
    FILE = open(pairspath, "wb")
    pickle.dump(blankPairs, FILE, protocol=2)
    FILE.close()
    FILE = open(namepath, "wb")
    pickle.dump(blankNames, FILE, protocol=2)
    FILE.close()    
    FILE = open(hairpinpath, "wb")
    pickle.dump(blankNames, FILE, protocol=2)
    FILE.close()    

    #Rebuilds only PDB data pulled from the PDB servers
    #IE Temp, pH, etc
    # allNames - list of all PDB names
    # pdbDataPath - path to pdbData.dat file
  def buildPDBData(self, allNames, pdbDataPath):
    pdbdata = {}
    count = 0
    for i in allNames:
      count += 1
      print "Now working on:", i
      print "  (Number", count, "of 5524)"
      pdbC = strucs.pdbDataC()
      pdbC = self.up.pPdbData(i, pdbC)
      pdbdata[i] = pdbC
    FILE = open(pdbDataPath, "wb")
    pickle.dump(pdbdata, FILE, protocol=2)
    FILE.close()

    #Builds a names.dat file using user-defined parameters
    #Users define a PDB search criterion (molecular structure, etc) and
    #a list of the PDB names that match this search are exported to names.dat
  def buildNames(self, namepath, resolution, upmethod):
    self.newPDBNames = self.up.pullPDBNames(upmethod, resolution)
      #Write out all the names
    FILE = open(namepath, "wb")
    pickle.dump(self.newPDBNames, FILE, protocol=2)
    FILE.close()

    #upmethod -- defines the update method, xray, nmr, or both  --or DEV
    #Upmethod gets passed to updateBPF, which tells it which PDB names to
    #retrieve from the PDB server
    #listpath - directory of list to be used to update
    #purge - bool to make program overwrite existing files
    #err_logp is the error log path to write errors to
  def fullUp(self, namepath, pairsdat, resolution, pdbdatapath, 
                    hairpinpath, upmethod, listpath, purge=False,
                    err_logp = "error.log"): 
      #Get list of previous PDB files
    FILE = open(namepath, "rb")
    self.nameDatPairs = pickle.load(FILE)
    FILE.close()

      #Get dictionary of previous hairpin classes
    FILE = open(hairpinpath, "rb")
    self.hairpins = pickle.load(FILE)
    FILE.close()

      #Get dictionary of previous PDB data files (IE. temp, etc)
    FILE = open(pdbdatapath, "rb")
    self.pdbDataObs = pickle.load(FILE)
    FILE.close()

      #Checks for PDB names of XRAY method with specified resolution
    self.newPDBNames = self.up.pullPDBNames(upmethod, resolution, listpath)

      #find the difference in new vs old, if any
    self.pdbDiffName = self.newPDBNames - self.nameDatPairs

    # Write out list of previously curated PDBs existing in database
    update_old_pdb = os.path.join(listpath, "Update-Previous_PDB_Database.txt")
    FILE = open(update_old_pdb, "wb")
    for l in sorted(list(self.nameDatPairs)):
      FILE.write(l + "\n")
    FILE.close()

    # Write out list of new PDBs to download and parse
    update_diff_pdb = os.path.join(listpath, "Update-Only_New_PDB_Requests.txt")
    FILE = open(update_diff_pdb, "wb")
    for l in sorted(list(self.pdbDiffName)):
      FILE.write(l + "\n")
    FILE.close()

      #Add PDB name here if you need it fixed selectively
    # self.pdbDiffName.add("1S72")
    # self.pdbDiffName.add("3UYD")

    #if the difference between the two is not zero, run the incremental update.
    if len(self.pdbDiffName) != 0:
      self.nameDatPairs |= self.newPDBNames
      #Download the difference first
      #Make a folder to keep the downloaded PDB's
      if not os.path.exists(self.tempPDB):
        os.makedirs(self.tempPDB)
      
      remove_names = set([])
      # Go through list of PDBs and download them as needed  
      for i in self.pdbDiffName:
        # Try to download any missing PDB files
        pdb_exists = self.up.curlOut(self.tempPDB, i, purge=purge, err_logp=err_logp)

        if pdb_exists:
          #Parse PDB data like temp, etc
          #This is separate from structure parsing. This just grabs
          #Experimental data from the PDB like pH and temp
          pdbC = strucs.pdbDataC()
          pdbC = self.up.pPdbData(i, pdbC)
          self.newpdbdata[i] = pdbC

            #Generate new hairpin classes, as needed
          hpC = strucs.hairpinC()
          self.newHairpins.append(hpC)
        else:
          # write out error
          erf = open(err_logp, "ab")
          erf.write('{:%Y-%m-%d %Hh%Mm%Ss}: '.format(datetime.datetime.now()))
          errstr = "PDB %s does not exist locally, cannot process it. Removing it from update.\n" % i
          erf.write(errstr)
          erf.close()
          print errstr
          remove_names.add(i) # set of names to remove from PBDDiffName
      # Remove PDB file from list
      self.pdbDiffName -= remove_names

      # Write out list of PDBs that are wanted but cannot get downloaded
      # or processed
      update_err_pdb = os.path.join(listpath, "Update-Error_PDB_Not_Local.txt")
      FILE = open(update_err_pdb, "wb")
      for l in sorted(list(remove_names)):
        FILE.write(l + "\n")
      FILE.close()

        #Run x3DNA and DSSR
      self.mbpDir, self.outpDir, self.dssrDir = self.r3dna.genFiles(self.tempPDB, self.pdbDiffName,
                                                                    purge=purge)
        #Parse the new files only to newpairs
      self.parse.parseOutp(self.outpDir, self.pdbDiffName, self.newpairs)

        #Generate a list of PDB names ONLY which have base pairs
      [self.tempnames.add(x.pName) for x in self.newpairs]

        #Now, assign some more descriptors to each base pair incl.
        #Unique ID, pair motif, chi descriptor, etc
      for i in self.newpairs:                
          #Assign subpair motifs
        if i.bAsubtype < i.bBsubtype:
          i.subPMotif = (i.bAsubtype + i.bBsubtype)
        else:
          i.subPMotif = (i.bBsubtype + i.bAsubtype)
          #Assign sorted pair motif
        if i.bAtype.upper() < i.bBtype.upper():       
          i.pairMotif = (i.bAtype + i.bBtype)
          self.allMotifs.add(i.bAtype.upper() + i.bBtype.upper())
        else:
          i.pairMotif = (i.bBtype + i.bAtype)
          self.allMotifs.add(i.bBtype.upper() + i.bAtype.upper())

          #Unique ID just for DSSR
        i.uniqueID = (i.pName + "_" + i.bAstrand + i.bAsubtype + str(i.bAnum) 
                              + i.bAnumX + "_" + i.bBstrand + i.bBsubtype + str(i.bBnum)
                              + i.bBnumX)

          #Write backbone torsion angles#
        if i.bbChiA != "---":
          if (math.fabs(float(i.bbChiA)) >= 90):
            i.bbChiAconf = "Anti"
          else:
            i.bbChiAconf = "Syn"
        if i.bbChiB != "---":
          if (math.fabs(float(i.bbChiB)) >= 90):
            i.bbChiBconf = "Anti"
          else:
            i.bbChiBconf = "Syn"

        #Temporary dictionary of new pairs based on uniqueID
      newPdict =  {}
      for i in self.newpairs:
        newPdict[i.uniqueID] = i

        #Parse .dssr files separately
        #DSSR files contain secondary structure, etc, information that is correlated
        #to the outp data. So to process DSSR files, only need to process those which
        #exist in outp files.      
      # self.parse.parseDSSR(self.dssrDir, (self.pdbDiffName & self.tempnames), newPdict,
      #                                      self.newHairpins)

      print "  [Serializing base-pairs]"

        #Create new dictionary to store re-assigned pair classes
      tempDict = {}

        #Reassign dictionary key and unique ID
      for i in newPdict:
           #Write readable format for pair ID
        newPdict[i].uniqueID = (newPdict[i].pName + "_" + newPdict[i].pairMotif + "_" +
                               newPdict[i].bAstrand + "_" + newPdict[i].bAsubtype + "_" + 
                               str(newPdict[i].bAnum) + newPdict[i].bAnumX + "_" +
                               newPdict[i].bBstrand + "_" + newPdict[i].bBsubtype + "_" + 
                               str(newPdict[i].bBnum) + newPdict[i].bBnumX)
        tempDict[newPdict[i].uniqueID] = newPdict[i]

        #Unassign newPdict in prep for name-reassignment
      newPdict =  {}

        #Combine the .dat base pairs with the new base pairs (tempDict)
      FILE = open(pairsdat, "rb")
      listAllPairs = dict((pickle.load(FILE).items()) + tempDict.items())
      FILE.close()

        #Write out the newly combined base pair data files
      FILE = open(pairsdat, "wb")
      pickle.dump(listAllPairs, FILE, protocol=2)
      FILE.close()

        #Write out all the names
      FILE = open(namepath, "wb")
      pickle.dump(self.nameDatPairs, FILE, protocol=2)
      FILE.close()

        #Update PDB data files with new ones before writing out
      self.pdbDataObs.update(self.newpdbdata)
      FILE = open(pdbdatapath, "wb")
      pickle.dump(self.pdbDataObs, FILE, protocol=2)
      FILE.close()

        #Update hairpins
      self.hairpins.update(self.newHairpins)
      FILE = open(hairpinpath, "wb")
      pickle.dump(self.hairpins, FILE, protocol=2)
      FILE.close()

    else: print "No need to update."