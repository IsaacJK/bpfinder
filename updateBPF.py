#!/usr/bin/env python

import sys, os, urllib, urllib2
import xml.etree.ElementTree as ET
# import pdbsplitter

class Utilities:
  def __init__(self):
    pass

  def checkStatus(self, path, name):
    # Check if database list exists.
    try:
      with open(path + name) as f: 
        #If it exists
        return(True)
    except IOError as e:
        #if it does not exist
      return(False)        

class UpdateBPF: 
  def __init__(self):

    # initialize the server
#    self.server = SOAPpy.SOAPProxy("http://www.pdb.org/pdb/services/pdbws")
    self.server2 = "http://www.rcsb.org/pdb/rest/search"

    # Directories
    self.curDir = os.getcwd()
    self.tmpPdbDir = os.path.join(self.curDir, "tempPDB/")

    # Lists of NMR strucs
    # self.nmr_RNA_multi = []
    # self.nmr_DNA_multi = []
    # self.nmr_HYB_multi = []
    # self.split_names = set([])
    
    #Combined Lists
    # self.all_multi = set([]) #nmr multimodel
    self.all_pdb = set([]) #List of all PDB names, no split files

    # define generic XML query for xray strucs
    self.xml = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>%s</containsProtein>
        <containsDna>%s</containsDna>
        <containsRna>%s</containsRna>
        <containsHybrid>%s</containsHybrid>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
        <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
        <refine.ls_d_res_high.min>%s</refine.ls_d_res_high.min>
        <refine.ls_d_res_high.max>%s</refine.ls_d_res_high.max>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>'''

    self.xmlNmr = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>%s</containsProtein>
        <containsDna>%s</containsDna>
        <containsRna>%s</containsRna>
        <containsHybrid>%s</containsHybrid>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <queryType>org.pdb.query.simple.ModelCountQuery</queryType>
        <mvStructure.modelCount.comparator>%s</mvStructure.modelCount.comparator>
        <mvStructure.modelCount.value>%s</mvStructure.modelCount.value>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>2</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
        <mvStructure.expMethod.value>SOLUTION NMR</mvStructure.expMethod.value>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>'''

    self.xmlNmrNoModels = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>%s</containsProtein>
        <containsDna>%s</containsDna>
        <containsRna>%s</containsRna>
        <containsHybrid>%s</containsHybrid>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
        <description>Experimental Method is SOLUTION NMR</description>
        <mvStructure.expMethod.value>SOLUTION NMR</mvStructure.expMethod.value>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>'''

  def parseXmlQuery(self, inArray):
    data = []
    for answer_dict in inArray:
      data.append(answer_dict)  
    return (data)

  # This removes duplicates from an input list (seq)
  def remDupes(self, seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]

  def fileExists(self, filePath):
    try:
      with open(filePath) as f: return(True)
    except IOError as e: return(False)

    #This function uses the program curl to download needed PDB files
    # Purge flag will over-write existing file exclusion
  def curlOut(self, path, fileName, purge=False):
    curlCmd = 'curl -o "%s%s.pdb" http://files.rcsb.org/view/%s.pdb'
    temp_pdb_path = os.path.join(path, fileName.strip() + ".pdb")
    if self.fileExists(temp_pdb_path) == True and purge == False:
      pass
    else:
      print "Downloading ( %s )" % fileName   
      os.system(curlCmd % (path, fileName.strip(),
                           fileName.strip()))

      #This function returns a list of PDB names that match the given search criteria.
      #Pass method: XRAY, NMR, BOTH, or DEV
      #DEV is a custom method that just returns HYB XRAY structures, to give
      #  the user a smaller pool to work with
      #Resolution : resolution value for XRAY structures
      #Function returns a combined set of PDB names (unique) that
      #  match the given search criterion for DNA, RNA and HYB macromolecules
  def pullPDBNames(self, method, resolution, listpath):
      #Want XRAY structures of resolution range (0 - RESOLUTION)
      #Requests RNA names, then DNA, the HYB.
      #Combines these names to a set, no duplicates.
      #Returns this list to caller
    if method == "XRAY":
        #Request for RNA-Xray
      req = urllib2.Request(self.server2, (self.xml % ("?", "?", "Y", "?", "0", resolution)))
      post = urllib2.urlopen(req)
      xRNA = [x.strip() for x in post]

        #Request for DNA-Xray 
      req = urllib2.Request(self.server2, (self.xml % ("?", "Y", "?", "?", "0", resolution)))
      post = urllib2.urlopen(req)
      xDNA = [x.strip() for x in post]

        #Request for HYB-Xray
      req = urllib2.Request(self.server2, (self.xml % ("?", "?", "?", "Y", "0", resolution)))
      post = urllib2.urlopen(req)
      xHYB = [x.strip() for x in post]

        #Parse the XML query and set the names equal to a list
      xRNA = self.parseXmlQuery(xRNA)
      xDNA = self.parseXmlQuery(xDNA)
      xHYB = self.parseXmlQuery(xHYB)

        #Combine the above lists, and cast them as a set to remove duplicates
      self.all_pdb |= set(xRNA + xDNA + xHYB)

      return(self.all_pdb)

      #Same as XRAY, but for NMR only
    elif method == "NMR":
        #Request for RNA-NMR
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "?", "Y", "?")))
      post = urllib2.urlopen(req)
      nRNA = [x.strip() for x in post]

        #Request for DNA-NMR 
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "Y", "?", "?")))
      post = urllib2.urlopen(req)
      nDNA = [x.strip() for x in post]

        #Request for HYB-NMR
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "?", "?", "Y")))
      post = urllib2.urlopen(req)
      nHYB = [x.strip() for x in post]

      nRNA = self.parseXmlQuery(nRNA)
      nDNA = self.parseXmlQuery(nDNA)
      nHYB = self.parseXmlQuery(nHYB)

      self.all_pdb |= set(nRNA + nDNA + nHYB)

      return(self.all_pdb)

      #Same as above, but for both XRAY and NMR structures
    elif method == "BOTH":
        #Request for RNA-Xray
      req = urllib2.Request(self.server2, (self.xml % ("?", "?", "Y", "?", "0", resolution)))
      post = urllib2.urlopen(req)
      xRNA = [x.strip() for x in post]

        #Request for DNA-Xray 
      req = urllib2.Request(self.server2, (self.xml % ("?", "Y", "?", "?", "0", resolution)))
      post = urllib2.urlopen(req)
      xDNA = [x.strip() for x in post]

        #Request for HYB-Xray
      req = urllib2.Request(self.server2, (self.xml % ("?", "?", "?", "Y", "0", resolution)))
      post = urllib2.urlopen(req)
      xHYB = [x.strip() for x in post]

      xRNA = self.parseXmlQuery(xRNA)
      xDNA = self.parseXmlQuery(xDNA)
      xHYB = self.parseXmlQuery(xHYB)

        #Request for RNA-NMR
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "?", "Y", "?")))
      post = urllib2.urlopen(req)
      nRNA = [x.strip() for x in post]

        #Request for DNA-NMR 
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "Y", "?", "?")))
      post = urllib2.urlopen(req)
      nDNA = [x.strip() for x in post]

        #Request for HYB-NMR
      req = urllib2.Request(self.server2, (self.xmlNmrNoModels % ("?", "?", "?", "Y")))
      post = urllib2.urlopen(req)
      nHYB = [x.strip() for x in post]

      nRNA = self.parseXmlQuery(nRNA)
      nDNA = self.parseXmlQuery(nDNA)
      nHYB = self.parseXmlQuery(nHYB)
        #Combine both XRAY and NMR lists
      self.all_pdb |= set(xRNA + xDNA + xHYB + nRNA + nDNA + nHYB)

      return(self.all_pdb)

      #Used for DEV mode. Returns a small list of Xray PDB names to be processed,
      #This gives the dev a much quicker load time
    elif method == "DEV":

        #Request for HYB-Xray
      req = urllib2.Request(self.server2, (self.xml % ("?", "?", "?", "Y", "0", resolution)))
      post = urllib2.urlopen(req)
      xHYB = [x.strip() for x in post]
      xHYB = self.parseXmlQuery(xHYB)
      self.all_pdb |= set(xHYB)
      return(self.all_pdb)

      #Takes the path to list file and generates a set of PDB names
      #Returns this list
    elif method == "LIST":
      FILE = open(listpath, "rU")
      for i in FILE:
        self.all_pdb.add(i.strip().upper())
        print i.strip().upper()
      FILE.close()

      return(self.all_pdb)
      sys.exit(0)

      #This function is used to split multiple model PDB files
      #before they are processed by 3DNA
  # def runSplitter(self):
  #   if len(self.all_multi) < 1: print "Multi-model PDB list is empty"
  #   else:
  #     split = pdbsplitter.PDBSplitter(self.all_multi, self.tmpPdbDir)
  #     split.runSplitter()
  #     self.split_names = set(split.split_names)

  def returnList(self, selPars, loRes, hiRes):
    seleNames = set([])
    seleNamesNMR = set([])
    tempList = []
      #Request for list structures
      #?1 ?2 ?3 ?4 ?5 ?6
      #Protein DNA RNA Hybrid MinRes MaxRes
    
    if selPars[4] == "N":
      req = urllib2.Request(self.server2,
                            (self.xml % (selPars[2], selPars[0], selPars[1], selPars[3], loRes, hiRes)))
      post = urllib2.urlopen(req)
      tempList = [x.strip() for x in post]
      tempList = self.parseXmlQuery(tempList)
      seleNames |= set(tempList)
      return(seleNames)
    
    elif selPars[4] == "Y":
        #Get crystal structure names, if any
      req = urllib2.Request(self.server2,
                            (self.xml % (selPars[2], selPars[0], selPars[1], selPars[3], loRes, hiRes)))
      post = urllib2.urlopen(req)
      tempList = [x.strip() for x in post]
      tempList = self.parseXmlQuery(tempList)
      seleNames |= set(tempList)

      req = urllib2.Request(self.server2,
                            (self.xmlNmrNoModels % (selPars[2], selPars[0], selPars[1], selPars[3])))
      post = urllib2.urlopen(req)
      tempList = [x.strip() for x in post]
      tempList = self.parseXmlQuery(tempList)
      seleNamesNMR |= set(tempList)
      seleNames = seleNames | seleNamesNMR
      return(seleNames)

  def searchStrucName(self, inputName):
    outNames = set([])
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.MoleculeNameQuery</queryType>
      <macromoleculeName>%s</macromoleculeName>
    </orgPdbQuery>''' % inputName

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def searchModResi(self, inpResi):
    outNames = set([])
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.NoModResQuery</queryType>
        <hasModifiedResidues>yes</hasModifiedResidues>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ChemCompNameQuery</queryType>
        <name>%s</name>
        <polymericType>Any</polymericType>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>''' % inpResi

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def searchLigands(self, ligandName):
    outNames = set([])
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbCompositeQuery version="1.0">
     <queryRefinement>
      <queryRefinementLevel>0</queryRefinementLevel>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.NoLigandQuery</queryType>
        <haveLigands>yes</haveLigands>
      </orgPdbQuery>
     </queryRefinement>
     <queryRefinement>
      <queryRefinementLevel>1</queryRefinementLevel>
      <conjunctionType>and</conjunctionType>
      <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ChemCompNameQuery</queryType>
        <comparator>Contains</comparator>
        <name>%s</name>
        <polymericType>Any</polymericType>
      </orgPdbQuery>
     </queryRefinement>
    </orgPdbCompositeQuery>''' % ligandName

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def generalSearch(self, searchGen):
    outNames = set([])
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
      <keywords>%s</keywords>
    </orgPdbQuery>''' % searchGen

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def tempCrystal(self, tLo, tHi):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.CrystalQuery</queryType>
      <exptl_crystal.density_Matthews.comparator>between</exptl_crystal.density_Matthews.comparator>
      <exptl_crystal.density_percent_sol.comparator>between</exptl_crystal.density_percent_sol.comparator>
      <exptl_crystal_grow.temp.comparator>between</exptl_crystal_grow.temp.comparator>
      <exptl_crystal_grow.temp.min>%s</exptl_crystal_grow.temp.min>
      <exptl_crystal_grow.temp.max>%s</exptl_crystal_grow.temp.max>
      <exptl_crystal_grow.pH.comparator>between</exptl_crystal_grow.pH.comparator>
    </orgPdbQuery>''' % (tLo, tHi)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def tempBeamline(self, tLo, tHi):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.XrayDiffrnSourceQuery</queryType>
      <diffrn_source.source.comparator>like</diffrn_source.source.comparator>
      <diffrn_source.pdbx_synchrotron_site.comparator>contains</diffrn_source.pdbx_synchrotron_site.comparator>
      <diffrn_source.pdbx_synchrotron_beamline.comparator>contains</diffrn_source.pdbx_synchrotron_beamline.comparator>
      <diffrn.ambient_temp.comparator>between</diffrn.ambient_temp.comparator>
      <diffrn.ambient_temp.min>%s</diffrn.ambient_temp.min>
      <diffrn.ambient_temp.max>%s</diffrn.ambient_temp.max>
    </orgPdbQuery>''' % (tLo, tHi)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def phCrystal(self, phLo, phHi):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.CrystalQuery</queryType>
      <exptl_crystal.density_Matthews.comparator>between</exptl_crystal.density_Matthews.comparator>
      <exptl_crystal.density_percent_sol.comparator>between</exptl_crystal.density_percent_sol.comparator>
      <exptl_crystal_grow.temp.comparator>between</exptl_crystal_grow.temp.comparator>
      <exptl_crystal_grow.pH.comparator>between</exptl_crystal_grow.pH.comparator>
      <exptl_crystal_grow.pH.min>%s</exptl_crystal_grow.pH.min>
      <exptl_crystal_grow.pH.max>%s</exptl_crystal_grow.pH.max>
    </orgPdbQuery>''' % (phLo, phHi)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def seqSearch(self, sequence):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.MotifQuery</queryType>
      <motif>%s</motif>
    </orgPdbQuery>''' % (sequence)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def searchChemID(self, chemical):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
      <version>head</version>
      <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
      <chemCompId>%s</chemCompId>
      <polymericType>Any</polymericType>
    </orgPdbQuery>''' % (chemical)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

      #pull PDB parameters for a specific PDB provided a given attribute
      #Example for 4HHB
        # status : CURRENT
        # nr_residues : 574
        # structureId : 4HHB
        # structure_authors : Fermi, G., Perutz, M.F.
        # pubmedId : 6726807
        # keywords : OXYGEN TRANSPORT
        # nr_atoms : 4779
        # nr_entities : 4
        # last_modification_date : 2011-07-13
        # replaces : 1HHB
        # citation_authors : Fermi, G., Perutz, M.F., Shaanan, B., Fourme, R.
        # title : THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 ANGSTROMS RESOLUTION
        # resolution : 1.74
        # release_date : 1984-07-17
        # expMethod : X-RAY DIFFRACTION
        # deposition_date : 1984-03-07
    def pullPDBParameter(self, pdbName, attribute):
      pdbDict = {}
      req = urllib2.Request(("http://www.rcsb.org/pdb/rest/describePDB?structureId=%s" %("4HHB")), None)
      post = urllib2.urlopen(req)
      pdbTree = ET.parse(post)
      root = pdbTree.getroot()
      for child in root:
        pdbDict = child.attrib

      return(pdbDict[attribute])

  def pruneRibosome(self):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.MoleculeNameQuery</queryType>
    <macromoleculeName>ribosomal</macromoleculeName>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.MolecularWeightQuery</queryType>
    <mvStructure.structureMolecularWeight.min>150000.0</mvStructure.structureMolecularWeight.min>
    <mvStructure.structureMolecularWeight.max>1.0E8</mvStructure.structureMolecularWeight.max>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>2</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <containsProtein>?</containsProtein>
    <containsDna>?</containsDna>
    <containsRna>Y</containsRna>
    <containsHybrid>?</containsHybrid>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>''' 

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def matchWeight(self, loWt, hiWt):
    outNames = set([])
    #just a temporary list, nothing to do with temperatures
    tempList = []
    searchParams = '''<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>
  <version>head</version>
  <queryType>org.pdb.query.simple.MolecularWeightQuery</queryType>
  <mvStructure.structureMolecularWeight.min>%s</mvStructure.structureMolecularWeight.min>
  <mvStructure.structureMolecularWeight.max>%s</mvStructure.structureMolecularWeight.max>
</orgPdbQuery>''' % (loWt, hiWt)

    req = urllib2.Request(self.server2, searchParams)
    post = urllib2.urlopen(req)
    tempList = [x.strip() for x in post]
    tempList = self.parseXmlQuery(tempList)
    outNames |= set(tempList)
    return(outNames)  

  def pPdbData(self, fileName, pdbCobj):
    req = urllib2.Request("http://www.rcsb.org/pdb/rest/customReport?pdbids=%s&customReportColumns=resolution,crystallizationTempK,phValue,pdbxDetails,pubmedId" % fileName)
    post = urllib2.urlopen(req)
    root = ET.fromstring(post.read())
      #For some reason some PDB's have no root, program will crash if don't ignore these
      #by excluding those with root's of 0 length
    if len(root) != 0:
      pdbCobj.pdbId = root[0][0].text
      pdbCobj.resolution = root[0][1].text
      pdbCobj.temp = root[0][2].text
      pdbCobj.phVal = root[0][3].text
      pdbCobj.conditions = root[0][4].text.replace(",", " |")
      pdbCobj.pubmedId = root[0][5].text
      # print (pdbCobj.pdbId, pdbCobj.resolution, pdbCobj.temp, 
                # pdbCobj.phVal, pdbCobj.pubmedId)
      return(pdbCobj)

    #Returns 2 data sets
  def matchSmiles(self, smiles):
    ligands = set([])
    smilelig = {}
    iddata = []
    try:
      req = urllib2.Request(
            "http://www.rcsb.org/pdb/rest/smilesQuery?smiles=%s&similarity=0.7&search_type=1"
            % smiles)
      post = urllib2.urlopen(req)
      for i in post:
        if "chemicalID" in i:
          iddata.append(i.strip())
        #Grab the ligand ID from the XML query result
      for i in iddata:
        ligands.add(i.split()[2].replace("\"","").split("=")[1])
      smilelig[smiles] = ligands
      return(ligands, smilelig)
      #Handle URL exceptions
    except urllib2.HTTPError, e:
      print "Post failed for SMILES", smiles
      print " Error code %s" % e.code
      return(None, None)


    #  #This is for counting ligands in the ID data
    # for i in ligands:
    #   count = 0
    #   for y in iddata:
    #     if i in y.split()[2].replace("\"","").split("=")[1]:
    #       count += 1
    #   print "Ligand,", i, "has %s hits" % count
