class pdbDataC:
  '''PDB Data like temp, pH, etc'''
  def __init__(self):
    self.phVal, self.temp, self.resolution = 0.0, 0.0, 0.0
    self.pdbId = "---"
    self.pubmedId = "---"
    self.conditions = "---"

class hairpinC:
  def __init__(self):
    self.pName = None
    self.uniqueID, self.size, self.idNum = None, None, None
    self.sType, self.lType = None, None


class pdbClass:
    '''PDB Data - Hairpins, loops, etc.'''
    def __init__(self):
        self.pdbName = None
          #Hairpin Loops
        self.hplSize, self.hplResi, self.hplNames = None, None, None
          #Internal Loops
        self.ilSize, self.ilMx = None, None
        self.ilResi, self.ilNames = None, None
          #Bulges
        self.blgSize, self.blgMx = None, None
        self.blgResi, self.blgNames = None, None

class BasepairC:
  """Basepair Data"""

  def __init__(self):
      ###Descriptors###

      #pName - PDB name (EX. 1S72)
      #pNumPairs - number of base pairs as dictated by the .outp file
    self.pName, self.pNumPairs = None, 0

      #bAtype/bBtype - Base truncated type (ie. G for G and DG)
      # lowercase values indicate a modified version (ie. c for CBR)
    self.bAtype, self.bBtype = None, None
      #bAsubtype/bBsubtype - full length name of base type
      # IE. CBR instead of c for brominated cytidine
    self.bAsubtype, self.bBsubtype = None, None
    # Dictionary of base types and subtypes
    self.btyped = {"bAtype": self.bAtype, "bAsubtype": self.bAsubtype,
                   "bBtype": self.bBtype, "bBsubtype": self.bBsubtype}
      #bAnum/bBnum - base pair number. IE 33 for G33
      #bNotn - Pair notation dictated by 3DNA (- in G-C, etc)
    self.bAnum, self.bBnum, self.bNotn = 0, 0, ""
      #bAnumX/bBnumX - addendum to base pair number.
      #  IE. "L" in G33L
    self.bAnumX, self.bBnumX = "", ""
      #bAstrand/bBstrand - the strand name of base A/B
    self.bAstrand, self.bBstrand = None, None
      #isWC - Boolean to descr. whether pair is watson-crick
    self.isWC = True
      #bondAtLHS/RHS - a list of h-bonding atoms in base A and B
    self.bondAtLHS, self.bondAtRHS = [], []
      #bondAtCombo - bonding atom combinations of a pair
      #  in alphabetical order. IE N1N17, etc
    self.bondAtCombo = ["---", "---", "---", "---"]
      #bondType - a list of bonding types det. by 3DNA
      #  IE. * in N1*N7
    self.bondType = []
      #bondDists - A list containing all h-bond distances
      #numBonds - How many HBonds this pair makes
    self.bondDists, self.numBonds = [], 0
      #puckerA/B - text descriptor of sugar pucker for base A/B
    self.puckerA, self.puckerB = None, None
      #uniqueID - a unique descriptor for this base pair
      #  widely used as an identifier. Composed of the following:
      #  PDB_[pair motif, alphabetical]_chainA_typeA_Anum_AnumX
      #  _chainB_typeB_Bnum_BnumX
      #pairMotif - alphabetical base pair motif
      #  IE. CBR - G would be cG
      #subPMotif - alphabetical base subtype motif
      # IE. CBR - G would be CBRG
    self.uniqueID, self.pairMotif, self.subPMotif = "", "", ""

      #Local BP params
    self.shear, self.stretch, self.stagger = 0.0, 0.0, 0.0
    self.buckle, self.propeller, self.opening = 0.0, 0.0, 0.0
        
    #   #Local BP step-params
    # self.shift, self.slide, self.rise  = 0.0, 0.0, 0.0
    # self.tilt, self.roll, self.twist = 0.0, 0.0, 0.0

      #DSSR Parameters
      #lwN - Leontis-Westhoff descriptors of a base
      #  IE. cW-W is a cis Watson-Watson facing bases
      #multbases - List of base ID's involved in multiplet
      #multNum - Multiplet Number of base pair
      #  IE. 2 is doublet, 3 is triplet, etc.
      #multType - multiplet type (IE. UUG, CCU, etc.)
    self.lwN, self.multbases, self.multNum, self.multType = "---", "---", 2, "---"

      #Virtual angles
    # self.lamda1, self.lamda2 = 0.0, 0.0
    self.c1c1, self.rn9yn1 = 0.0, 0.0
        #Temp removed this to save space
        #If adding them back, make sure to uncomment parse.py
    #   #Local helical params
    # self.xDisp, self.yDisp, self.hRise = 0.0, 0.0, 0.0
    # self.hIncl, self.hTip, self.hTwist = 0.0, 0.0, 0.0

      #Backbone torsions
    self.bbAlphaA, self.bbBetaA, self.bbGammaA = 0.0, 0.0, 0.0
    self.bbDeltaA, self.bbEpsilonA, self.bbZetaA = 0.0, 0.0, 0.0
    self.bbChiA, self.bbChiAconf = 0.0, "---"
    self.bbAlphaB, self.bbBetaB, self.bbGammaB = 0.0, 0.0, 0.0
    self.bbDeltaB, self.bbEpsilonB, self.bbZetaB = 0.0, 0.0, 0.0
    self.bbChiB, self.bbChiBconf = 0.0, "---"

      #Sugar conformational params
    self.vee0a, self.vee0b, self.vee1a, self.vee1b = 0.0, 0.0, 0.0, 0.0
    self.vee2a, self.vee2b, self.vee3a, self.vee3b = 0.0, 0.0, 0.0, 0.0
    self.vee4a, self.vee4b = 0.0, 0.0
    # self.tAmpa, self.tAmpb, self.pPhaseA, self.pPhaseB = 0.0, 0.0, 0.0, 0.0

    #   #Shift mean and standard deviations
    # self.C1p, self.C2, self.C2p, self.C3p = 0.0, 0.0, 0.0, 0.0
    # self.C4, self.C4p, self.C5, self.C5p = 0.0, 0.0, 0.0, 0.0
    # self.C6, self.C8, self.H1, self.H1p = 0.0, 0.0, 0.0, 0.0
    # self.H2, self.H2p, self.H21, self.H22 = 0.0, 0.0, 0.0, 0.0
    # self.H3, self.H3p, self.H4p, self.H41 = 0.0, 0.0, 0.0, 0.0
    # self.H42, self.H5, self.H5p, self.H5pp = 0.0, 0.0, 0.0, 0.0
    # self.H6, self.H61, self.H62, self.H8 = 0.0, 0.0, 0.0, 0.0
    # self.HO2p, self.N1, self.N2, self.N4 = 0.0, 0.0, 0.0, 0.0
    # self.N6, self.N7, self.N9, self.P = 0.0, 0.0, 0.0, 0.0
