import getopt, sys, os, shutil, urllib, gc #SOAPpy
from sets import Set
import strucs

class ParserC:
  """Parse data"""
  
  def __init__(self):
    self.allpairs = set([])
    self.nopairs = set([])
    self.badpairs = set([])

  def fileExists(self, filePath):
    try:
      with open(filePath) as f: return(True)
    except IOError as e: return(False)

    #FileNames are those that are stated to have base pairs
    #As was dictated by the .outp files
    #Pairs is all the pair data already deposited by the parseOutp function
  def parseDSSR(self, dssrDir, fileNames, pairs, hairpins): 
      #Collects the number of base pairs dictated in each .dssr file
    for ii in fileNames:
        #Raw data for each file
      rD = []
      tnumps, ct = 0, -1
        #mnumps - number of multiplets
        #numhpl - number of hairpin loops
        #numil - number of internal loops
      mnumps, numhpl, numil = 0, 0, 0
        #Pair ID range
      pidR = []
        #Multiplet Range
        #hplR - hairpin loop Range
      multR, hplR = [], []
      try:
        with open(dssrDir + ii + ".dssr", "r") as f:
          rD = [line for line in f]
      except IOError as e: print "No file found for:", ii + ".dssr"
      for x in range(len(rD)):
        if "List of" in rD[x] and "base pair(s)" in rD[x]:
          tnumps = int(rD[x].strip().split()[2])
          # print "PDB:", ii, "contains (", tnumps, ") base pairs."
          pidR = range(x+1, x+1+(5*tnumps), 5)

        if "List of" in rD[x] and "multiplet(s)" in rD[x]:
          mnumps = int(rD[x].strip().split()[2])
          # print "  and (", mnumps, ") multiplets."
          multR = range(x+1, x+1+(mnumps))

        if "List of" in rD[x] and "hairpin loop(s)" in rD[x]:
          numhpl = int(rD[x].strip().split()[2])
          hplR = range(x+1, (x+1+(2*numhpl)), 2)

        #Get Leontis-Westhof nomenclature for bp
      for jj in pidR:
        temp = rD[jj].strip().split()
        # print temp
        combo = (ii + "_" + temp[1].replace("/","").replace("^","") + "_" 
                        + temp[2].replace("/","").replace("^", ""))
        combo = combo.replace(".","")
        for yy in pairs:
          if combo in yy:
            pairs[yy].lwN = temp[-1]

        # Parse multiplets
      for jj in multR:
        temp = rD[jj].strip().split()
        tempPlist = temp[2].split("+")
        tempPlist = [x.replace(".", "") for x in tempPlist]

        for kk in tempPlist:
          for nn in pairs:
            if (kk == (pairs[nn].bAstrand + pairs[nn].bAsubtype + str(pairs[nn].bAnum)
                + pairs[nn].bAnumX) and pairs[nn].pName == ii):
              pairs[nn].multbases = temp[2].strip()
              pairs[nn].multType = temp[-1].replace("[","").replace("]", "")
              pairs[nn].multNum = len(tempPlist)

            elif (kk == (pairs[nn].bBstrand + pairs[nn].bBsubtype + str(pairs[nn].bBnum)
                + pairs[nn].bBnumX) and pairs[nn].pName == ii):      
              pairs[nn].multbases = temp[2].strip()
              pairs[nn].multType = temp[-1].replace("[","").replace("]", "")
              pairs[nn].multNum = len(tempPlist)

      hairpinResi = []
      hairpinType = []
      hairpinSize = []

      if len(hplR) != 0:
        for jj in hplR:
          temp = rD[jj + 1].strip().split()
          tempPlist = temp[0].split("+")
          tempPlist = [x[1:-1].replace(".","_") for x in tempPlist]
          hpType = temp[1].replace("[","").replace("]","")
          hpSize = len(hpType)
          hairpinResi.append(tempPlist)
          hairpinType.append(hpType)
          hairpinSize.append(hpSize-2)
      else:
        hairpinResi, hairpinType, hairpinSize = None, None, None
      if hairpinSize != None:
        print ii, len(hairpinSize)
      else: print ii

      #This function parses the outp files it is given from
      #the passed variable fileNames, in the folder "outpDir",
      #and stores the associated data in the "pairs" list
  def parseOutp(self, outpDir, fileNames, pairs):
    rawData = []
    for i in fileNames: 
        #Temp bp class object to be appended in to primary bp list
      numPairs = 0
      numBonds = 0
      ct = -1
      sp = []
      boolSugar = False
      boolBB = False
      try:
        with open(outpDir + i + ".outp", "r") as f:
          # count += 1
          # sys.stdout.write("\r" + "  [Parsing %d PDB files]" % count)
          # sys.stdout.flush()
          rawData = [line for line in f]
          for x in range(len(rawData)):      
            if rawData[x].find("Number of base-pairs:") >= 0:
              numPairs = int(rawData[x][22:].strip())
            if rawData[x].find("Strand I                    Strand II") >= 0:
              for bpN in range((x+1), (x+numPairs+1)):
                tC = strucs.BasepairC()
                sp.append(tC)
                ct = ct + 1
                sp[ct].pName = i
                sp[ct].pNumPairs = numPairs
                  #Strand name/num
                sp[ct].bAstrand = rawData[bpN][20].strip()
                sp[ct].bBstrand = rawData[bpN][52].strip()
                  #Base num
                sp[ct].bAnum = int(rawData[bpN][22:26].replace(".", "").strip())
                sp[ct].bBnum = int(rawData[bpN][46:50].replace(".", "").strip())
                  #Base num addition
                if rawData[bpN][26] != "_":
                  sp[ct].bAnumX = rawData[bpN][26]
                if rawData[bpN][50] != "_":
                  sp[ct].bBnumX = rawData[bpN][50]
                  #Base type
                sp[ct].bAtype = rawData[bpN][33]
                sp[ct].bBtype = rawData[bpN][39]
                  #Base sub-type
                sp[ct].bAsubtype = rawData[bpN][29:32].replace(".", "")
                sp[ct].bBsubtype = rawData[bpN][41:44].replace(".", "")
                  #Bonding notation, ie WC, etc
                sp[ct].bNotn = rawData[bpN][34:39]

                if sp[ct].bNotn != "-----": sp[ct].isWC = False
                else: 
                  sp[ct].isWC = True

            ct = -1
            if rawData[x].find("Detailed H-bond") >= 0:
              for bpN in range((x+1), (x+numPairs+1)):
                ct = ct + 1
                  #Number of bonding atoms in this pair
                sp[ct].numBonds = int(rawData[bpN][15])
                  #Temp fix for 3R8T_A (G 2857) -- (G 2859) A
                if sp[ct].numBonds == 5:
                  sp[ct].numBonds = 4
                if sp[ct].numBonds == 1:
                  sp[ct].bondAtLHS.append(rawData[bpN][19:22].strip())
                  sp[ct].bondType.append(rawData[bpN][22])
                  sp[ct].bondAtRHS.append(rawData[bpN][24:27].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][28:32].strip()))
                  for i in range(3):
                    sp[ct].bondAtLHS.append("---")
                    sp[ct].bondType.append("---")
                    sp[ct].bondAtRHS.append("---")
                    sp[ct].bondDists.append(0.0)                    
                  if sp[ct].bondAtLHS[0] < sp[ct].bondAtRHS[0]:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtLHS[0] + sp[ct].bondAtRHS[0]
                  else:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtRHS[0] + sp[ct].bondAtLHS[0]               
                elif sp[ct].numBonds == 2:
                  sp[ct].bondAtLHS.append(rawData[bpN][19:22].strip())
                  sp[ct].bondType.append(rawData[bpN][22])
                  sp[ct].bondAtRHS.append(rawData[bpN][24:27].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][28:32].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][34:37].strip())
                  sp[ct].bondType.append(rawData[bpN][37])
                  sp[ct].bondAtRHS.append(rawData[bpN][39:42].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][43:47].strip()))
                  for i in range(2):
                    sp[ct].bondAtLHS.append("---")
                    sp[ct].bondType.append("---")
                    sp[ct].bondAtRHS.append("---")
                    sp[ct].bondDists.append(0.0)                    
                  if sp[ct].bondAtLHS[0] < sp[ct].bondAtRHS[0]:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtLHS[0] + sp[ct].bondAtRHS[0]
                  else:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtRHS[0] + sp[ct].bondAtLHS[0]
                  if sp[ct].bondAtLHS[1] < sp[ct].bondAtRHS[1]:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtLHS[1] + sp[ct].bondAtRHS[1]
                  else:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtRHS[1] + sp[ct].bondAtLHS[1]
                elif sp[ct].numBonds == 3:
                  # if sp[ct].pName == "1S72":
                  #   print "---------",bpN,"------------"
                  #   print "First H-Bond:", rawData[bpN][19:22].strip(), rawData[bpN][22], rawData[bpN][24:27].strip(), rawData[bpN][28:32].strip()
                  #   print "Second H-Bond:", rawData[bpN][34:37].strip(), rawData[bpN][37], rawData[bpN][39:42].strip(), rawData[bpN][43:47].strip()
                  #   print "Third H-Bond:", rawData[bpN][49:52].strip(), rawData[bpN][52], rawData[bpN][54:57].strip(), rawData[bpN][58:62].strip()
                  sp[ct].bondAtLHS.append(rawData[bpN][19:22].strip())
                  sp[ct].bondType.append(rawData[bpN][22])
                  sp[ct].bondAtRHS.append(rawData[bpN][24:27].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][28:32].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][34:37].strip())
                  sp[ct].bondType.append(rawData[bpN][37])
                  sp[ct].bondAtRHS.append(rawData[bpN][39:42].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][43:47].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][49:52].strip())
                  sp[ct].bondType.append(rawData[bpN][52])
                  sp[ct].bondAtRHS.append(rawData[bpN][54:57].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][58:62].strip()))
                  sp[ct].bondAtLHS.append("---")
                  sp[ct].bondType.append("---")
                  sp[ct].bondAtRHS.append("---")
                  sp[ct].bondDists.append(0.0)                    
                  if sp[ct].bondAtLHS[0] < sp[ct].bondAtRHS[0]:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtLHS[0] + sp[ct].bondAtRHS[0]
                  else:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtRHS[0] + sp[ct].bondAtLHS[0]
                  if sp[ct].bondAtLHS[1] < sp[ct].bondAtRHS[1]:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtLHS[1] + sp[ct].bondAtRHS[1]
                  else:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtRHS[1] + sp[ct].bondAtLHS[1]
                  if sp[ct].bondAtLHS[2] < sp[ct].bondAtRHS[2]:
                    sp[ct].bondAtCombo[2] = sp[ct].bondAtLHS[2] + sp[ct].bondAtRHS[2]
                  else:
                    sp[ct].bondAtCombo[2] = sp[ct].bondAtRHS[2] + sp[ct].bondAtLHS[2]
                elif sp[ct].numBonds == 4:
                  sp[ct].bondAtLHS.append(rawData[bpN][19:22].strip())
                  sp[ct].bondType.append(rawData[bpN][22])
                  sp[ct].bondAtRHS.append(rawData[bpN][24:27].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][28:32].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][34:37].strip())
                  sp[ct].bondType.append(rawData[bpN][37])
                  sp[ct].bondAtRHS.append(rawData[bpN][39:42].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][43:47].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][49:52].strip())
                  sp[ct].bondType.append(rawData[bpN][52])
                  sp[ct].bondAtRHS.append(rawData[bpN][54:57].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][58:62].strip()))
                  sp[ct].bondAtLHS.append(rawData[bpN][64:67].strip())
                  sp[ct].bondType.append(rawData[bpN][67])
                  sp[ct].bondAtRHS.append(rawData[bpN][69:72].strip())
                  sp[ct].bondDists.append(float(rawData[bpN][73:77].strip()))
                  if sp[ct].bondAtLHS[0] < sp[ct].bondAtRHS[0]:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtLHS[0] + sp[ct].bondAtRHS[0]
                  else:
                    sp[ct].bondAtCombo[0] = sp[ct].bondAtRHS[0] + sp[ct].bondAtLHS[0]
                  if sp[ct].bondAtLHS[1] < sp[ct].bondAtRHS[1]:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtLHS[1] + sp[ct].bondAtRHS[1]
                  else:
                    sp[ct].bondAtCombo[1] = sp[ct].bondAtRHS[1] + sp[ct].bondAtLHS[1]
                  if sp[ct].bondAtLHS[2] < sp[ct].bondAtRHS[2]:
                    sp[ct].bondAtCombo[2] = sp[ct].bondAtLHS[2] + sp[ct].bondAtRHS[2]
                  else:
                    sp[ct].bondAtCombo[2] = sp[ct].bondAtRHS[2] + sp[ct].bondAtLHS[2]
                  if sp[ct].bondAtLHS[3] < sp[ct].bondAtRHS[3]:
                    sp[ct].bondAtCombo[3] = sp[ct].bondAtLHS[3] + sp[ct].bondAtRHS[3]
                  else:
                    sp[ct].bondAtCombo[3] = sp[ct].bondAtRHS[3] + sp[ct].bondAtLHS[3]

            ct = -1
            if rawData[x].find("Local base-pair parameters") >= 0:              
              for bpN in range((x+2), (x+numPairs+2)):
                ct = ct + 1
                sp[ct].shear = float(rawData[bpN][13:20].strip())
                sp[ct].stretch = float(rawData[bpN][23:30].strip())
                sp[ct].stagger = float(rawData[bpN][33:40].strip())
                sp[ct].buckle = float(rawData[bpN][43:50].strip())
                sp[ct].propeller = float(rawData[bpN][53:60].strip())
                sp[ct].opening = float(rawData[bpN][63:70].strip())
                ##Local base pair step parameters.
            # ct = -1
            # if rawData[x].find("Local base-pair step parameters") >= 0:
            #   for bpN in range((x+2), (x+numPairs+1)):
            #     ct = ct + 1
            #     sp[ct].shift = float(rawData[bpN][13:20].strip())
            #     sp[ct].slide = float(rawData[bpN][23:30].strip())
            #     sp[ct].rise = float(rawData[bpN][33:40].strip())
            #     sp[ct].tilt = float(rawData[bpN][43:50].strip())
            #     sp[ct].roll = float(rawData[bpN][53:60].strip())
            #     sp[ct].twist = float(rawData[bpN][63:70].strip())
            ##DONT REMOVE ME
            # ct = -1
            # if rawData[x].find("Local base-pair helical parameters") >= 0:
            #   for bpN in range((x+2), (x+numPairs+1)):
            #     ct = ct + 1
            #     sp[ct].xDisp = float(rawData[bpN][13:20].strip())
            #     sp[ct].yDisp = float(rawData[bpN][23:30].strip())
            #     sp[ct].hRise = float(rawData[bpN][33:40].strip())
            #     sp[ct].hIncl = float(rawData[bpN][43:50].strip())
            #     sp[ct].hTip = float(rawData[bpN][53:60].strip())
            #     sp[ct].hTwist = float(rawData[bpN][63:70].strip())
            ##DONT REMOVE ME
            ct = -1
            if rawData[x].find("lambda: virtual angle between C1") >= 0:
              for bpN in range((x+8), (x+numPairs+8)):
                ct = ct + 1
                ##DONT REMOVE ME
                # sp[ct].lamda1 = rawData[bpN][13:19].strip()
                # sp[ct].lamda2 = rawData[bpN][23:29].strip()
                ##DONT REMOVE ME
                sp[ct].c1c1 = rawData[bpN][33:39].strip()
                sp[ct].rn9yn1 = rawData[bpN][43:49].strip()
            ct = -1
            if rawData[x].find("Main chain and chi torsion angles") >= 0:
              boolBB = True
              for bpN in range((x+14), (x+numPairs+14)):
                ct = ct + 1
                  #Strand A
                sp[ct].bbAlphaA = rawData[bpN][9:15].strip()
                sp[ct].bbBetaA = rawData[bpN][17:23].strip()
                sp[ct].bbGammaA = rawData[bpN][25:31].strip()
                sp[ct].bbDeltaA = rawData[bpN][33:39].strip()
                sp[ct].bbEpsilonA = rawData[bpN][41:47].strip()
                sp[ct].bbZetaA = rawData[bpN][49:55].strip()
                sp[ct].bbChiA = rawData[bpN][57:63].strip()            
                  #Strand B
                sp[ct].bbAlphaB = rawData[bpN+numPairs+3][9:15].strip()
                sp[ct].bbBetaB = rawData[bpN+numPairs+3][17:23].strip()
                sp[ct].bbGammaB = rawData[bpN+numPairs+3][25:31].strip()
                sp[ct].bbDeltaB = rawData[bpN+numPairs+3][33:39].strip()
                sp[ct].bbEpsilonB = rawData[bpN+numPairs+3][41:47].strip()
                sp[ct].bbZetaB = rawData[bpN+numPairs+3][49:55].strip()
                sp[ct].bbChiB = rawData[bpN+numPairs+3][57:63].strip()    
            ct = -1
            if rawData[x].find("P:  the phase angle of pseudorotation") >= 0:
              boolSugar = True
              for bpN in range((x+4), (x+numPairs+4)):
                ct = ct + 1
                  #Strand A
                sp[ct].vee0a = rawData[bpN][9:15].strip()
                sp[ct].vee1a = rawData[bpN][17:23].strip()
                sp[ct].vee2a = rawData[bpN][25:31].strip()
                sp[ct].vee3a = rawData[bpN][33:39].strip()
                sp[ct].vee4a = rawData[bpN][41:47].strip()
                ##DONT REMOVE ME
                # sp[ct].tAmpa = rawData[bpN][49:55].strip()
                # sp[ct].pPhaseA = rawData[bpN][57:63].strip()
                sp[ct].puckerA = rawData[bpN][67:75].strip()
                  #Strand B
                sp[ct].vee0b = rawData[bpN+numPairs+3][9:15].strip()
                sp[ct].vee1b = rawData[bpN+numPairs+3][17:23].strip()
                sp[ct].vee2b = rawData[bpN+numPairs+3][25:31].strip()
                sp[ct].vee3b = rawData[bpN+numPairs+3][33:39].strip()
                sp[ct].vee4b = rawData[bpN+numPairs+3][41:47].strip()
                ##DONT REMOVE ME
                # sp[ct].tAmpb = rawData[bpN+numPairs+3][49:55].strip()
                # sp[ct].pPhaseB = rawData[bpN+numPairs+3][57:63].strip()
                sp[ct].puckerB = rawData[bpN+numPairs+3][67:75].strip()
          del rawData[:len(rawData)]
      except IOError as e: self.nopairs.add(i)
      # except IOError as e: print "    No base pairs found for ", i
      if (boolBB == True and boolSugar == True):
        [pairs.append(x) for x in sp]       
      else: self.badpairs.add(i)      
      boolBB = False 
      boolSugar = False
    self.badpairs = self.badpairs.difference(self.nopairs)
    print "  [Parsed %d PDBs]" % (len(fileNames) - len(self.nopairs))
