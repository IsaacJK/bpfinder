def writeProgram(filename, outFolder, expansion, filetype, realpath, workingpath):
  FILE = open("tempBpfPym.py", "wb")
  FILE.write('''
import colorsys,sys, csv, os
from pymol import cmd

def testInt(incVal):
  try:
    int(incVal)
    return True
  except ValueError:
    return False

def exportFiles(filename, outFolder, expansion, filetype):
  array = []
  count = -1
  FILE = open(filename, "rU")
  f = csv.reader(FILE, delimiter=',')
  pdbname, chain1, type1, resi1, chain2, type2, resi2 = "", "", "", "", "", "", ""
  for i in f:
    count += 1
    if count != 0:
      pdbname = i[0]
      ident = i[3]
      chain1 = i[4]
      type1 = i[8]
      resi1 = i[10]
      chain2 = i[5]
      type2 = i[9]
      resi2 = i[11]
      string = "string"

      cmd.fetch(pdbname)
      cmd.util.cbag('all')

      if testInt(resi1) == True and testInt(resi2) == True:
        if int(resi1) < 0 and int(resi2) >= 0:
          strSel = ("(chain %s and resi \%s)+(chain %s and resi %s)"
                    % (chain1, resi1, chain2, resi2))        
        elif int(resi2) < 0 and int(resi1) >= 0:
          strSel = ("(chain %s and resi %s)+(chain %s and resi \%s)"
                    % (chain1, resi1, chain2, resi2))   
        elif int(resi2) < 0 and int(resi1) < 0:
          strSel = ("(chain %s and resi \%s)+(chain %s and resi \%s)"
                    % (chain1, resi1, chain2, resi2))   
        else:
          strSel = ("(chain %s and resi %s)+(chain %s and resi %s)"
                    % (chain1, resi1, chain2, resi2))   
      else:
        strSel = ("(chain %s and resi %s)+(chain %s and resi %s)"
                  % (chain1, resi1, chain2, resi2))     
      cmd.select(strSel)
      cmd.center("sele")
      cmd.util.cbac('sele')
      strObj1 = ("%s_%s_%s-%s_%s"
                % (pdbname, chain1, resi1, chain2, resi2))
      cmd.create(strObj1, "byres sele expand %s" % expansion)
      cmd.delete(pdbname)

      strObj2 = ("%s_%s (%s %s) -- (%s %s) %s"
                 % (pdbname, chain1, type1, resi1, type2, resi2, chain2))
      if filetype == "PSE":
        strSave = outFolder + ident + ".pse"
        cmd.save(strSave, strObj1, 0, 'pse')
      elif filetype == "PDB":
        strSave = outFolder + ident + ".pdb"
        cmd.save(strSave, strObj1, 0, 'pdb')
      cmd.delete(strObj1)
    ''')

  FILE.write('''
if not os.path.exists("%s"): os.makedirs("%s")
os.chdir("%s")
    ''' % (realpath, realpath, realpath))

  FILE.write('''
exportFiles("%s", "%s", %s, "%s")
    ''' % (filename, outFolder, expansion, filetype))

  FILE.write('''
os.chdir("%s")
    ''' % (workingpath))

  FILE.close

def writeClusts(filename, outname, outFolder, workingpath, numclusts):
  FILE = open("tempClust.r", "wb")
  FILE.write('''
#Must run 64-bit R
#>R --arch x86_64

#Change directory to folder working in
setwd("%s")

library(mclust)
myData = read.csv("%s.csv", sep=",", header=TRUE)

#Scale data matrix / standardize & center
#Use 13:18 for all, 13:15 for shear stretch stagger, 16:18 for buckle/propeller/open
myStdData = scale(myData[,13:18])

myClust = Mclust(myStdData, G=1:%s)
summary(myClust, paramters=TRUE)
pdf("%s_plots.pdf")
plot(myClust)
dev.off()
myData$EMClusters = myClust$classification
write.csv(myData, file="%s_clusters.csv", row.names=FALSE)
means = myClust$parameters$mean
write.csv(means, file="%s_means.csv")
variance = myClust$parameters$variance
write.csv(variance[4], file="%s_variance.csv")
#Change directory back again
setwd("%s")
    ''' % (outFolder, filename, numclusts, outname, 
            outname, outname, outname, workingpath))