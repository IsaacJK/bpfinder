import os, sys, string, csv

def sysHelp():
  print "Renames appropriately titled structure files with their cluster ID"
  print "Usage:"
  print " >clustStruc -[mv/sub] [in .CSV file] [file folder]"

def writeClusts(filename, outname, outFolder, Rfilename, numclusts):
  FILE = open(Rfilename + ".r", "wb")
  FILE.write('''
#Must run 64-bit R
#>R --arch x86_64

#Change directory to folder working in
setwd("%s")

library(mclust)
myData = read.csv("%s.csv", sep=",", header=FALSE)

#Scale data matrix / standardize & center
#Use 13:18 for all, 13:15 for shear stretch stagger, 13:18 for buckle/propeller/open
myStdData = scale(myData[,5:7])

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
            outname, outname, outname, outFolder))

argc = len(sys.argv)
curDir = os.getcwd() + "/"

if argc != 4:
  sysHelp()

elif sys.argv[1] == "-mv":
  if curDir in sys.argv[2]:
    filename = sys.argv[2]
  else:
    filename = curDir + sys.argv[2]

  if curDir in sys.argv[3]:
    folder = sys.argv[3]
  else:
    folder = curDir + sys.argv[3]

  FILE = open(filename, "rU")
  f = csv.reader(FILE, delimiter=',')
  f.next()

  data = []
  clustdata = {}
  for i in f:
    data.append(i)

  numclust = set([])
  clusters = {}
  for i in data:
    clusters[i[3]] = i[70]
    numclust.add(i[70])

  for i in clusters:
  #   # print "mv %s %s" % ((folder + i + ".pse"), 
  #     #                               (folder + clusters[i] + "_" + i + ".pse"))
  #   # os.system("mv %s %s" % ((folder + i + ".pse"), 
  #   #                   (folder + clusters[i] + "_" + i + ".pse")))

    os.system("mv %s %s" % ((folder + i + ".pse"), 
                      (folder + clusters[i] + "_" + i + ".pse")))

elif sys.argv[1] == "-sub":
  if curDir in sys.argv[2]:
    filename = sys.argv[2]
  else:
    filename = curDir + sys.argv[2]

  if curDir in sys.argv[3]:
    folder = sys.argv[3]
  else:
    folder = curDir + sys.argv[3]

  FILE = open(filename, "rU")
  f = csv.reader(FILE, delimiter=',')
  f.next()

  data = []
  clustdata = {}
  for i in f:
    data.append(i)

  numclust = set([])
  clusters = {}
  for i in data:
    clusters[i[3]] = i[70]
    numclust.add(i[70])

  sortclust = []

  for i in range(len(numclust)):
    sortclust.append(i + 1)

  for y in sortclust:
    filename = ("clust_%s" % str(y))
    FILE = open((filename + ".csv"), "wb")
    for i in data:
      if y == int(i[70]):
        FILE.write(i[3] + ","  + i[12] + "," + i[13] + "," + i[14] + "," + i[15] + "," + i[16] + "," + i[17]
                         + "," + str(y) + "," + "\n")
    FILE.close
    writeClusts(filename, 
                        filename,
                         curDir, 
                         filename,
                         20)

  for y in range(len(sortclust) + 1):
      print "NOW WORKING ON: clust_", y+1
      os.system("R --arch x86_64 --silent < clust_%s.r --no-save > clustOut_%s.txt" 
                        % (str(y+1), str(y+1)))
      print "**************"
