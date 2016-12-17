import numpy as np
# from numpy import (mean, array, concatenate, column_stack, isnan
#                                     as np)
          
import matplotlib.pyplot as plt
import sys, csv, math
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.mlab import PCA
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D

class pairStats:
  def __init__ (self): None

  def plotStats(self, filename, hasShifts):
    pdf = PdfPages(filename.strip(".csv") + ".pdf")

    ## Put the variable names here.
    VN = ["Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Opening"]
    WN = ["Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta"]
    XN = ["V0", "V1", "V2", "V3", "V4"]
    X = np.genfromtxt((filename), delimiter=',', skip_header=1)
    pMotif, sugarP = [], []
    xvals, yvals, xlabels = [], [], []

    xLen = X.shape

    if len(xLen) > 1:
      if xLen[0] > 6:
        pairPars = np.array([[]])
        pairPars = X[:,12:18]
        pairPars = pairPars[~np.isnan(pairPars).any(1)]   

        bbTorsions = np.array([[]])
            #Concatenate lhs torsions and rhs torsions to one 2D array
        bbTorsions = np.concatenate(((X[:,34:40]),(X[:,40:46])), axis=0)

          #Not stacked
        bbTorsions2 = np.array([[]])
            #Concatenate lhs torsions and rhs torsions to one 2D array
        bbTorsions2 = X[:,34:46]

        chiAngles = np.array([[]])
             #Concatenate lhs torsions and rhs torsions to one 1D array
        chiAngles = np.concatenate(((X[:,30]),(X[:,31])), axis=0)

        sugar = np.array([[]])
            #Concatenate lhs torsions and rhs torsions to one 2D array
        sugar = np.concatenate(((X[:,20:25]),(X[:,25:30])), axis=0)

        sugar2 = np.array([[]])
            #Concatenate lhs torsions and rhs torsions to one 2D array
        sugar2 = X[:,20:30]

         #combo backbone and local pair parameters
        sugBB = np.array([[]])
        sugBB = np.column_stack((bbTorsions, sugar))

         #combo backbone and local pair parameters
        sugBB2 = np.array([[]])
        sugBB2 = np.column_stack((bbTorsions2, sugar2))

          #C1-C1 virtual distances
        c1c1 = np.array([])
        c1c1 = X[:,46]

        #Remove any row that has Nan values (ie "---")
        sugar = sugar[~np.isnan(sugar).any(1)]
        sugar2 = sugar2[~np.isnan(sugar2).any(1)]
        bbTorsions = bbTorsions[~np.isnan(bbTorsions).any(1)]
        bbTorsions2 = bbTorsions2[~np.isnan(bbTorsions2).any(1)]
        sugBB = sugBB[~np.isnan(sugBB).any(1)]
        sugBB2 = sugBB2[~np.isnan(sugBB2).any(1)]
        chiAngles = chiAngles[~np.isnan(chiAngles)]   
        c1c1 = c1c1[~np.isnan(c1c1)]

        # pParsBB = pParsBB[~np.isnan(pParsBB).any(1)]
        # for i in pParsBB:
        #     print i

        if hasShifts == 1:
            h6shifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            h6shifts = np.concatenate(((X[:,93]),(X[:,167])), axis=0)
                #Remove any row that has Nan values (ie "---")
            h6shifts = h6shifts[~np.isnan(h6shifts)]   

            h8shifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            h8shifts = np.concatenate(((X[:,96]),(X[:,170])), axis=0)
                #Remove any row that has Nan values (ie "---")
            h8shifts = h8shifts[~np.isnan(h8shifts)]   

            h1pshifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            h1pshifts = np.concatenate(((X[:,80]),(X[:,154])), axis=0)
                #Remove any row that has Nan values (ie "---")
            h1pshifts = h1pshifts[~np.isnan(h1pshifts)]   

            c8shifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            c8shifts = np.concatenate(((X[:,78]),(X[:,152])), axis=0)
                #Remove any row that has Nan values (ie "---")
            c8shifts = c8shifts[~np.isnan(c8shifts)]   
            c6shifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            c6shifts = np.concatenate(((X[:,77]),(X[:,151])), axis=0)
                #Remove any row that has Nan values (ie "---")
            c6shifts = c6shifts[~np.isnan(c6shifts)]   

            c8shifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            c8shifts = np.concatenate(((X[:,78]),(X[:,152])), axis=0)
                #Remove any row that has Nan values (ie "---")
            c8shifts = c8shifts[~np.isnan(c8shifts)]   

            c1pshifts = np.array([[]])
                 #Concatenate lhs torsions and rhs torsions to one 1D array
            c1pshifts = np.concatenate(((X[:,69]),(X[:,143])), axis=0)
                #Remove any row that has Nan values (ie "---")
            c1pshifts = c1pshifts[~np.isnan(c1pshifts)]   

      infile = open((filename), "rU")
      reader = csv.reader(infile, delimiter=",")

      counter = 1
      for i in reader:
        if counter != 1:
          pMotif.append(i[2])
          sugarP.append(i[18])
          sugarP.append(i[19])
        counter += 1
      infile.close

      cnt = Counter(pMotif)
      for i in cnt:
        xlabels.append(i)
        yvals.append(cnt[i])

      for i in range(len(xlabels)):
        xvals.append(i +1)
      plt.clf()
      plt.bar(xvals, yvals, align="center", color="#d10000")
      plt.xticks(xvals, xlabels, rotation="vertical")
      plt.ylabel("Frequency (n=%s)" % len(pMotif), size=17)
      plt.xlabel("Pair motif", size=17)
      pdf.savefig()

      xvals, yvals, xlabels = [], [], []

      cnt = Counter(sugarP)
      for i in cnt:
        xlabels.append(i)
        yvals.append(cnt[i])

      for i in range(len(xlabels)):
        xvals.append(i +1)
      plt.clf()
      plt.bar(xvals, yvals, align="center", color="#8400d0")
      plt.xticks(xvals, xlabels, rotation="vertical", size=8)
      plt.ylabel("Frequency (n=%s)" % len(sugarP), size=17)
      plt.xlabel("Pair motif", size=17)
      pdf.savefig()

      # Imagine the points belong to three groups.
      IG = np.floor(3*np.random.uniform(size=len(X)))
      IG = [np.flatnonzero(IG==j) for j in 0,1,2]

      def iqr(x):
          v = np.sort(x)
          return v[0.75*len(v)] - v[0.25*len(v)]

      if hasShifts == 1 and xLen[0] > 6:
              #Plot H1' shifts
          if len(h1pshifts) != 0:
              if 2*iqr(h1pshifts[:]/len(h1pshifts)**(1/3.0)) != 0.0:
                plt.clf()
                h = 2*iqr(h1pshifts[:])/len(h1pshifts)**(1/3.0)
                nbins = np.ceil((max(h1pshifts[:]) - min(h1pshifts[:]))/h)
                plt.hist(h1pshifts[:], color='#0330b4', alpha=0.7, bins=nbins)
                plt.xlabel(("H1' (ppm)"), size=17)
                plt.ylabel("Frequency (n=%s)" % len(h1pshifts), size=17)
                pdf.savefig()

              else:
                plt.clf()
                plt.bar(h1pshifts[0], len(h1pshifts), align="center", color="#0330b4")
                plt.ylabel("Frequency (n=%s)" % len(h1pshifts), size=17)
                plt.xlabel("H1' (ppm)", size=17)
                pdf.savefig()     

              #Plot H6 shifts
          if len(h6shifts) != 0:
              if 2*iqr(h6shifts[:]/len(h6shifts)**(1/3.0)) != 0.0:
                plt.clf()
                h = 2*iqr(h6shifts[:])/len(h6shifts)**(1/3.0)
                nbins = np.ceil((max(h6shifts[:]) - min(h6shifts[:]))/h)
                plt.hist(h6shifts[:], color='#0330b4', alpha=0.7, bins=nbins)
                plt.xlabel(("H6 (ppm)"), size=17)
                plt.ylabel("Frequency (n=%s)" % len(h6shifts), size=17)
                pdf.savefig()

              else:
                plt.clf()
                plt.bar(h6shifts[0], len(h6shifts), align="center", color="#0330b4")
                plt.ylabel("Frequency (n=%s)" % len(h6shifts), size=17)
                plt.xlabel("H6 (ppm)", size=17)
                pdf.savefig()                


          if len(h8shifts) != 0:
                #Plot of H8 shifts
              if 2*iqr(h8shifts[:]/len(h8shifts)**(1/3.0)) != 0.0:
                  plt.clf()
                  h = 2*iqr(h8shifts[:])/len(h8shifts)**(1/3.0)
                  nbins = np.ceil((max(h8shifts[:]) - min(h8shifts[:]))/h)
                  plt.hist(h8shifts[:], color='#0330b4', alpha=0.7, bins=nbins)
                  plt.xlabel(("H8 (ppm)"), size=17)
                  plt.ylabel("Frequency (n=%s)" % len(h8shifts), size=17)
                  pdf.savefig()

              else:
                plt.clf()
                plt.bar(h8shifts[0], len(h8shifts), align="center", color="#0330b4")
                plt.ylabel("Frequency (n=%s)" % len(h8shifts), size=17)
                plt.xlabel("H8 (ppm)", size=17)
                pdf.savefig()      

              #Plot C1' shifts
          if len(c1pshifts) != 0:
              if 2*iqr(c1pshifts[:]/len(c1pshifts)**(1/3.0)) != 0.0:
                plt.clf()
                h = 2*iqr(c1pshifts[:])/len(c1pshifts)**(1/3.0)
                nbins = np.ceil((max(c1pshifts[:]) - min(c1pshifts[:]))/h)
                plt.hist(c1pshifts[:], color='#0230a4', alpha=0.7, bins=nbins)
                plt.xlabel(("C1' (ppm)"), size=17)
                plt.ylabel("Frequency (n=%s)" % len(c1pshifts), size=17)
                pdf.savefig()

              else:
                plt.clf()
                plt.bar(c1pshifts[0], len(c1pshifts), align="center", color="#0230a4")
                plt.ylabel("Frequency (n=%s)" % len(c1pshifts), size=17)
                plt.xlabel("C1' (ppm)", size=17)
                pdf.savefig()     

              #Plot C6 shifts
          if len(c6shifts) != 0:
              if 2*iqr(c6shifts[:]/len(c6shifts)**(1/3.0)) != 0.0:
                plt.clf()
                h = 2*iqr(c6shifts[:])/len(c6shifts)**(1/3.0)
                nbins = np.ceil((max(c6shifts[:]) - min(c6shifts[:]))/h)
                plt.hist(c6shifts[:], color='#0230a4', alpha=0.7, bins=nbins)
                plt.xlabel(("C6 (ppm)"), size=17)
                plt.ylabel("Frequency (n=%s)" % len(c6shifts), size=17)
                pdf.savefig()

              else:
                plt.clf()
                plt.bar(c6shifts[0], len(c6shifts), align="center", color="#0230a4")
                plt.ylabel("Frequency (n=%s)" % len(c6shifts), size=17)
                plt.xlabel("C6 (ppm)", size=17)
                pdf.savefig()                


          if len(c8shifts) != 0:
                #Plot of C8 shifts
              if 2*iqr(c8shifts[:]/len(c8shifts)**(1/3.0)) != 0.00:
                  plt.clf()
                  h = 2*iqr(c8shifts[:])/len(c8shifts)**(1/3.0)
                  nbins = np.ceil((max(c8shifts[:]) - min(c8shifts[:]))/h)
                  plt.hist(c8shifts[:], color='#0230a4', alpha=0.7, bins=nbins)
                  plt.xlabel(("C8 (ppm)"), size=17)
                  plt.ylabel("Frequency (n=%s)" % len(c8shifts), size=17)
                  pdf.savefig()

              else:
                plt.clf()
                plt.bar(c8shifts[0], len(c8shifts), align="center", color="#0230a4")
                plt.ylabel("Frequency (n=%s)" % len(c8shifts), size=17)
                plt.xlabel("C8 (ppm)", size=17)
                pdf.savefig()      

      if xLen[0] > 6:
          #Plot of chi angles
        if len(chiAngles) != 0:
          if 2*iqr(chiAngles[:]/len(chiAngles)**(1/3.0)) != 0.00:
            plt.clf()
            h = 2*iqr(chiAngles[:])/len(chiAngles)**(1/3.0)
            nbins = np.ceil((max(chiAngles[:]) - min(chiAngles[:]))/h)
            plt.hist(chiAngles[:], color='#0416a4', alpha=0.7, bins=nbins)
            plt.xlabel(("Chi Angle"), size=17)
            plt.ylabel("Frequency (n=%s)" % len(chiAngles), size=17)
            pdf.savefig()

          #Plot of C1-C1 distance
        if len(c1c1) != 0:
          if 2*iqr(c1c1[:]/len(c1c1)**(1/3.0)) != 0.00:
            plt.clf()
            h = 2*iqr(c1c1[:])/len(c1c1)**(1/3.0)
            nbins = np.ceil((max(c1c1[:]) - min(c1c1[:]))/h)
            plt.hist(c1c1[:], color='#0260a4', alpha=0.7, bins=nbins)
            plt.xlabel(("C1-C1 Distance (Angstroms)"), size=17)
            plt.ylabel("Frequency (n=%s)" % len(c1c1), size=17)
            pdf.savefig()

        #Plot histograms of backbone torsions
        for k in range(6):
            plt.clf()
            h = 2*iqr(bbTorsions[:,k])/len(bbTorsions)**(1/3.0)
            nbins = np.ceil((max(bbTorsions[:,k]) - min(bbTorsions[:,k]))/h)
            plt.hist(bbTorsions[:,k], color='#23ef88', alpha=0.7, bins=nbins)
            plt.xlabel(("Backbone torision " + WN[k]), size=17)
            plt.ylabel("Frequency (n=%s)" % len(bbTorsions), size=17)
            pdf.savefig()

        #Plot histograms of local pair parameters
        for k in range(6):
            plt.clf()
            h = 2*iqr(pairPars[:,k])/len(pairPars)**(1/3.0)
            nbins = np.ceil((max(pairPars[:,k]) - min(pairPars[:,k]))/h)
            plt.hist(pairPars[:,k], color='#ffbd2e', alpha=0.7, bins=nbins)
            plt.xlabel(VN[k], size=17)
            plt.ylabel("Frequency (n=%s)" % len(pairPars), size=17)
            pdf.savefig()

        #Plot histograms of sugar dihedrals
        for k in range(5):
            plt.clf()
            h = 2*iqr(sugar[:,k])/len(sugar)**(1/3.0)
            nbins = np.ceil((max(sugar[:,k]) - min(sugar[:,k]))/h)
            plt.hist(sugar[:,k], color='#800000', alpha=0.7, bins=nbins)
            plt.xlabel(("Sugar dihedral " + XN[k]), size=17)
            plt.ylabel("Frequency (n=%s)" % len(sugar), size=17)
            pdf.savefig()

        ##################PCA of local pair geometry##################
        XS = (pairPars - pairPars.mean(0))/pairPars.std(0)
        # R = np.corrcoef(bbTorsions.T)
        # u,s,vt = np.linalg.svd(R)
        # Q = np.dot(XS, u)

        shape = XS.shape
        if shape[0] > shape[1]:
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Local Pair Geometry Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of local pair geometry##################

        ##################PCA of Backbone Torsions##################
        XS = (bbTorsions - bbTorsions.mean(0))/bbTorsions.std(0)
        # R = np.corrcoef(bbTorsions.T)
        # u,s,vt = np.linalg.svd(R)
        # Q = np.dot(XS, u)

        shape = XS.shape
        if shape[0] > shape[1]:
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Intrabase Backbone Torsion Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of Backbone Torsions##################

        ##################PCA of Backbone Torsions2##################
        XS = (bbTorsions2 - bbTorsions2.mean(0))/bbTorsions2.std(0)
        # R = np.corrcoef(bbTorsions.T)
        # u,s,vt = np.linalg.svd(R)
        # Q = np.dot(XS, u)
        shape = XS.shape
        if shape[0] > shape[1]:
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Interbase Backbone Torsion Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of Backbone Torsions2##################
        
        ##################PCA of sugar dihedrals##################
        XS = (sugar - sugar.mean(0))/sugar.std(0)
        shape = XS.shape
        if shape[0] > shape[1]:    

          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Intrabase Sugar Dihedral Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of sugar dihedrals##################

        ##################PCA of sugar dihedrals##################
        XS = (sugar2 - sugar2.mean(0))/sugar2.std(0)
        shape = XS.shape
        if shape[0] > shape[1]:
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Interbase Sugar Dihedral Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of sugar dihedrals##################

        ##################PCA of sugar + BB##################
        XS = (sugBB - sugBB.mean(0))/sugBB.std(0)
        shape = XS.shape
        if shape[0] > shape[1]:    
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Intrabase Backbone Torsions - Sugar Dihedral Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of sugar + BB##################

        ##################PCA of sugar + BB##################
        XS = (sugBB2 - sugBB2.mean(0))/sugBB2.std(0)
        shape = XS.shape
        if shape[0] > shape[1]:    
          plt.clf()
          results = PCA(XS) 

          #this will return an array of variance percentages for each component
          results.fracs

          #this will return a 2d array of the data projected into PCA space
          results.Y 

          x = []
          y = []
          z = []
          for item in results.Y:
            x.append(item[0])
            y.append(item[1])
            z.append(item[2])

          plt.close('all') # close all latent plotting windows
          fig1 = plt.figure() # Make a plotting figure
          ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
          pltData = [x,y,z] 
          ax.scatter(pltData[0], pltData[1], pltData[2], 'bo') # make a scatter plot of blue dots from the data
           
          # make simple, bare axis lines through space:
          xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
          ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
          yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
          ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
          zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
          ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
           
          # label the axes 
          str1 = str(round(results.fracs[0] * 100, 2))
          str2 = str(round(results.fracs[1] * 100, 2))
          str3 = str(round(results.fracs[2] * 100, 2))
          ax.set_xlabel("PC 1 (" + str1 + " %)") 
          ax.set_ylabel("PC 2 (" + str2 + " %)") 
          ax.set_zlabel("PC 3 (" + str3 + " %)") 
          ax.set_title('''Interbase Backbone Torsions - Sugar Dihedral Correlations
      Principal Component Analysis (n=%s)''' % (len(XS)))
          pdf.savefig()
        ##################PCA of sugar + BB##################
      pdf.close()