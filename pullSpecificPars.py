import urllib, urllib2
import unicodedata
import xml.etree.ElementTree as ET
from xml.dom.minidom import parseString

def parseXmlQuery(inArray):
  data = []
  for answer_dict in inArray:
    data.append(answer_dict)  
  return (data)
pdbArray = ["1IK5", "1S72", "131D", "1D48", "1ENN", "1O56", "292D", "2DCG", "3U89"]
for xx in pdbArray:
  pdbId = xx
  req = urllib2.Request("http://www.rcsb.org/pdb/rest/customReport?pdbids=%s&customReportColumns=resolution,crystallizationTempK,phValue,pdbxDetails,pubmedId" % pdbId)
  post = urllib2.urlopen(req)
  root = ET.fromstring(post.read())
dicta = {}
#pdb id
print root[0][0].text
#resolution
print root[0][1].text
#temperature
print root[0][2].text
#ph
print root[0][3].text
#conditions
str4 = root[0][4].text.replace(",", " |")
print str4
dicta[root[0][0].text] = root[0][4].text
print ""
print dicta
