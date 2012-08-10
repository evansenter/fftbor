#! /usr/bin/env python

# createFFTborAnimationUsingJavascript.py
# P.Clote

import sys,os
from os.path import splitext,basename

MARKER = "+++"

def aux(rnaSeqStrFileName):
  file = open(rnaSeqStrFileName)
  line = file.readline()
  while line.strip()[0]==">":
    line = file.readline()
  rna    = line.strip().upper()
  secStr = file.readline().strip()
  file.close()
  print rna
  print secStr
  return rna,secStr

def rnaLenForGifFiles(gifFileNames):
  #determine min length of RNA in gifFileNames
  L = []
  for filename in gifFileNames:
    L.append( int(splitext(basename(filename))[0]) )
  return L


def main(rnaSeqStrFileName,gifFileNames):
  rna,secStr = aux(rnaSeqStrFileName)
  text = ""
  file = open("indexFFTborAnimationTemplate.html")
  #search for +++1
  line = file.readline()
  while line and line[:3]!=MARKER:
    text += line
    line = file.readline()
  text += "%s\n" % rna
  #search for +++2
  line = file.readline()
  while line and line[:3]!=MARKER:
    text += line
    line = file.readline()
  text += "%s\n" % secStr
  text += "<p/><hr/><p/>\n"
  text += "<h2>Individual gif files to peruse manually</h2>\n"
  text += "<a href='indexGifFiles.html'>indexGifFiles.html</a>\n"
  text += "<h2>Animation of gif files below</h2>\n"
  #now create indexGifFiles.html
  outfile = open('indexGifFiles.html','w')
  outtext = "<html>\n"
  outtext += "<body bgcolor='white'>\n"
  rnaLengths = rnaLenForGifFiles(gifFileNames)
  rnaLengths.sort(); gifFileNames.sort()
  outtext = "<ul>\n" 
  assert( len(rnaLengths) == len(gifFileNames) )
  N  = len(rnaLengths)
  for k in range(N):
    rnaLen      = rnaLengths[k]
    filename    = gifFileNames[k]
    outtext += "<li><b>%d</b>\n" % rnaLen
    outtext += "<img src='%s' alt='%s'></img></li>\n" % (filename,filename)
  text += "</ol>\n"
  outtext += "</body>\n"
  outtext += "</html>\n"
  outfile.write(outtext)
  cmd = "chmod og+r indexGifFiles.html"
  os.system(cmd)

  #search for +++3
  line = file.readline()
  while line and line[:3]!=MARKER:
    text += line
    line = file.readline()
  for filename in gifFileNames:
    text += "<li><img src='%s' alt='%s'></img></li>\n" % (filename,filename)
  text += "<!-- END of graphics files -->\n"
  #now read until end of file
  line = file.readline()
  while line:
    text += line
    line = file.readline()
  file.close()
  file = open("indexFFTborOutputAnimation.html",'w')
  file.write(text)
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 4:
    text = """Usage: %s rnaSecStrFilename gifFilename(s) 
    (1) first file contains rna seq and initial sec str
    (2) gifFileNames are possibly in a subdirectory such as GIF or OUTPUT
Use * to group all files of form XXX.gif where XXX are digits"""
    print text % sys.argv[0]
    sys.exit(1)
  rnaSeqStrFileName = sys.argv[1]
  gifFileNames      = sys.argv[2:]
  main(rnaSeqStrFileName,gifFileNames)


