#! /usr/bin/env python

# runFFTborOnIncreasingRNAprefixesUsingMFEstr.py
# P.Clote

# Program takes an input file name, which contains
# a single RNA sequence all on one line, then runs
# FFTbor on all prefixes of the RNA, outputting files
# with the naming convention rnaXXX.txt, where XXX is
# the length of the RNA prefix. Then program runs FFTbor
# on rnaXXX.txt and outputs the result in rnaXXX.txt.output.
# This latter file contains (1) FASTA comment, (2) RNA sequence
# which FFTbor is run on (i.e. a prefix of the original sequence),
# (3) MFE structure of RNA sequence on which FFTbor is run, then
# the Boltzmann probability, followed by 2 summary lines at end.

# There are CONSTANTS that define the starting length of minimum length
# prefix, and the increment amount of each new prefix.

# Program trimFFTBORoutput.py then removes the top 3 lines and last
# 2 lines of this output for graphical processing.

import sys,os
from os.path import splitext, curdir


OUTPUTDIR = "OUTPUT" #directory for output files

NumDigits = 3  
  #output files of form XXX where string representing integer has len NumDigits
  #value 3 should suffice, but if RNA len is > 999 then set to 4
def checkPresenceOfEnergyParameterFile():
  #check if energy.par is in current directory
  directoryContents = os.listdir(curdir)
  if 'energy.par' not in directoryContents:
    print "File energy.par not in current directory!"
    sys.exit(1)
  return 


def aux(filename):
  #code assumes that filename contains a possible FASTA comment, then
  #a SINGLE RNA sequence, possibly broken into separate lines, and nothing
  #else. For instance input could be FASTA format with fixed line length of 60
  rna0 = ""
  file = open(filename)
  line = file.readline()
  #skip FASTA comment(s)
  while line.strip()[0]=='>':
    line = file.readline()
  #read RNA sequence into variable rna0
  while line:
    rna0 += line.strip()
    line  = file.readline()
  file.close()
  rna0 = rna0.upper()
  return rna0

def MFEstr(rna):
  #computes MFE structure using Vienna RNA package RNAfold
  cmd    = "echo %s | RNAfold" % rna
  file   = os.popen(cmd)
  line   = file.readline()
  assert( rna == line.strip() )
  secStr = file.readline().split()[0]
  return secStr

def reformatFFTborOutputFile(filename):
  #delete first 3 lines and last 2 lines to get only the Boltz probabilities
  n    = 0; L = []
  file = open(filename)
  line = file.readline()
  while line:
    n += 1
    L.append(line) #line terminates with \n
    line  = file.readline()
  file.close()
  L = L[3:-2]
  return L


def main(filename,startPos,windowIncrement):
  #check if OUTPUTDIR in current directory, and if not create this directory
  directoryContents = os.listdir(curdir)
  if OUTPUTDIR not in directoryContents:
    os.mkdir(OUTPUTDIR)
  #get input RNA sequence
  rna0   = aux(filename)
  length = startPos
  rna    = rna0[:length]
  #generate prefixes of the input RNA sequence and run FFTbor
  while length<=len(rna0):
    outfilename = "%d.txt" % length
    #prepend zeros to output file name so integer has NumDigits digits
    #this is done so that lexicographic order corresponds to integer order
    #when calling convert to produce ImageMagick animation
    lenDiff     = NumDigits - len(str(length))
    for k in range(lenDiff):
      outfilename = '0'+outfilename
    outfilename = OUTPUTDIR + os.sep + outfilename
    file        = open(outfilename,'w')
    txt         = "> len is %d\n" % length
    file.write(txt)
    file.write("%s\n" % rna)
    secStr      = MFEstr(rna)
    file.write("%s\n" % secStr)
    file.close() 
    fftoutfilename = splitext(outfilename)[0] + ".output0"
    cmd  = "FFTbor %s > %s" % (outfilename,fftoutfilename)
    os.system(cmd)
    L    = reformatFFTborOutputFile(fftoutfilename)
    file = open(splitext(outfilename)[0]+".output1",'w')
    for line in L:
      file.write(line)
    file.close() 
    length     += windowIncrement
    rna         = rna0[:length]


if __name__ == '__main__':
  checkPresenceOfEnergyParameterFile()
  if len(sys.argv) < 2:
    text = "Usage: %s filename [startPos [windowIncrement]]\n" % sys.argv[0]
    text+= "\nProgram computes prefixes of RNA sequence given in input file\n"
    text+= "and then runs FFTbor on each prefix, naming the output files\n"
    text+= "as rnaXXX.txt for prefix and rnaXXX.output for FFTbor output\n"
    text+= "Here XXX may range from 000 to 999.\n"
    text+= "\t(1) startPos is minimum prefix length (default = 10)\n"
    text+= "\t(2) startPos is sliding window advancement length (default 1)\n"
    print text
    sys.exit(1)
  filename = sys.argv[1]
  if len(sys.argv)>2:
    startPos = int(sys.argv[2])
  else:
    startPos = 10
  if len(sys.argv)>3:
    windowIncrement = int(sys.argv[3])
  else:
    windowIncrement = 1
  main(filename,startPos,windowIncrement)

