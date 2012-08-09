#! /usr/bin/env python

# graphFFTborOutput.py
# P. Clote

# GNUPLOT can be used in batch mode, where you type
# "gnuplot infile", where infile contains gnuplot commands.
# Using os.system(cmd), where
#               cmd = "gnuplot " + sys.argv[1]
# and sys.argv[1] contains the name (say "infile") containing
# the batch commands, you can invoke gnuplot from within a
# python script. A better way is to use the module tempfile
# to create a temporary file name and write the batch commands
# to that file, rather than by hand preparing the file "infile".
# This is done below.

import sys,os,tempfile
from os.path import splitext

def aux0(infilename,xrange,yrange):
  infilename0 = splitext(infilename)[0]
  filename    = tempfile.mktemp()
  file        = open(filename,"w")
  text        = """set noautoscale
  set xrange [0:%d]
  set yrange [0:%f]
  set terminal postscript eps
  set output '%s.eps'
  plot '%s' with linespoints """ % (xrange,yrange,infilename0,infilename)
  file.write(text)
  file.flush()
  cmd = "gnuplot %s" % filename
  os.system(cmd)
  cmd = "convert %s.eps %s.gif" % (infilename0,infilename0)
  os.system(cmd)
#  cmd = "rm -f %s.eps" % infilename
#  os.system(cmd)
#  cmd = "convert %s.eps %s.pdf" % (infilename,infilename)
#  os.system(cmd)
  return


def main(filenames,xrange,yrange):
  text  = "convert -delay %d -loop 0 " % delay
  for infilename in filenames:
    aux0(infilename,xrange,yrange)
    text += "%s.gif " % splitext(infilename)[0]
  text += "fftborAnimated.gif"
  print text
  os.system(text)
  cmd   = "chmod og+r fftborAnimated.gif"
  os.system(cmd)
    

if __name__ == '__main__':
  if len(sys.argv) < 8:
    print "Usage: %s -x xrange -y yrange -d delay graphicsFile(s)" % sys.argv[0]
    print "\t(1) xrange is (integer) max RNA sequence length from FFTbor output"
    print "\t(2) yrange is (float <= 1) max probability in FFTbor output"
    print "\t(3) delay (integer): num of 1/100 sec for ImageMagick animation"
    print "\t(4) graphics file(s), each histogram, or single column of numbers"
    sys.exit(1)
  assert( sys.argv[1]=='-x' )
  xrange     = int(sys.argv[2])
  assert( sys.argv[3]=='-y' )
  yrange     = float(sys.argv[4])
  assert( sys.argv[5]=='-d' )
  delay      = int(sys.argv[6])
  filenames  = sys.argv[7:]
  main(filenames,xrange,yrange)


