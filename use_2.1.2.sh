#!/bin/bash
FFTBOR2D_VERSION="2.1.2"
git checkout fftbor2d_vienna_2.x
ruby -r fileutils -e "puts Dir['libRNA*'].map { |i| %x|readlink #{i}|.include?('$FFTBOR2D_VERSION') ? (FileUtils.mv(i, 'libRNA.a') unless i == 'libRNA.a') : FileUtils.mv(i, %x|readlink #{i}|.chomp.gsub(/\//, '_')) }"
make clean
make
make install