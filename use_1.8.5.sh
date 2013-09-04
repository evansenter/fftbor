#!/bin/bash
FFTBOR2D_VERSION="1.8.5"
git checkout fftbor2d
ruby -r fileutils -e "puts Dir['libRNA*'].map { |i| %x|readlink #{i}|.include?('$FFTBOR2D_VERSION') ? (FileUtils.mv(i, 'libRNA.a') unless i == 'libRNA.a') : FileUtils.mv(i, %x|readlink #{i}|.chomp.gsub(/\//, '_')) }"
make clean
make
make install