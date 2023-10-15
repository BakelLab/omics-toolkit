#!/bin/sh

#################################
# CABOG (WGS) COMPILATION NOTES #
#################################

# Make a place to build the assembler. All further steps should start 
# from within this directory.
mkdir wgs
cd wgs

# Download and compile 'Figaro'. The cabog compilation expects to find 
# the figaro installation in ../figaro/Linux-amd64
# tested:   gcc 4.1.2  => works (compiled on Redhat ES 5.5)
#           gcc 4.4.0  => fails with declaration errors
#           icc 11.1   => works
wget http://surfnet.dl.sourceforge.net/project/amos/Figaro/1.05/Figaro-1.05.tar.gz
tar zxvf Figaro-1.05.tar.gz
mkdir -p figaro/Linux-amd64
mv Figaro-1.05/* figaro/Linux-amd64
rmdir Figaro-1.05
cd figaro/Linux-amd64
gmake
gmake install
cd ../../

# Download and compile 'kmer'. Once again, the cabog compilation expects
# to find the kmer installation in ../kmer/Linux-amd64
# tested:   gcc 4.1.2  => works (compiled on Redhat ES 5.5)
#           gcc 4.4.0  => fails with declaration errors  (add: #include <stdlib.h> to kmerdist.cc)
#           icc 11.1   => works (following instructions below)
svn co https://svn.code.sf.net/p/kmer/code/trunk kmer
cd kmer
rm -rf scripts seagen seatac sim4db sim4dbutils snapper tapper trie leaff developer-doc atac-driver ESTmapper libsim4
perl -pi -e 's/gcc/icc/' configure.sh
perl -pi -e 's/g\+\+/icpc/' configure.sh
perl -pi -e 's/-O3 /-O3 -xHost /g' configure.sh
perl -pi -e 's/-fexpensive-optimizations //g' configure.sh
perl -pi -e 's/-Wno-char-subscripts //g' configure.sh
sh configure.sh
gmake
gmake install
cd ..

# Download and compile wgs. This will only work if the previous
# compiles finished successfully
# tested:   gcc 4.1.2  => works (compiled on Redhat ES 5.5)
#           gcc 4.4.0  => not tested
#           icc 11.1   => works
svn checkout svn://svn.code.sf.net/p/wgs-assembler/svn/trunk/src src
cd src
perl -pi -e 's/gcc/icc/' c_make*
perl -pi -e 's/g\+\+/icpc/' c_make*
perl -pi -e 's/-O2 /-O2 -xHost /g' c_make*
perl -pi -e 's/-O4 /-O3 -xHost /g' c_make*
perl -pi -e 's/-Wno-char-subscripts//g' c_make*
perl -pi -e 's/-fexpensive-optimizations//g' c_make*
gmake
cd ..

# After the last compilation, the suite of cabog programs will be located in the 
# wgs/Linux-amd64 folder
