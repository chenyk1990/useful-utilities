#!/bin/bash

#see: ~/.bashrc
source /etc/profile

module list

## cross-compiling:

## compilation for xcreate_header_file
## to run on login node
echo
echo "compiling xcreate_header_file..."
echo

# changes target rules
if [ -e src/create_header_file/rules.mk ]; then
  # new makefile - rules
  sed -i 's:@-rm -f $@.*:@-rm -f $@; cp bin.header_file/xcreate_header_file $E/:' src/create_header_file/rules.mk
else
  # old makefile
  sed -i "s:(cd ../create_header_file; make):#(cd ../create_header_file; make):g" src/specfem3D/Makefile
  sed -i "s:(cd ../create_header_file; make):#(cd ../create_header_file; make):g" src/auxiliaries/Makefile
fi

##
## compiles xcreate_header_file
## (for static compilation)
module switch craype-interlagos craype-target-native
echo
module list
echo

make clean
make xcreate_header_file
#if [ ! -e bin/xcreate_header_file ]; then exit; fi
rm -rf bin.header_file
cp -rp bin bin.header_file

module switch craype-target-native craype-interlagos
echo
module list
echo


# forward simulations
#
# only needed to produce synthetics for flexwin/measurements
#
echo
echo "compiling forward simulations..."
echo
make clean
sed -i "s/.*READ_ADJSRC_ASDF                =.*/READ_ADJSRC_ASDF                =.true./" DATA/Par_file

rm -rf OUTPUT_FILES/*
rm -rf bin.forward

mkdir -p bin
cp bin.header_file/xcreate_header_file bin/
make -j 8 xmeshfem3D
make -j 8 xspecfem3D

cp -rp bin bin.forward
cp OUTPUT_FILES/* bin.forward/
