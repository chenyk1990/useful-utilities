#!/bin/bash

module list

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
#sed -i "s:(cd ../create_header_file; make):#(cd ../create_header_file; make):g" src/specfem3D/Makefile
#sed -i "s:(cd ../create_header_file; make):#(cd ../create_header_file; make):g" src/auxiliaries/Makefile

##
## compiles xcreate_header_file
## (for static compilation)
#module switch xtpe-interlagos xtpe-istanbul
module switch craype-interlagos craype-istanbul
#module switch craype-interlagos craype-target-native
echo
module list
echo

make clean
make xcreate_header_file
rm -rf bin.header_file
cp -rp bin bin.header_file

module switch craype-istanbul craype-interlagos
echo
module list
echo

make clean

rm -rf OUTPUT_FILES/*
rm -rf bin.forward

./change_simulation_type.pl -f

mkdir -p bin
cp bin.header_file/xcreate_header_file bin/
make -j 8 xmeshfem3D
make -j 8 xspecfem3D 
make -j 8 xcombine_vol_data_vtk_adios

#if [ ! -e bin.forward/xspecfem3D ]; then exit; fi

