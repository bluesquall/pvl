#!/bin/bash

# build-script
#   a short BASH script to build original PVL from fortran source
#   TODO use GNU development tools/autotools in the future

rm ./*.o

echo "old modules cleaned, building modules"
gfortran -c splutil.f90
gfortran -c Wrench.f90
gfortran -c Simeqn.f90
gfortran -c Pvlmod.f90
gfortran -c Forces.f90
gfortran -c Pvl.f90

echo "modules built--linking pvl executable"
gfortran *.o -O2 -Wall -o pvl

echo "pvl built"
echo "try running: "
echo "echo 'example.inp' | ./pvl"
