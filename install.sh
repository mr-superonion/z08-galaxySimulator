#!/bin/sh

# compile the library
f2py -m GalSP -c ./lib/*.f95 ./lib/*.f90 ./lib/*.f 

# put the library in the specific PYTHONPATH
cp GalSP.so /home/xiangchong/Mylib
