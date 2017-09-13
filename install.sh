#!/bin/sh

# compile the library
f2py -m GalSP -c ./lib/*.f95 ./lib/*.f90 ./lib/*.f 

mv GalSP.so lib
