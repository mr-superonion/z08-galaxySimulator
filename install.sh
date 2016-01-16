#!/bin/sh

f2py -m GalSP -c ./lib/*.f95 ./lib/*.f90 ./lib/*.f 

cp GalSP.so /home/xiangchong/Mylib
