#!/bin/sh

f2py -m GalSP -c ./lib/*.f95 ./lib/*.f > ./examples/GalSP.so

