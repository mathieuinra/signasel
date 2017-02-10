#!/bin/bash

# Error control.
set -o nounset
set -o errexit


# Detecting libgsl
gsl=$( whereis libgsl.a | cut -f 2 -d ' ' )
if [ -z $gsl ]
then 
    exit 1
fi


# Creating the package
R -e "library(Rcpp);Rcpp.package.skeleton(\"signasel\", example_code = F)" &> log
echo -e "PKG_CXXFLAGS=-I/usr/include/ --std=c++1z -static\nPKG_LIBS=$gsl" > Makevars
cp max_f.hpp signasel.cpp Makevars signasel/src/

# Compiling 
R -e "library(Rcpp);compileAttributes(\"signasel\")" &>> log
R -e "library(devtools);install(\"signasel\")" &>> log



