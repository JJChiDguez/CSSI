#!/bin/sh

magma ./magma/set.m
cat ./magma/Arith.01 Arith.new ./magma/Arith.02 > Arith.S
mv ./setup_*.h ./inc
mv Arith.S ./lib
rm Arith.new
echo "We have updated Arith.S, setup_FF.h, and setup_EC.h!"
