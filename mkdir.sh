#!/bin/sh

mkdir block
mkdir result
mkdir trunc
cat>Parameter<<EOF
nmax= 6
D= 21
LatticeSize= 20
gr= 1.25
gcr= 0
Jr= 0.005
Jcr= 0
EOF
