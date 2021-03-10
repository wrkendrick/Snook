Snook
=====

LeastSquaresFit.h : the VectorPostProcessor that I hijacked to run HTPIPE within MOOSE. It's corresponding header is LeastSquaresFit.h.

htpipe.cpp and htpipe_v2.cpp : c++ files that run HTPIPE, with the former having a corresponding header file "include_htpipe.h".

input.i and input_v2.i : MOOSE input files, only use v2.

snook-opt : MOOSE executable.

heatflux_tester/: suite of test cases where I learn via toy problem.
  /input_v6.i : testing the HTPIPE c++ integration into the toy problem.
