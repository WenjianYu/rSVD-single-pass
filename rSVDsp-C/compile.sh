#!/bin/bash
icc -O3 -mkl -openmp rsvd.c rsvd_test.c matrix_funs_intel_mkl.c memory_usage.c -o rSVDsp_test
