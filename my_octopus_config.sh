#########################################################################
# File Name: myconfig.sh
# Author:kuangxiong 
# mail:kuangxiong@lsec.cc.ac.cn 
# Created Time: 2018年01月15日 星期一 19时59分41秒
#########################################################################
#!/bin/bash

./configure --prefix=/share/home/kuangx/opt/octopus-4.1.2 --with-blas=/share/home/kuangx/opt/lapack/lib/libblas.a --with-lapack=/share/home/kuangx/opt/lapack/lib/liblapack.a CC=mpicc FC=mpif90 CPP=cpp FCCPP="/lib/cpp -ansi -C -ffreestanding"  --with-gsl-prefix=/share/home/kuangx/octopus/gsl-1.9/ --with-scalapack=/share/home/kuangx/opt/scalapack-2.0.2/lib/libscalapack.a --with-libxc-prefix=/share/home/kuangx/opt/libxc-2.1.2 --with-fft-lib=/share/home/kuangx/opt/fftw/lib/libfftw3.a
