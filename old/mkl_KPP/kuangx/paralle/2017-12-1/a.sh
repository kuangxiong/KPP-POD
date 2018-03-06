#########################################################################
# File Name: make.sh
# Author:kuangxiong 
# mail:kuangxiong@lsec.cc.ac.cn 
# Created Time: 2017年12月18日 星期一 22时48分02秒
#########################################################################
#!/bin/bash
make clean
make main
mpirun -np 4 ./main
