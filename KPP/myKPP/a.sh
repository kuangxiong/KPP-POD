#########################################################################
# File Name: a.sh
# Author:kuangxiong 
# mail: kuangxiong@lsec.cc.ac.cn
# Created Time: 2017年12月27日 星期三 09时42分01秒
#########################################################################
#!/bin/bash
make clean
make main
mpirun -np 4 ./main

