#########################################################################
# File Name: copy.sh
# Author:kuangxiong 
# mail:kuangxiong@lsec.cc.ac.cn 
# Created Time: 2018年01月15日 星期一 15时00分13秒
#########################################################################
#!/bin/bash
scp -r kuangx@10.4.3.15:mkl_KPP ./
git add .
git commit -m "KPP-POD201801124"
git push origin master
git log
