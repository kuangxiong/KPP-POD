#########################################################################
# File Name: copy.sh
# Author:kuangxiong 
# mail:kuangxiong@lsec.cc.ac.cn 
# Created Time: 2018年01月15日 星期一 15时00分13秒
#########################################################################
#!/bin/bash
scp -r lssc4:KPP ./
git add .
git commit -m "KPP-POD201801118"
git push origin master
git log
