#!/bin/bash
#
# Issue multiple queries for the MultiDark Database
#
# This script splits the following query:
#
#     select * from MDR1..BDMV where snapnum = 85
#
# which would take about 4 min, into smaller queries which
# are executed in a loop.
# Also see http://www.multidark.org/MultiDark/Help?page=timeout
# for more explanations
# 
# NOTE: You need to add your own 
# 	username and password 
# for the MultiDark Database below!
#
# The MultiDark Database Team, latest update: Feb. 2012

username="robertrose"	# <--- your database user name
userpwd="adv32323"	# <--- your password


snapnum=85
plus="%2B"
times="%2A"

for (( i=0; i<20; i++ ))
do
   # fofId ranges from ($snapnum*10 + level) * 1e8 + rank in file  
   start=$i
   end=$(($i+1))

   sql="select * from MDR1..BDMV where bdmId between $snapnum $times 1e8 $plus $start $times 1e6 and $snapnum $times 1e8 $plus $end $times 1e6 - 1 and snapnum = $snapnum"

   wget --http-user=$username --http-passwd=$userpwd --cookies=on --keep-session-cookies --save-cookies=cookie.txt --load-cookies=cookie.txt -O result_$i.csv "http://wget.multidark.org/MyDB?action=doQuery&SQL=$sql" 


done
# -- takes 18 s for first bin, in total 14 bins

exit
