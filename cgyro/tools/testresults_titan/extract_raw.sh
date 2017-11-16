#!/bin/bash
files=`/bin/ls */*/out.cgyro.timing | grep -v '_s'`
for f in $files; do 
grep -A 4 str_com $f > a.txt; lines=`cat a.txt | wc -l`; 
if [ $lines -gt 2 ]; then 
  el=`echo $f | awk -F "/" '{print $2}'`
  if [ $lines -gt 3 ]; then 
     echo -n "$el " ; tail -3 a.txt | head -1
     echo -n "$el " ; tail -2 a.txt | head -1
     echo -n "$el " ; tail -1 a.txt
   elif [ $lines -gt 3 ]; then
     echo -n "$el " ; tail -2 a.txt | head -1
     echo -n "$el " ; tail -1 a.txt
   else
     echo -n "$el " ;   tail -1 a.txt; 
  fi; 
fi; 
done
