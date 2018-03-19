#!/bin/bash
files=`/bin/ls */*/out.cgyro.timing`
for f in $files; do 
grep -A 12 str_com $f > a.txt; lines=`cat a.txt | wc -l`; 
if [ $lines -gt 11 ]; then 
  el=`echo $f | awk -F "/" '{print $2}'`
  echo -n "$el " ;   tail -5 a.txt | head -1
  echo -n "$el " ;   tail -4 a.txt | head -1
  echo -n "$el " ;   tail -3 a.txt | head -1
  echo -n "$el " ;   tail -2 a.txt | head -1
  echo -n "$el " ;   tail -1 a.txt 
fi; 
done
