#!/bin/sh

samtools view -F 260 $1 | cut -f 3 | sort | uniq -c > uniq
awk '{printf("%s\t%s\n", $2, $1)}'  uniq > counts.txt
awk -F' ' "{ for(i=1; i<${2};i++) print ${1}" counts.txt | awk '{print ${1} }' | xargs -n 1 -I '{}' grep -A 1 {} $2 > $3
