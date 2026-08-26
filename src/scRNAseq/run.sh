#!/bin/bash
mapfile -t arr < a.txt
for i in ${arr[@]}
do
    echo $i
    conda remove $i --yes
done