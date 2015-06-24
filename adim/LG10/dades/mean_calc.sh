#!/bin/bash
$(rm Lmeans.dat)
for i in $(seq 1 4); do
    fase=$(echo "$i*0.125" | bc)
    name1a=$(echo "$fase*1000" | bc)
    name2a=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" ${name1a})
    name4="$name2a-L.dat"
    nlines=$(wc -l "$name4" | awk '{print $1;}')
    echo "$fase $(./main_Lm $name4 $nlines)" >> "Lmeans.dat"
done
