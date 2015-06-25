#!/bin/bash
fase=$1
name1a=$(echo "$fase*1000" | bc)
name2a=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" ${name1a})
name4="dades/$name2a-L.dat"
$(rm $name4)
for i in $(LC_NUMERIC="en_US.UTF-8" seq 2.6 -0.1 0.9); do
    name1b=$(echo "$i*10" | bc)
    name2b=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f" ${name1b})
    name3="${name2a}-${name2b}.dat"
    if [ ! -f "$name3" ]; then
        $(./wavefunction 0.01 25 "$fase" 0 "$i" > "$name3")
    fi
    nlines=$(wc -l "$name3" | awk '{print $1;}')
    L="$(./main_L $name3 $nlines)"
    echo "$i $L" >> $name4
done
