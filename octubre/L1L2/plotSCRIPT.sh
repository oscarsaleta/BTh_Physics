#!/bin/bash
for i in $(seq 0 $2); do
    gnuplot -e "name='$1';k=$i" plotWF.plt
done
