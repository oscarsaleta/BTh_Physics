# $ gnuplot -e "name=filename_no_extension" plotTR.plt
reset
set encoding utf8
set terminal pngcairo

set output name.".png"
plot name.'.dat' using 1:2 w l

set output name."-radi.png"
plot name.'.dat' using 3:4 w l
