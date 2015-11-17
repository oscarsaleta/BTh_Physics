# $ gnuplot -e "name=filename_no_extension; k=index" plotWF.plt
set encoding utf8
set terminal pngcairo

set view map
set pm3d
unset surface

set output name.".png";
splot name.'.dat' using 1:2:3 index k
