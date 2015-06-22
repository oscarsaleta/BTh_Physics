reset
set encoding utf8

set terminal postscript eps color enhanced
set output 'tunnel_effect.eps'

set view map
set pm3d
unset surface
set cbrange [0:0.35]
set xrange [-8:8]
set yrange [-8:8]
set size ratio -1

set object 1 circle at first -1.5,0 size first 1.5
set object 2 circle at first -1.5,0 size first 1.5

set multiplot layout 2,2 title "Tunneling effect in a double harmonic trap"

#0,0
set label 't=0' at 5.8,6 front tc ls 1
splot 'tunnel.dat' index 0 u 1:2:3
unset label

#0,1
set label 't=
