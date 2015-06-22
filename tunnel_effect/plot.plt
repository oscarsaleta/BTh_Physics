reset
set encoding utf8

set terminal postscript eps color enhanced
set output 'tunnel_effect.eps'

set view map
set pm3d
unset surface
unset key
set cbrange [0:0.35]
set xrange [-8:8]
set yrange [-8:8]
#set size ratio -1

set object 1 circle at first -1.5,0 size first 1.5 front lc rgb "green"
set object 2 circle at first -1.5,0 size first 1.5 front lc rgb "green"

set multiplot layout 2,2 title "Tunneling effect in a double harmonic trap"

#0,0
set label 't=0.0' at 5.3,6 front tc ls 1
splot 'tunnel.dat' index 0 u 1:2:3
unset label

#0,1
set label 't=4.5' at 5.3,6 front tc ls 1
splot 'tunnel.dat' index 45 u 1:2:3
unset label

#1,0
set label 't=9.0' at 5.3,6 front tc ls 1
splot 'tunnel.dat' index 90 u 1:2:3
unset label

#0,1
set label 't=13.5' at 5.0,6 front tc ls 1
splot 'tunnel.dat' index 135 u 1:2:3
unset label

system('epstopdf tunnel_effect.eps')
