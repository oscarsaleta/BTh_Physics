reset
set encoding utf8

set terminal pngcairo
set output "test.png"

set multiplot layout 1,2 title "Rotació d'àtoms freds"

set grid
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
unset key
set size ratio -1
#
set bmargin 5
set title "Clockwise: |10>+i|01>"
plot "f1.dat" w l lw 2, "<echo '0 1'" pt 7 ps 2
#
set title "Counterclockwise: |10>-i|01>"
plot "f-1.dat" w l lw 2, "<echo '0 1'" pt 7 ps 2
#
unset multiplot
