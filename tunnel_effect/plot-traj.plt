reset
set encoding utf8

set terminal postscript eps color enhanced solid
set output 'tunnel-tr.eps'
set key outside
set xlabel 'Time'
set xrange [0:20]
set ylabel 'Position'

set grid

plot 'tunnel-tr-10.dat' u 3:1 title 'x_{(-2.5,0)}' w l lw 2, '' u 3:2 title 'y_{(-2.5,0)}' w l lw 2,\
     'tunnel-tr00.dat' u 3:1 title 'x_{(-1.5,0)}' w l lw 2, '' u 3:2 title 'y_{(-1.5,0)}' w l lw 2,\
     'tunnel-tr10.dat' u 3:1 title 'x_{(-0.5,0)}' w l lw 2, '' u 3:2 title 'y_{(-0.5,0)}' w l lw 2

system('epstopdf tunnel-tr.eps')
