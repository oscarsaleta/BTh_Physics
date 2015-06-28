reset
set encoding utf8
set terminal postscript eps color enhanced
set output 'kepler.eps'

set title "Trajectories for a -1/r potential"
set xlabel "x"
set ylabel "y"
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]

set key outside spacing 1.3

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind

set grid lc rgb "#808080" lt 0 lw 1

plot 'tr10.dat' u 1:2 every 5::::625 w p pt 6 title 'r_0=(1.0,0.0)',\
     'tr10-fit.dat' u 2:3 w l lw 2 title 'Fitted r_0',\
     'tr050.dat' u 1:2 every 3::::154 w p pt 6 title 'r_1=(0.5,0.0)',\
     'tr050-fit.dat' u 2:3 w l lw 2 title 'Fitted r_1'

system('epstopdf "kepler.eps" && rm *.eps')
