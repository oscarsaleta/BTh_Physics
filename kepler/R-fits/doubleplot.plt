reset
set encoding utf8
set terminal postscript eps color solid enhanced
set output 'kepler-wftr.eps'

set view map
set pm3d
unset colorbox

set title "Trajectories for V=-1/r"
set xlabel "x"
set ylabel "y"

set xrange [-2.24:2.24]
set yrange [-2.24:2.24]

set size ratio -1

set key outside spacing 1.3

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind

spl '../kepler-wf.dat' u 1:2:3 index 0 w pm3d notitle,\
        'tr10.dat' u 1:2:(0.0) every 10::::625 w p pt 6 title 'r_1=(1.0,0.0)',\
        'tr10-fit.dat' u 2:3:(0.0) w l lw 2 title 'Fitted r_1',\
        'tr20.dat' u 1:2:(0.0) every 20::::2500 w p pt 6 title 'r_2=(2.0,0.0)',\
        'tr20-fit.dat' u 2:3:(0.0) w l lw 2 lc rgb "red" title 'Fitted r_2'
system('epstopdf kepler-wftr.eps && rm *.eps')
#'tr050.dat' u 1:2:(0.0) every 3::::154 w p pt 6 title 'r_1=(0.5,0.0)',\
#'tr050-fit.dat' u 2:3:(0.0) w l lw 2 lc rgb "blue" title 'Fitted r_1',\
