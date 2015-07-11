reset
set terminal postscript eps color enhanced
set output 'kepler-wftr-nosim.eps'

set view map
set key outside spacing 1.3
unset colorbox
set size ratio -1

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind


set xrange [-3.072:3.072]
set yrange [-3.072:3.072]

splot '3ab.dat' index 0 u 1:2:3 w pm3d notitle,\
    'ab-tr10.dat' u 1:2:(0.0) every 10::::1890 w p pt 6 lc rgb "green" title "r_3=(1,0)",\
    'ab-tr10-fit.dat' u 2:3:(0.0) w l lw 2 lc rgb "red" title "Fitted r_3",\
    'ab-tr11.dat' u 1:2:(0.0) every 10::::3900 w p pt 6 title "r_4=(1,1)",\
    'ab-tr11-fit.dat' u 2:3:(0.0) w l lw 2 lc rgb "orange-red" title "Fitted r_4"

system('epstopdf kepler-wftr-nosim.eps && rm *.eps')

