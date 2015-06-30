reset
set encoding utf8
set terminal postscript eps color enhanced
set output 'kepler-area-nosim.eps'


set grid lc rgb "#808080" lt 0 lw 1

set xrange [0:50]
set xlabel 't'

set multiplot layout 1,2 margins 0.08,0.98,0.1,0.9 spacing 0.09,0.15

#set yrange [:0.517]
set title 'Area swept per time unit'
plot 'ab-tr10-area.dat' u 3:($4<10 ? $4 : 1/0) w l title 'r_1',\
     'ab-tr11-area.dat' u 3:($4<10 ? $4 : 1/0) w l title 'r_2'

#set yrange[:0.001]
unset key
set title 'Energy'
plot 'ab-tr10-area.dat' u 3:($5<10 ? $5 : 1/0) w l title 'r_1',\
     'ab-tr11-area.dat' u 3:($5<10 ? $5 : 1/0) w l title 'r_2'


unset multiplot

system('epstopdf kepler-area-nosim.eps && rm *.eps')
