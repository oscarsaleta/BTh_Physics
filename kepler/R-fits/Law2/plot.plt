reset
set encoding utf8
set terminal postscript eps color enhanced
set output 'kepler.eps'


set grid lc rgb "#808080" lt 0 lw 1

set xrange [0:50]
set xlabel 't'

set multiplot layout 1,2 margins 0.08,0.98,0.1,0.9 spacing 0.09,0.15

#set yrange [:0.517]
set title 'Area swept per unit time'
plot 'tr10-area.dat' u 3:($4<1 ? $4 : 1/0) w l title 'r_1',\
     'tr20-area.dat' u 3:($4<1 ? $4 : 1/0) w l title 'r_2'
#'tr050-area.dat' u 3:($4<1 ? $4 : 1/0) w l title 'r_1',\

#set yrange[:0.001]
unset key
set title 'Energy'
plot 'tr10-area.dat' u 3:($5<1 ? $5 : 1/0) w l title 'r_1',\
     'tr20-area.dat' u 3:($5<1 ? $5 : 1/0) w l title 'r_2'
#'tr050-area.dat' u 3:($5<0 ? $5 : 1/0) w l title 'r_1',\


unset multiplot

system('epstopdf kepler.eps && rm *.eps')
