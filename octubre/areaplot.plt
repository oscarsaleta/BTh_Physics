reset
set encoding utf8
set terminal postscript eps color enhanced
set output 'kepler-area.eps'


set grid lc rgb "#808080" lt 0 lw 1

set xrange [0:5]
set xlabel 't'

set multiplot layout 1,2 margins 0.08,0.98,0.1,0.9 spacing 0.09,0.15

set yrange [:3.3]
set title 'Area swept per time unit'
plot 'dades/r1-area.dat' u 3:($4<1 ? $4 : 1/0) w l title 'r_1'#,\
#'r2-area.dat' u 3:($4<1 ? $4 : 1/0) w l title 'r_2'
unset key
set yrange [-2.8:0]
set title 'Energy'
plot 'dades/r1-area.dat' u 3:($5<1 ? $5 : 1/0) w l title 'r_1'#,\
#     'r2-area.dat' u 3:($5<1 ? $5 : 1/0) w l title 'r_1'


unset multiplot

system('epstopdf kepler-area.eps && rm *.eps')
