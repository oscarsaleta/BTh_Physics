reset
set encoding utf8
set terminal pdfcairo enhanced
set output 'L20.pdf'

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind
set grid
set key outside
set yrange [1.9:2.1]
set xlabel "Initial position y_0"
set ylabel "Angular momentum (ħ)"

plot '500-L.dat' w l title 'LG_{2,0}'
