reset
set encoding utf8
set terminal pdfcairo enhanced
set output 'L.pdf'

files=system('ls *-L.dat')

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind
set grid
set key outside
set ytics 0.1
set yrange [1.8:2.3]

#plot for [file in files] file w l title file
plot '500-L.dat' w l title 'LG_{2,0}'
