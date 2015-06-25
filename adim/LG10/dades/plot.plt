reset
set encoding utf8
set terminal pdfcairo enhanced
set output 'L.pdf'

set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "#FFFFFF" behind

files=system('ls *-L.dat')

set grid
set key outside
set ytics 0.1
set yrange [-1.1:0.1]

#plot for [file in files] file w l title file
plot '125-L.dat' w l title '1/8',\
     '250-L.dat' w l title '1/4',\
     '375-L.dat' w l title '3/8',\
     '500-L.dat' w l title '1/2'
