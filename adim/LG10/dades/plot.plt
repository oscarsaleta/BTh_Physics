reset
set encoding utf8
set terminal pdfcairo enhanced
set output 'L.pdf'

files=system('ls *-L.dat')

set grid
set key box outside
set ytics 0.1
set yrange [-1.1:1.1]

#plot for [file in files] file w l title file
plot '125-L.dat' w l title '1/8',\
     '250-L.dat' w l title '1/4',\
     '375-L.dat' w l title '3/8',\
     '500-L.dat' w l title '1/2',\
     '625-L.dat' w l title '5/8',\
     '750-L.dat' w l title '3/4',\
     '875-L.dat' w l title '7/8'
