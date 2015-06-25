reset
set terminal postscript eps color enhanced solid
set output "plot-500.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-3.1:3.1]
set yrange [-3.1:3.1]
#set palette grey negative
#set palette defined (0 "white", 0.15 "orange", 0.3 "red")
splot "500-wf.dat" u 1:2:3 w pm3d,\
      "500-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "500-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "500-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "500-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "500-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "500-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "500-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "500-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "500-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "500-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "500-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)",\
      "500-20.dat" u 1:2:(0.0) w l title "(0.0,0.20)",\
      "500-21.dat" u 1:2:(0.0) w l title "(0.0,0.21)",\
      "500-22.dat" u 1:2:(0.0) w l title "(0.0,0.22)",\
      "500-23.dat" u 1:2:(0.0) w l title "(0.0,0.23)",\
      "500-24.dat" u 1:2:(0.0) w l title "(0.0,0.24)",\
      "500-25.dat" u 1:2:(0.0) w l title "(0.0,0.25)",\
      "500-26.dat" u 1:2:(0.0) w l title "(0.0,0.26)"

system("epstopdf plot-500.eps")

