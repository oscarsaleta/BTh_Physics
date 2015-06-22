reset
set terminal postscript eps color enhanced solid
set output "plot-250.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
#set palette grey negative
#set palette defined (0 "white", 0.15 "orange", 0.3 "red")
splot "250-wf.dat" u 1:2:3 w pm3d,\
      "250-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "250-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "250-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "250-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "250-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "250-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "250-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "250-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "250-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "250-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "250-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "250-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "250-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "250-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "250-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "250-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "250-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-250.eps")

