reset
set terminal postscript color enhanced solid
set output "plot-125.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
splot "125-wf.dat" u 1:2:3 w pm3d,\
      "125-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "125-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "125-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "125-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "125-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "125-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "125-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "125-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "125-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "125-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "125-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "125-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "125-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "125-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "125-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "125-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "125-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-125.eps")
