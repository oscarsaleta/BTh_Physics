reset
set terminal postscript color enhanced
set output "plot-500.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
splot "500-wf.dat" u 1:2:3 w pm3d,\
      "500-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "500-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "500-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "500-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "500-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "500-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
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
      "500-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-500.eps")

