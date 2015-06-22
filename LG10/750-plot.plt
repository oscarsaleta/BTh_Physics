reset
set terminal postscript color enhanced
set output "plot-750.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
splot "750-wf.dat" u 1:2:3 w pm3d,\
      "750-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "750-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "750-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "750-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "750-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "750-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "750-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "750-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "750-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "750-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "750-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "750-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "750-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "750-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "750-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "750-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "750-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-750.eps")
