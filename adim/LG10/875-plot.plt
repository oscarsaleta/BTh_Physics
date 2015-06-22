reset
set terminal postscript color enhanced solid
set output "plot-875.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
splot "875-wf.dat" u 1:2:3 w pm3d,\
      "875-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "875-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "875-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "875-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "875-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "875-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "875-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "875-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "875-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "875-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "875-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "875-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "875-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "875-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "875-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "875-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "875-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-875.eps")

