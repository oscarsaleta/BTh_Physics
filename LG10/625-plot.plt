reset
set terminal postscript color enhanced solid
set output "plot-625.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.1:2.1]
set yrange [-2.1:2.1]
splot "625-wf.dat" u 1:2:3 w pm3d,\
      "625-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "625-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "625-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "625-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "625-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "625-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "625-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "625-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "625-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "625-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "625-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "625-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "625-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "625-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "625-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "625-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "625-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-625.eps")

