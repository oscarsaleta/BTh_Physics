reset
set terminal postscript eps color enhanced solid
set output "plot-375.eps"

set view map
set size ratio -1
set key textcolor rgb "white"
unset key
set xrange [-2.24:2.24]
set yrange [-2.24:2.24]
splot "375-wf.dat" u 1:2:3 w pm3d,\
      "375-3.dat" u 1:2:(0.0) w l title "(0.0,0.03)",\
      "375-4.dat" u 1:2:(0.0) w l title "(0.0,0.04)",\
      "375-5.dat" u 1:2:(0.0) w l title "(0.0,0.05)",\
      "375-6.dat" u 1:2:(0.0) w l title "(0.0,0.06)",\
      "375-7.dat" u 1:2:(0.0) w l title "(0.0,0.07)",\
      "375-8.dat" u 1:2:(0.0) w l title "(0.0,0.08)",\
      "375-9.dat" u 1:2:(0.0) w l title "(0.0,0.09)",\
      "375-10.dat" u 1:2:(0.0) w l title "(0.0,0.10)",\
      "375-11.dat" u 1:2:(0.0) w l title "(0.0,0.11)",\
      "375-12.dat" u 1:2:(0.0) w l title "(0.0,0.12)",\
      "375-13.dat" u 1:2:(0.0) w l title "(0.0,0.13)",\
      "375-14.dat" u 1:2:(0.0) w l title "(0.0,0.14)",\
      "375-15.dat" u 1:2:(0.0) w l title "(0.0,0.15)",\
      "375-16.dat" u 1:2:(0.0) w l title "(0.0,0.16)",\
      "375-17.dat" u 1:2:(0.0) w l title "(0.0,0.17)",\
      "375-18.dat" u 1:2:(0.0) w l title "(0.0,0.18)",\
      "375-19.dat" u 1:2:(0.0) w l title "(0.0,0.19)"

system("epstopdf plot-375.eps")
