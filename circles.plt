# define x0,y0,name
set xrange [-5:5]
set yrange [-5:5]
set grid

ini = "<echo '".x0." ".y0."'"

plot name w l title "Trajectory",\
        "<echo '-1.9 0 0.5'" w circles lc rgb "blue" fs transparent solid 0.15 noborder notitle,\
        "<echo '1.9 0 0.5'" w circles lc rgb "blue" fs transparent solid 0.15 noborder notitle,\
        ini pt 7 title "Initial position"

