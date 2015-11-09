# Llamar el programa como
# $ gnuplot -e "N=número; dt=número; x0=x_inicial; xf=x_final;" gnupscript.plt
#set terminal postscript color enhanced "Helvetica" 16
#set output 'test2.eps'
set encoding utf8
set terminal pngcairo
unset key
set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
#set nocbtics
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set xrange [x0:xf] noreverse nowriteback
set yrange [y0:yf] noreverse nowriteback
set cbrange [0:0.55] noreverse nowriteback
set palette rgbformulae 7, 5, 15


do for [i=1:N]{
    name = sprintf("%05d",i);
#    name2 = sprintf("%05d",i-1);
    
    # Se cocatenan cadenas con .
    set output name.".png";
     
    # poner la etiqueta con t=(i-1)*dt
#    set label 1 't=%3.5g',i*dt at x0-0.1*x0,yf-0.1*yf;

#   if (i % 2 == 1) {
#       plot name.'dat' using 1:2:5 with image,\
#            name'out.dat' using 1:2:(1) with 
#   }
    plot name.'.dat' using 1:2:5 with image
}
