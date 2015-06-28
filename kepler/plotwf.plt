reset
set encoding utf8
set terminal postscript eps color solid enhanced
set output 'keplerwf.eps'

set view map
set pm3d
unset surface
unset key

set size ratio -1

splot 'kepler-wf.dat' u 1:2:3 index 0

system('epstopdf keplerwf.eps && rm *.eps')
