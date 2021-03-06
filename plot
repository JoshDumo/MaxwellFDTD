#!/bin/bash
for i in `seq 0 49`;
do
name="Iteration-${i}"
gnuplot << plotting_many_files
unset key
set xlabel "x"
set ylabel "y"
set xrange [0:50]
set yrange [0:50]
set zrange [0:0.4]
set cbrange[0:0.4]
set pm3d interpolate 1,1
set dgrid3d 100,100 qnorm 2
set contour base
#set cntrlabel  format '%8.3g' font ',7' start 5 interval 20
set cntrparam order 8
set cntrparam bspline
#set style data lines
xyz="./$name"
set datafile separator ","
set terminal pngcairo
set output "$name.png"
set title "Propagation of a Gaussian pulse initiated in the middle"
splot xyz with pm3d
plotting_many_files
done
