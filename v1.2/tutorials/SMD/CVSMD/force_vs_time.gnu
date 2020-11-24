set encoding iso_8859_1

set terminal pngcairo size 2400,1500 enhanced font 'Helvetica,30'
set output 'Force_time.png'

set key spacing 2.0
set key on at 700,15 box
set key width 5
set xrange[0:1500]
set yrange[-10:15]
set ylabel "Pulling force (kcal/mol.\305^2)" offset -1.5, 0
set xlabel'Time (ps)' offset 0, -0.6
set grid
set xlabel font "Helvetica,30"
set ylabel font "Helvetica,30"
set xtics font "Helvetica,30"
set ytics font "Helvetica,30"

set style line 1 lt 1 lc rgb "red" 
set style line 2 lt 1 lc rgb "green" 

f(x) = 0
xdir = -0.628
ydir = 0.717
zdir = -0.302
plot 'SMD_output.dat' u 2:(xdir*$6+ydir*$7+zdir*$8) w l ls 1 linewidth 5.0 notitle, \
     f(x) w l ls 2 linewidth 5.0 notitle
set output


