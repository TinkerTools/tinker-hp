set encoding iso_8859_1

set terminal pngcairo size 2400,1500 enhanced font 'Helvetica,30'
set output 'Time_distance.png'

set key spacing 2.0
set key on at 700,15 box
set key width 5
set xrange[0:1500]
set yrange[0:16]
set ylabel "End-to-end distance (\305)" offset -1.5, 0
set xlabel'Time (ps)' offset 0, -0.6
set grid
set xlabel font "Helvetica,30"
set ylabel font "Helvetica,30"
set xtics font "Helvetica,30"
set ytics font "Helvetica,30"

set style line 1 lt 1 lc rgb "red" 
set style line 2 lt 1 lc rgb "green" 

xCOM = -7.879
yCOM = 5.403
zCOM = 1.177
plot 'SMD_output.dat' u 2:(sqrt(($3-xCOM)**2+($4-yCOM)**2+($5-zCOM)**2)) w l ls 1 linewidth 5.0 notitle
set output


