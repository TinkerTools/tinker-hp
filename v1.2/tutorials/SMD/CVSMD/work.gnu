set encoding iso_8859_1

set terminal pngcairo size 2400,1500 enhanced font 'Helvetica,30'
set output 'work.png'

set key spacing 2.0
set key on at 700,15 box
set key width 5
set xrange[0:15]
set yrange[-5:40]
set ylabel "Pulling work (kcal/mol)" offset -1.5, 0
set xlabel "Distance (\305)" offset 0, -0.6
set grid
set xlabel font "Helvetica,30"
set ylabel font "Helvetica,30"
set xtics font "Helvetica,30"
set ytics font "Helvetica,30"

set style line 1 lt 1 lc rgb "red" 
set style line 2 lt 1 lc rgb "green" 

plot "work.dat" w l ls 1 linewidth 5.0 notitle
set output


