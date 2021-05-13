#!/usr/bin/gnuplot -persist
reset
set terminal epslatex size 12cm,16cm color colortext standalone header \
   "\\newcommand{\\ft}[0]{\\footnotesize}"
set output 'scan.tex'
# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror # define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
#set grid back ls 12

# line styles
set style line 1 lt 1 lw 3 lc rgb '#66C2A5' # teal
set style line 2 lt 1 lw 3 lc rgb '#FC8D62' # orange
set style line 3 lt 1 lw 3 lc rgb '#8DA0CB' # lilac
set style line 4 lt 1 lw 3 lc rgb '#E78AC3' # magentat
set style line 5 lt 1 lw 3 lc rgb '#A6D854' # lime green
set style line 6 lt 1 lw 3 lc rgb '#FFD92F' # banana
set style line 7 lt 1 lw 3 lc rgb '#E5C494' # tan
set style line 8 lt 1 lw 3 lc rgb '#B3B3B3' # grey
# color definitions
#set style line 1 lc rgb '#8b1a0e' pt 7 ps 0.5 lt 1 dt 1 lw 2 # --- red
#set style line 2 lc rgb '#5e9c36' pt 7 ps 1 lt 1 dt 2 lw 2 # --- green
#set style line 3 lc rgb '#000000' pt 7 ps 0.5 lt 1 dt 1 lw 1 # --- green
#set style line 4 lc rgb '#5e369c' pt 7 ps 0.5 lt 1dt 1 lw 2 # --- green
#set style line 5 lc rgb '#365e9c' pt 7 ps 1 lt 1 lw 2 # --- green
#set style line 6 lc rgb '#1a0e8b' pt 7 ps 1 lt 1 lw 2 # --- green

unset tics
#plot "fanglei.png" binary filetype=png w rgbimage
unset key
# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
#set grid back ls 12
#set xrange [0.6:1.9]
#set yrange [1e7:1e22]

set grid
#set xlabel '$\nu_0$'
#et ylabel '$P_{eq}, \%$'
#set format y "$10^{%T}$"
#set label 1 '$n_{\varphi} = \hskip 8$' at 50,1e-3  tc ls 1 #   rotate by  78.5 center
#set title "Convergence"
#set logscale y
g = 6
# define the functions depending on the current number
#fstr(N) = sprintf("f%d(x) = a%d*x**(%d)", N, N, 1)
#
## The fitting string for a specific file and the related function
#fitstr(N) = sprintf("fit f%d(x) 'errtm%dn54.dat'   u 1:4 every %d via a%d", N, N, g, N)
#
#N = 3
#
#do for [t=1:N] {
#  eval(fstr(t))
#  eval(fitstr(t))
#}

set multiplot layout 2,1
set key r t
#set key width +3 title '$\sigma_0 \approx 1.65\times10^{-4}$' enhanced
set key box lc rgb '#808080' lt 1
#set term x11
#set xtics 12.0,0.1,36.0
set ytics 0.0,10.0,100.0
set grid
set ylabel 'Polarization ($\%$)'
set xlabel 'a$\gamma$'
#et xlabel 'Energy, GeV'
set title 'Equilibrium polarizations with perfect alignment'
plot [:][0:100] \
  "pol.out" u 2:($4 * 1) title "Total  Polarization"  w l ls 1,\
  "pol.out" u 2:($3 * 1) title "S-T  Polarization"  w l ls 2,\
  "mctdep.out" u 2:(84.*($3/(25. +$3)))  title "M-C poln (total)"  w l lt 3

#set output "EIC.ps"
#test

#et term x11

set title "PLOT 36: ZDR,RQ,0.3mrad 3s,coupl.res., standard SLICK (D.P. Barber)"
plot [:][0:100] \
  "pol.out" u 2:($5 * 1) title "x     Polarization"   w l ls 1,\
  "pol.out" u 2:($6 * 1) title "y     Polarization"   w l ls 2,\
  "pol.out" u 2:($7 * 1) title "s     Polarization"   w l ls 3
