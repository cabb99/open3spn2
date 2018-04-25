#!/usr/bin/gnuplot 

set term postscript enhanced color eps "Times-Roman,24"

set encoding iso_8859_1
set output "all_splines.eps"
set size 3.0,6.0
#set grid
set multiplot

set xlabel "{/Times-Roman-Italic r}  [\305]" font "Times-Roman,36"
set xtics font ",30"
set xrange [1:12]

set ylabel "U(r) {/Times-Roman [kCal/mol]}" font "Times-Roman-Italic,36"
set ytics font ",30"
set yrange [-2:7]

# P-P
set origin 0.0,5.0
set size 1.0,1.0
set title "P-P" font ",36"

plot \
"knots/P-P.param" u 1:2 w points title "Knots" lw 10 pt 7 lc 1, \
"splines/P-P.table" u 2:3 w lines title "U_{corr}" lw 8 lt 3 lc 0, \
"total_potential/P-P.table" u 2:3 w lines title "U_{eff}" lw 8 lt 1 lc 1

# P-Na
set origin 1.0,5.0
set size 1.0,1.0
set title "P-Na" font ",36"

plot \
"knots/P-Na.param" u 1:2 w points notitle lw 10 pt 7 lc 2, \
"splines/P-Na.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/P-Na.table" u 2:3 w lines notitle lw 8 lt 1 lc 2

# P-Mg
set origin 2.0,5.0
set size 1.0,1.0
set title "P-Mg" font ",36"

plot \
"knots/P-Mg.param" u 1:2 w points notitle lw 10 pt 7 lc 3, \
"splines/P-Mg.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/P-Mg.table" u 2:3 w lines notitle lw 8 lt 1 lc 3

# P-Cl
set origin 0.0,4.0
set size 1.0,1.0
set title "P-Cl" font ",36"

plot \
"knots/P-Cl.param" u 1:2 w points notitle lw 10 pt 7 lc 4, \
"splines/P-Cl.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/P-Cl.table" u 2:3 w lines notitle lw 8 lt 1 lc 4

# P-N
set origin 1.0,4.0
set size 1.0,1.0
set title "P-N" font ",36"

plot \
"knots/P-N.param" u 1:2 w points notitle lw 10 pt 7 lc 5, \
"splines/P-N.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/P-N.table" u 2:3 w lines notitle lw 8 lt 1 lc 5


# Na-Na
set origin 2.0,4.0
set size 1.0,1.0
set title "Na-Na" font ",36"

plot \
"knots/Na-Na.param" u 1:2 w points notitle lw 10 pt 7 lc 6, \
"splines/Na-Na.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Na-Na.table" u 2:3 w lines notitle lw 8 lt 1 lc 6

# Na-Mg
set origin 0.0,3.0
set size 1.0,1.0
set title "N-Mg" font ",36"

plot \
"knots/Na-Mg.param" u 1:2 w points notitle lw 10 pt 7 lc 7, \
"splines/Na-Mg.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Na-Mg.table" u 2:3 w lines notitle lw 8 lt 1 lc 7


# Na-Cl
set origin 1.0,3.0
set size 1.0,1.0
set title "Na-Cl" font ",36"

plot \
"knots/Na-Cl.param" u 1:2 w points notitle lw 10 pt 7 lc 8, \
"splines/Na-Cl.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Na-Cl.table" u 2:3 w lines notitle lw 8 lt 1 lc 8

# Na-N
set origin 2.0,3.0
set size 1.0,1.0
set title "Na-N" font ",36"

plot \
"knots/Na-N.param" u 1:2 w points notitle lw 10 pt 7 lc 9, \
"splines/Na-N.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Na-N.table" u 2:3 w lines notitle lw 8 lt 1 lc 9

# Mg-Mg
set origin 0.0,2.0
set size 1.0,1.0
set title "Mg-Mg" font ",36"

plot \
"knots/Mg-Mg.param" u 1:2 w points notitle lw 10 pt 7 lc 10, \
"splines/Mg-Mg.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Mg-Mg.table" u 2:3 w lines notitle lw 8 lt 1 lc 10


# Mg-Cl
set origin 1.0,2.0
set size 1.0,1.0
set title "Mg-Cl" font ",36"

plot \
"knots/Mg-Cl.param" u 1:2 w points notitle lw 10 pt 7 lc 11, \
"splines/Mg-Cl.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Mg-Cl.table" u 2:3 w lines notitle lw 8 lt 1 lc 11

# Mg-N
set origin 2.0,2.0
set size 1.0,1.0
set title "Mg-N" font ",36"

plot \
"knots/Mg-N.param" u 1:2 w points notitle lw 10 pt 7 lc 12, \
"splines/Mg-N.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Mg-N.table" u 2:3 w lines notitle lw 8 lt 1 lc 12


# Cl-Cl
set origin 0.0,1.0
set size 1.0,1.0
set title "Cl-Cl" font ",36"

plot \
"knots/Cl-Cl.param" u 1:2 w points notitle lw 10 pt 7 lc 13, \
"splines/Cl-Cl.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Cl-Cl.table" u 2:3 w lines notitle lw 8 lt 1 lc 13

# Cl-N
set origin 1.0,1.0
set size 1.0,1.0
set title "Cl-N" font ",36"

plot \
"knots/Cl-N.param" u 1:2 w points notitle lw 10 pt 7 lc 14, \
"splines/Cl-N.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/Cl-N.table" u 2:3 w lines notitle lw 8 lt 1 lc 14

# N-N
set origin 2.0,1.0
set size 1.0,1.0
set title "N-N" font ",36"

plot \
"knots/N-N.param" u 1:2 w points notitle lw 10 pt 7 lc 15, \
"splines/N-N.table" u 2:3 w lines notitle lw 8 lt 3 lc 0, \
"total_potential/N-N.table" u 2:3 w lines notitle lw 8 lt 1 lc 15

# Adding bonded parameters

# A1
set origin 0.0,0.0
set size 1.0,1.0
set ylabel "U({/Symbol-Oblique q}) {/Times-Roman [kCal/mol]}" font "Times-Roman-Italic,36"
set xlabel "{/Symbol-Oblique q}" font "Times-Roman-Italic,36"
set title "N1-N2-N3" font ",36"
set xrange[0:180]
plot \
"knots/A1.param" u 1:2 w points notitle lw 10 pt 7 lc 16, \
"splines/A1.table" u 2:3 w lines notitle lw 8 lt 3 lc 0

set ylabel "U(r) {/Times-Roman [kCal/mol]}" font "Times-Roman-Italic,36"
set xlabel "{/Times-Roman-Italic r}  [\305]" font "Times-Roman,36"
set xrange [3:7]
# B1
set origin 1.0,0.0
set size 1.0,1.0
set title "N1-N2" font ",36"

plot \
"knots/B1.param" u 1:2 w points notitle lw 10 pt 7 lc 17, \
"splines/B1.table" u 2:3 w lines notitle lw 8 lt 3 lc 0

# B2
set origin 2.0,0.0
set size 1.0,1.0
set title "N2-N3" font ",36"

plot \
"knots/B2.param" u 1:2 w points notitle lw 10 pt 7 lc 18, \
"splines/B2.table" u 2:3 w lines notitle lw 8 lt 3 lc 0



