# program made by Nikola Vinko
# contact : unispieler@gmail.com
# _______________________________________________________________________#

# LymphoVisual - Clonal Rearrangement Visualisation Tool helps you to visualise IGH Fr1, IGH Fr2, IGK and TRG Clonal Rearrangements
#Copyright (C) GPL v3  NikicaJEa

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.


set terminal pngcairo  transparent enhanced background rgb 'white' font "arial,20" fontscale 1.0 size 1920, 2160 
set output 'TCRG_VJ_Sequence_Frequencies.png'
set title font "Helvetica,37" 
set title 'LymphoVisual TRG Assay - V - J Sequence Frequencies : Top 200 Sequences'
set grid y
set ylabel "% Total Reads"
set style data histograms
set style histogram columnstacked
set xtics rotate
set boxwidth 0.9 #relative
#set style fill solid border -1
set style fill solid 1.0 noborder

#set terminal svg enhanced background rgb 'black'
set style rectangle fillstyle #noborder

# set nokey # disable legend
set key outside top font ",30"

plot 'TCRG_VJ_Sequence_Frequencies.tsv' using 2 ti col, \
	'' using 3 ti col , \
	'' using 4 ti col , \
	'' using 5 ti col , \
	'' using 6 ti col , \
	'' using 7 ti col , \
	'' using 8 ti col , \
	'' using 9 ti col , \
	'' using 10 ti col, \
	'' using 11 ti col, \
	'' using 12 ti col, \
	'' using 13 ti col, \
	'' using 14 ti col, \
	'' using 15 ti col, \
	'' using 16 ti col, \
	'' using 17 ti col, \
	'' using 18 ti col, \
	'' using 19 ti col, \
	'' using 20 ti col, \
	'' using 21 ti col, \
	'' using 22 ti col, \
	'' using 23 ti col, \
	'' using 24 ti col, \
	'' using 25 ti col, \
	'' using 26 ti col, \
	'' using 27 ti col, \
	'' using 28 ti col, \
	'' using 29 ti col, \
	'' using 30 ti col, \
	'' using 31 ti col, \
	'' using 32 ti col, \
	'' using 33 ti col, \
	'' using 34 ti col, \
	'' using 35 ti col, \
	'' using 36 ti col, \
	'' using 37 ti col, \
	'' using 38 ti col, \
	'' using 39 ti col, \
	'' using 40 ti col, \
	'' using 41 ti col, \
	'' using 42 ti col, \
	'' using 43 ti col, \
	'' using 44 ti col, \
	'' using 45 ti col, \
	'' using 46 ti col, \
  



#, '' using 3:xticlabels(1) t "Var 2"
