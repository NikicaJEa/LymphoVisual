#!/usr/bin/env bash

# program made by Nikola Vinko
#contact : unispieler@gmail.com
# _______________________________________________________________________
# Clonal Rearrangement Visualisation Tool helps you to visualise IGH Fr1, IGH Fr2, IGK and TRG Clonal Rearrangements
#Copyright (C) GPL v3 2019 Nikola Vinko
#
#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



echo 'Gnomovision version 69, Copyright (C) year name of author
      Gnomovision comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome
      to redistribute it under certain conditions. See: "Copyright.txt" for details.'
echo

sleep 5s


ANALYSIS_DIR=$(pwd)
SCRIPT_BASE="$( cd "$(dirname "$0")" ; pwd -P )"
SCRIPT_NAME=$(basename $0)


if [ "$ANALYSIS_DIR" = "$SCRIPT_BASE" ]; then
	echo "do not run in script folder."
	echo "this MUST be started in the folder containing \"*_output_v1.3.3_automated\" folders!"
	exit 1
fi


source $SCRIPT_BASE/venv/bin/activate


python $SCRIPT_BASE/LymphoVisual_main.py

python $SCRIPT_BASE/LymphoVisual_comparison.py

python $SCRIPT_BASE/Scripts/IGH_FR1_table.py
python $SCRIPT_BASE/Scripts/IGH_FR2_table.py
python $SCRIPT_BASE/Scripts/IGK_table.py
python $SCRIPT_BASE/Scripts/TCRG_table.py

python $SCRIPT_BASE/Scripts/GNUPLOT/IGH_FR1/gnuplot_tsv_generate.py
python $SCRIPT_BASE/Scripts/GNUPLOT/IGH_FR2/gnuplot_tsv_generate.py
python $SCRIPT_BASE/Scripts/GNUPLOT/IGK/gnuplot_tsv_generate.py
python $SCRIPT_BASE/Scripts/GNUPLOT/TCRG/gnuplot_tsv_generate.py



deactivate


