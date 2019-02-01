# program made by Nikola Vinko
# contact : unispieler@gmail.com
# _______________________________________________________________________#

# Clonal Rearrangement Visualisation Tool helps you to visualise IGH Fr1, IGH Fr2, IGK and TRG Clonal Rearrangements
#Copyright (C) GPL v3  NikicaJEa

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

import matplotlib
import six
matplotlib.use('Agg')
import inspect
import os
import sys
import plotly
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from plotly.graph_objs import Bar, Figure, Layout
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from shutil import copyfile
import timeit
import random
import re


def TCRG(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
         region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
         srr3, srr4):
		 
	# opening the csv file
    with open(loc_save + "index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the V segs, J segs and total_reads values
        V_segs = fastq_data.iloc[svr:svr + 10, svc]  
        J_segs = fastq_data.iloc[sjr:sjr + 10, sjc]
        total_reads = fastq_data.iloc[srr:srr + 10, sjc + 1]
        total_count = fastq_data.iloc[0, 2]
        merge_count = fastq_data.iloc[svr:svr + 10, svc -1]
    total = str(total_count)
    slash = "/"
    lista = [total]
    for i in range(len(list(V_segs)) - 1):
        lista.append(slash)
    total_reads_rounded = [round(float(elem), 2) for elem in total_reads]

    if len(V_segs) <= 0:
        return

    df = pd.DataFrame()
    df['Total Count'] = lista
    df['V_gene'] = list(V_segs)
    df['J-gene'] = list(J_segs)
    df['Merge Count'] = list(merge_count)
    df['% of Total Reads'] = list(total_reads_rounded)

    def TCRG_table(data, col_width=10.0, row_height=0.625, font_size=14,
                   header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                   bbox=[0, 0, 1, 1], header_columns=0,
                   ax=None, **kwargs):
        if ax is None:
            size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            ax.axis('off')

        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)

        for k, cell in six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0] % len(row_colors)])
        fig.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" +patient_id + "_TCRG_Top10_Merged_seq_table.png")
        return ax

    TCRG_table(df, header_columns=0, col_width=5.0)


def main(argv):
    basepath = os.getcwd()

    if not os.path.exists(basepath):
        print("ERROR: Path \"" + basepath + "\" does not exist!")
        exit(1)

    os.chdir(basepath)

    mypath = "./TCRG_output_v1.3.3_automated"
    if os.path.exists(mypath):
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # filter patient ids
        lista_id = []
        for file in onlyfiles:
            match = re.search(r"index([0-9]+_[0-9]+)\.", file, re.IGNORECASE)
            if match:
                lista_id.append(match.group(1))

        lista_id = set(lista_id)

        for patient_id in lista_id:
            for folder in onlyfiles:
                if folder.startswith("index" + patient_id):
                    print("Working on: " + folder)

                    # needs to be done multiple times, since gnuplot also changes dir
                    os.chdir(basepath)

                    # PARAMETERS for TCRG

					# location of the tsv files
                    loc = "."  
					# location where to save the csv files
                    loc_save = "./CSV_files_TCRG/"  
                    region = "/TCRG_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.tsv"
                    dic = {}
                    svr = 1  # selects V segs , rows
                    svc = 7  # sleects V segs, column
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # VJ parametars
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage.tsv"
                    svc2 = 7
                    sjc2 = 8
                    rows = 9
                    columns = 5
                    srr1 = 8  # for percen
                    srr2 = 13
                    srr3 = 1  # for raw
                    srr4 = 6

                    TCRG(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                         region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                         srr3, srr4)

                    break


if __name__ == "__main__":
    main(sys.argv[1:])
