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

import matplotlib
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
from shutil import copyfile, move
import timeit
import random
import re

def IGK(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
        region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
        srr3, srr4):
    print("Running IGK for patient: " + patient_id)

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

    tsv_file = loc + region + region_graph
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)

    print("Running IGK for file: " + tsv_file)

    try:
        csv_table = pd.read_table(tsv_file, sep='\t')
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

    try:
        print("Saving data to: " + loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    with open(loc_save + "index" + patient_id + ".fastq_read_summary_family.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
        V_segs = fastq_data.iloc[svr:, svc]
        J_segs = fastq_data.iloc[sjr:, sjc]
        total_reads = fastq_data.iloc[srr:, src]


    for i in range(1, len(V_segs) + 1):
        name = V_segs[i] + "-" + J_segs[i]  # selects the V and J segemnt names
        dic[name].append(total_reads[i])

    data = []
    max_el = 0
    for key in dic:  # calculating the max number of elements in a list to know how many bars i have to create for a loop
        if len(dic[key]) > max_el:
            b = len(dic[key])
            max_el = b

    for key in dic:
        if len(dic[key]) < max_el:
            for l in range(max_el):
                dic[key].append(0)
                if len(dic[key]) == max_el:
                    break
                else:
                    continue

    num_stacked_items = max_el
    num_elements = 45
    ind = map(hex, range(0, max_el))
    width = 0.5

    timestamp_start = timeit.default_timer()
    lista = []

    if num_stacked_items > 0:

        with open(folder_name + "/IGK_VJ_Sequence_Frequencies.tsv", "w") as f:
            f.write("empty"+"\t" + "\t".join(key for key in dic) + "\n")  # write header line

            for i in range(0, num_stacked_items):
                lista.append([])
                for value in dic:
                    lista[i].append(dic[value][i])
                f.write("A" + str(i) + "\t" + "\t".join(str(x) for x in lista[i]) + "\n")

        copyfile(os.path.dirname(os.path.abspath(inspect.stack()[0][1])) + "/IGK_VJ_Sequence_Frequencies.gnuplot",
                 folder_name + "/IGK_VJ_Sequence_Frequencies.gnuplot")

        os.chdir(folder_name)

        print("Executing gnuplot")
        os.system("gnuplot IGK_VJ_Sequence_Frequencies.gnuplot")

        try:
            move("./IGK_VJ_Sequence_Frequencies.png", "./" + patient_id + "_IGK_VJ_Sequence_Frequencies.png")

            os.remove("./IGK_VJ_Sequence_Frequencies.gnuplot")
            os.remove("./IGK_VJ_Sequence_Frequencies.tsv")

        except Exception:
            print("Unexpected error:", sys.exc_info()[0])
            pass
    else:
        print("No data found - skipping")

    timestamp_stop = timeit.default_timer()
    print("overall", " time: {:.3f} sec ".format(timeit.default_timer() - timestamp_start),
          "time: {:.3f} min".format((timeit.default_timer() - timestamp_start) / 60))


def main(argv):
    basepath = os.getcwd()

    if not os.path.exists(basepath):
        print("ERROR: Path \"" + basepath + "\" does not exist!")
        exit(1)

    os.chdir(basepath)

    mypath = "./IGK_output_v1.3.3_automated"
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

                    # read frequency
                    loc = "."  # location of the tsv files
                    loc_save = "./CSV_files_IGK/"  # location where to save the csv files
                    region = "/IGK_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary_family.tsv"
                    dic = {"V1-J1": [], "V1-J2": [], "V1-J3": [], "V1-J4": [], "V1-J5": [], "V1-IGKDEL": [], "V1-none": [],
                           "V2-J1": [], "V2-J2": [], "V2-J3": [], "V2-J4": [], "V2-J5": [], "V2-IGKDEL": [], "V2-none": [],
                           "V3-J1": [], "V3-J2": [], "V3-J3": [], "V3-J4": [], "V3-J5": [], "V3-IGKDEL": [], "V3-none": [],
                           "V4-J1": [], "V4-J2": [], "V4-J3": [], "V4-J4": [], "V4-J5": [], "V4-IGKDEL": [], "V4-none": [],
                           "V5-J1": [], "V5-J2": [], "V5-J3": [], "V5-J4": [], "V5-J5": [], "V5-IGKDEL": [], "V5-none": [],
                           "V6-J1": [], "V6-J2": [], "V6-J3": [], "V6-J4": [], "V6-J5": [], "V6-IGKDEL": [], "V6-none": [],
                           "V7-J1": [], "V7-J2": [], "V7-J3": [], "V7-J4": [], "V7-J5": [], "V7-IGKDEL": [], "V7-none": [],
                           "IGKINTR-J1": [], "IGKINTR-J2": [], "IGKINTR-J3": [], "IGKINTR-J4": [], "IGKINTR-J5": [], "IGKINTR-IGKDEL": [], "IGKINTR-none": [],
                           "none-J1": [], "none-J2": [], "none-J3": [], "none-J4": [], "none-J5": [], "none-IGKDEL": [], "none-none": [],

                           }
                    svr = 1  # selects V segs , rows
                    svc = 7  # sleects V segs, column
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # VJ
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage_family.tsv"
                    svc2 = 8
                    sjc2 = 9
                    rows = 9
                    columns = 7
                    srr1 = 9  # for percen
                    srr2 = 16
                    srr3 = 1  # for raw
                    srr4 = 8

                    try:
                        IGK(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                            region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                            srr3, srr4)
                    except Exception:
                        pass

                    break


if __name__ == "__main__":
    main(sys.argv[1:])
