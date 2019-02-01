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

def IGH_FR1(folder_name,loc,loc_save,region,patient_id,region_graph,dic,svr,svc,sjr,sjc,srr,src,
       region_graph2,svc2,sjc2,rows,columns,srr1,srr2,
       srr3,srr4):
    
    print("Running IGH_FR1 for patient: " + patient_id)
   
 
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        
    if not os.path.exists(loc_save):
            os.makedirs(loc_save)
    
    
    tsv_file=loc+region+region_graph
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)
        
    print("Running IGH_FR1 for file: " + tsv_file)
    
    try:
        csv_table=pd.read_table(tsv_file,sep='\t')
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])
        
    print("Data read done ")
        
    try:
        print("Saving data to: " + loc_save+'index'+patient_id+".fastq_read_summary_family.csv") 
            
        csv_table.to_csv(loc_save+'index'+patient_id+".fastq_read_summary_family.csv",index=False)
    except:
        print("Oops! " + sys.exc_info()[0])
        
    print("Data save done !")
    
    with open(loc_save+"index"+patient_id+".fastq_read_summary_family.csv") as f:
        fastq_data=pd.read_csv(f,header=None)
        V_segs=fastq_data.iloc[svr:,svc]
        J_segs=fastq_data.iloc[sjr:,sjc]
        total_reads=fastq_data.iloc[srr:,src]
    
    
    
    for i in range(1,len(V_segs)+1):
        name=V_segs[i]+"-"+J_segs[i] # selects the V and J segemnt names
        dic[name].append(total_reads[i])

    data=[]
    max_el=0
    for key in dic: # calculating the max number of elements in a list to know how many bars i have to create for a loop
          if len(dic[key]) > max_el:
            b= len(dic[key])
            max_el=b

    for key in dic:
        if len(dic[key]) < max_el:
            for l in range(max_el):
                dic[key].append(0)
                if len(dic[key]) == max_el:
                    break
                else:
                    continue
                    
    num_stacked_items = max_el
    num_elements = 56
    ind = map(hex, range(0,max_el))
    width = 0.5
    
   
    lista=[]
    
    if num_stacked_items > 0:

        with open(folder_name + "/IGH_FR1_VJ_Sequence_Frequencies.tsv", "w") as f:
            f.write("empty" + "\t" + "\t".join(key for key in dic) + "\n")  # write header line

            for i in range(0, num_stacked_items):
                lista.append([])
                for value in dic:
                    lista[i].append(dic[value][i])
                f.write("A" + str(i) + "\t" + "\t".join(str(x) for x in lista[i]) + "\n")

        copyfile(os.path.dirname(os.path.abspath(inspect.stack()[0][1])) + "/IGH_FR1_VJ_Sequence_Frequencies.gnuplot",
                 folder_name + "/IGH_FR1_VJ_Sequence_Frequencies.gnuplot")

        os.chdir(folder_name)

        print("Executing gnuplot")
        os.system("gnuplot IGH_FR1_VJ_Sequence_Frequencies.gnuplot")

        try:
            move("./IGH_FR1_VJ_Sequence_Frequencies.png", "./" + patient_id + "_IGH_FR1_VJ_Sequence_Frequencies.png")

            os.remove("./IGH_FR1_VJ_Sequence_Frequencies.gnuplot")
            os.remove("./IGH_FR1_VJ_Sequence_Frequencies.tsv")

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

    mypath = "./IGH_FR1_output_v1.3.3_automated"
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


                    # PARAMETERS for IGH_FR1
                    #freqeuncy 
                    loc="." #location of the tsv files
                    loc_save="./CSV_files_IGH_FR1/" # location where to save the csv files
                    region="/IGH_FR1_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph="/index"+patient_id+".fastq_read_summary_family.tsv"
                    dic={"VH1-J1":[],"VH1-J2":[],"VH1-J3":[],"VH1-J4":[],"VH1-J5":[],"VH1-J6":[],"VH1-none":[], # dictionary with all possible combinations of V-J segs
                         "VH2-J1":[],"VH2-J2":[],"VH2-J3":[],"VH2-J4":[],"VH2-J5":[],"VH2-J6":[],"V2-none":[],
                         "VH3-J1":[],"VH3-J2":[],"VH3-J3":[],"VH3-J4":[],"VH3-J5":[],"VH3-J6":[],"VH3-none":[],
                         "VH4-J1":[],"VH4-J2":[],"VH4-J3":[],"VH4-J4":[],"VH4-J5":[],"VH4-J6":[],"VH4-none":[],
                         "VH5-J1":[],"VH5-J2":[],"VH5-J3":[],"VH5-J4":[],"VH5-J5":[],"VH5-J6":[],"VH5-none":[],
                         "VH6-J1":[],"VH6-J2":[],"VH6-J3":[],"VH6-J4":[],"VH6-J5":[],"VH6-J6":[],"VH6-none":[],
                         "VH7-J1":[],"VH7-J2":[],"VH7-J3":[],"VH7-J4":[],"VH7-J5":[],"VH7-J6":[],"VH7-none":[],
                         "none-J1":[],"none-J2":[],"none-J3":[],"none-J4":[],"none-J5":[],"none-J6":[],"none-none":[]

                        }
                    svr=1 #selects V segs , rows
                    svc=7 # sleects V segs, column
                    sjr=1
                    sjc=8
                    srr=1
                    src=9

                    #VJ parametars
                    region_graph2="/index"+patient_id+".fastq_VJ_usage_family.tsv"
                    svc2=8
                    sjc2=9
                    rows=8
                    columns=7
                    srr1=9 #for percen
                    srr2=16
                    srr3=1 # for raw
                    srr4=8
		
                    try:
                       IGH_FR1(folder_name,loc,loc_save,region,patient_id,region_graph,dic,svr,svc,sjr,sjc,srr,src,
                         region_graph2,svc2,sjc2,rows,columns,srr1,srr2,srr3,srr4)
                    except Exception:
                        pass
			
                        break

if __name__ == "__main__":
    main(sys.argv[1:]) 
