# program made by Nikola Vinko
# contact : unispieler@gmail.com
# _______________________________________________________________________#

# LymphoVisual helps you to visualise IGH Fr1, IGH Fr2, IGK and TRG Clonal Rearrangements
# Copyright (C) GPL v3  

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


import matplotlib
matplotlib.use('Agg')
import os
import sys
import plotly
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from plotly.graph_objs import Bar, Figure, Layout
import pandas
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import time
import re




def IGH_FR1(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
            region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
            srr3, srr4, ):
    print("Running IGH_FR1 for patient: " + patient_id)

	# creating a patients directory in ./plots directory
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)  
    
	# creating ./CSV_files/ directory
    if not os.path.exists(loc_save):  
        os.makedirs(loc_save)
    
	# selecting the tsv file
    tsv_file = loc + region + region_graph  
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)

    print("Running IGH_FR1 for file: " + tsv_file)
	
	# converting tsv file into csv file
    try:
        csv_table = pandas.read_table(tsv_file, sep="\t")  
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")
	
	# saving the csv file into ./CSV_files/
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv", index=False)  
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

	# VJ Usage percent Plot
    print("Plot VJ Usage percent Plot")

    tsv_file = tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)

    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)
	
	# opening the tsv file
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
		# selecting the V_segs
        V_segs = data1.iloc[:, svc2]  

    z1 = []
    for row in range(rows):
        z2 = []
		# selects the tsv file columns by index
        for col in range(srr1, srr2):
			# selects the values from the column
            val = data1.iloc[row, col]  
            z2.append(val)
        z1.append(z2)
	
	# creating an array that will be used for plotting
    data = np.array(z1)
	
	# selecting  J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]  
    
        x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)
	
	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)

    # setting the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('% Total Reads')
	
	# setting up the x axis
    ticksx = np.arange(0.5, columns, 1)  
    plt.xticks(ticksx, column_names)
	
	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)  
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Percent Graph', fontsize=16)
	
	# saving the VJ Usage Percent Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR1 VJ Usage Percent Graph.png", dpi=200)
    plt.close(fig)

    
	# VJ Usage RAW count graph
    print("Plot VJ Usage RAW Count Graph")

    # selecting tsv file
    tsv_file = tsv_file = loc + region + region_graph2  
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
	
	# saving csv file
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)  
	
	# selecting V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2]  

    z1 = []
	# selects the column by index
    for row in range(svc2):
        z2 = []
        for col in range(srr3, srr4):
			# selects the values from cells  of the current column
            val = data1.iloc[row, col]  
            z2.append(val)
        z1.append(z2)  
	
	# creating an array that will be used for plotting
    data = np.array(z1)
	
	# selecting J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]  
    
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)
	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)

    # setting up the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('Raw Count')

    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Raw Graph', fontsize=16)

	# saving the VJ Usage Raw Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR1 VJ Usage Raw Graph.png", dpi=200)

    
	# Pie Chart top 10 sequences from top 500 Reads
    print("Plot Pie Chart top 10 sequences from top 500 Reads")
	
	# selecting tsv file
    tsv_file = loc + region + "/index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.tsv"
    csv_table = pd.read_table(tsv_file, sep='\t')
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv", index=False)
	
	# selecting V and J segs
    with open(loc_save + "index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
        V_segs = fastq_data.iloc[1:, 7]
        J_segs = fastq_data.iloc[1:, 8]
        total_reads = fastq_data.iloc[1:, 9]

    dic = {}
    x=[]
    y=[]

    for i in range(1, len(V_segs) + 1):
        name = V_segs[i] + "-" + J_segs[i]
        x.append(name)
        y.append(total_reads[i])

    colors = ['gold', 'lawngreen', 'aqua', 'orange', "black", "magenta"
        , "dodgerblue", "maroon", "violet", "pink"]

    # plotting
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.set_position([0.5, 0.1, 0.5, 0.8])
    patches, texts = plt.pie(y, colors=colors, shadow=True, )
    labels = ['{0} - {1:.2f} %'.format(i, float(j)) for i, j in zip(x, y)]
    sort_legend = True
    plt.title("Top 10 Merged Sequences")
    plt.legend(patches, labels, loc='best', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)

    plt.axis('equal')
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR1 Pie Chart.png", dpi=300)
    plt.close(fig)


def IGH_FR2(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
            region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
            srr3, srr4):
    print("Running IGH_FR2 for patient: " + patient_id)
	
	# creating a patients directory in ./plots directory
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
	
	# creating ./CSV_files/ directory
    if not os.path.exists(loc_save):
        os.makedirs(loc_save)
	
	# selecting the tsv file
    tsv_file = loc + region + region_graph
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)
    print("Running IGH_FR2 for file: " + tsv_file)
	
	# converting tsv file into csv file
    try:
        csv_table = pd.read_table(tsv_file, sep='\t')
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("IGH_FR2 Data read done ")

	# saving the csv file into ./CSV_files/
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv")
        csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv", index=False)

    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    # VJ Usage percent Plot
    print("Plot VJ Usage percent plot")

    tsv_file = tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)

	# selecting the V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2]
    
    z1 = []
    for row in range(rows):
        z2 = []
		# selects the tsv file columns by index
        for col in range(srr1, srr2):
			# selects the values from the column
            val = data1.iloc[row, col]  
            z2.append(val)
        z1.append(z2)

	# creating an array that will be used for plotting
    data = np.array(z1)

	# selecting  J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)

	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)  

	# setting the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('% Total Reads')

	# setting up the x axis
    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Percent Graph', fontsize=16)

	# saving the VJ Usage Percent Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR2 VJ Usage Percent Graph.png",
                dpi=200)  
    plt.close(fig)

    # VJ Usage RAW count graph
    print("Plot VJ Usage RAW count graph")
	
	# selecting tsv file
    tsv_file = tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)

	# selecting V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2]

    z1 = []
	# selects the column by index
    for row in range(svc2):
        z2 = []
        for col in range(srr3, srr4): 
			# selects the values from cells  of the current column
            val = data1.iloc[row, col]  
            z2.append(val)
        z1.append(z2)
		
	# creating an array that will be used for plotting
    data = np.array(z1)

	# selecting J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)

	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)

    # setting up the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('Raw Count')

    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Raw Graph', fontsize=16)

	# saving the VJ Usage Raw Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR2 VJ Usage Raw Graph.png", dpi=200)

    # Pie Chart top 10 sequences from top 500 Reads
    print("Plot Pie Chart top 10 from top 500 Reads")

	# selecting tsv file
    tsv_file = loc + region + "/index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.tsv"
    csv_table = pd.read_table(tsv_file, sep='\t')
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv", index=False)

	# selecting V and J segs
    with open(loc_save + "index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
        V_segs = fastq_data.iloc[1:, 7]
        J_segs = fastq_data.iloc[1:, 8]
        total_reads = fastq_data.iloc[1:, 9]

    x=[]
    y=[]

    for i in range(1, len(V_segs) + 1):
        name = V_segs[i] + "-" + J_segs[i]
        x.append(name)
        y.append(total_reads[i])

    colors = ['gold', 'lawngreen', 'aqua', 'orange', "black", "magenta"
        , "dodgerblue", "maroon", "violet", "pink"]

	# plotting
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.set_position([0.5, 0.1, 0.5, 0.8])
    patches, texts = plt.pie(y, colors=colors, shadow=True, )
    labels = ['{0} - {1:.2f} %'.format(i, float(j)) for i, j in zip(x, y)]
    sort_legend = True
    plt.title("Top 10 Merged Sequences")
    plt.legend(patches, labels, loc='best', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)

    plt.axis('equal')
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGH_FR2 Pie Chart.png", dpi=300)
    plt.close(fig)


def IGK(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
        region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
        srr3, srr4):
    print("Running IGK for patient: " + patient_id)
	
	# creating a patients directory in ./plots directory
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
	
	# creating ./CSV_files/ directory
    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

	
	# selecting the tsv file
    tsv_file = loc + region + region_graph
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)

    print("Running IGK for file: " + tsv_file)

	# converting tsv file into csv file
    try:
        csv_table = pd.read_table(tsv_file, sep='\t')

    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")
	
	# saving the csv file into ./CSV_files/
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv")
        csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    # VJ Usage percent Plot
    print("Plot VJ Usage percent Plot")

    tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)

	# opening the tsv file
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
		# selecting the V_segs
        V_segs = data1.iloc[:, svc2]

    z1 = []
    for row in range(rows):
        z2 = []
		# selects the tsv file columns by index
        for col in range(srr1, srr2):
			# selects the values from the column
            val = data1.iloc[row, col]  
            z2.append(val)
        z1.append(z2)

	# creating an array that will be used for plotting
    data = np.array(z1)

	# selecting  J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)	
	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', "g"]
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))  
    hopa = np.zeros(columns)  

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):  
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)  

    # setting the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('% Total Reads')  

	# setting up the x axis
    ticksx = np.arange(0.5, columns, 1)  
    plt.xticks(ticksx, column_names)

	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)  
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Percent Graph', fontsize=16)

	# saving the VJ Usage Percent Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGK VJ Usage Percent Graph.png", dpi=200)

    # VJ Usage RAW count graph
    print("Plot RAW")

	# selecting tsv file
    tsv_file = loc + region + region_graph2  
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage_family.tsv", index=False)

	# selecting V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2] 
    z1 = []
    for row in range(rows):
        z2 = []
		# selects the column by index
        for col in range(srr3, srr4):  
		# selects the values from cells  of the current column
            val = data1.iloc[row, col]
            z2.append(val)
        z1.append(z2)

	# creating an array that will be used for plotting
    data = np.array(z1)

	# selecting J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage_family.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)

	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', "g"]
	# Convert positions to 1D array
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))  
    hopa = np.zeros(columns)  

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows): 
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)  

    # setting up the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('% Raw Count') 

        # setting up the x axis
    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Raw Graph', fontsize=16)

	# saving the VJ Usage Raw Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGK VJ Usage Raw Graph.png", dpi=200)  
    plt.close(fig)

    # Pie Chart top 10 from top 500 Reads
    print("Plot Pie Chart top 10 from top 500 Reads")

	# selecting tsv file
    tsv_file = loc + region + "/index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.tsv"
    csv_table = pd.read_table(tsv_file, sep='\t')
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv", index=False)

	# selecting V and J segs
    with open(loc_save + "index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
        V_segs = fastq_data.iloc[1:, 7]
        J_segs = fastq_data.iloc[1:, 8]
        total_reads = fastq_data.iloc[1:, 9]

    x=[]
    y=[]

    for i in range(1, len(V_segs) + 1):
        name = V_segs[i] + "-" + J_segs[i]
        x.append(name)
        y.append(total_reads[i])

    colors = ['gold', 'lawngreen', 'aqua', 'orange', "black", "magenta"
        , "dodgerblue", "maroon", "violet", "pink"]

	# plotting
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.set_position([0.55, 0.1, 0.5, 0.8])
    patches, texts = plt.pie(y, colors=colors, shadow=True, )
    labels = ['{0} - {1:.2f} %'.format(i, float(j)) for i, j in zip(x, y)]
    sort_legend = True
    plt.title("Top 10 Merged Sequences")
    plt.legend(patches, labels, loc='best', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)

    plt.axis('equal')
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_IGK Pie Chart.png", dpi=300)
    plt.close(fig)


def TCRG(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
         region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
         srr3, srr4):
    print("Running TCRG for patient: " + patient_id)

	# creating a patients directory in ./plots directory
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

	# creating ./CSV_files/ directory
    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

	# selecting the tsv file
    tsv_file = loc + region + region_graph
    if not os.path.exists(tsv_file):
        print("ERROR: file not found: " + tsv_file)

    print("Running TCRG for file: " + tsv_file)

	# converting tsv file into csv file
    try:
        csv_table = pd.read_table(tsv_file, sep='\t')
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

	# saving the csv file into ./CSV_files/
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id + ".fastq_read_summary.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    
    # VJ Usage percent Plot
    print("Plot VJ Usage percent Plot")

	# opening the tsv file
    tsv_file = tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage.tsv", index=False)

	# selecting V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage.tsv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2]

    z1 = []
    for row in range(rows):
        z2 = []
		# selects the tsv file columns by index
        for col in range(srr1, srr2):
			# selects the values from the column
            val = data1.iloc[row, col]
            z2.append(val)
        z1.append(z2)
		
	# creating an array that will be used for plotting
    data = np.array(z1)

	# selecting  J_segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage.tsv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)

	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'g']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns) 

    # setting the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('% Total Reads')

	# setting up the x axis
    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Percent Graph', fontsize=16)

	# saving the VJ Usage Percent Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_TCRG VJ Usage Percent Graph.png", dpi=200)
    plt.close(fig)

    # VJ Usage RAW count graph
    print("VJ Usage RAW count graph plot")

	# selecting tsv file
    tsv_file = tsv_file = loc + region + region_graph2
    csv_table = pd.read_table(tsv_file, sep='\t', skiprows=2)
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_VJ_usage.csv", index=False)

	# selecting V segs
    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage.csv") as f:
        data1 = pd.read_csv(f)
        V_segs = data1.iloc[:, svc2]

    z1 = []
    for row in range(rows):
        z2 = []
		# selects the column by index
        for col in range(srr3, srr4):
            val = data1.iloc[row, col]
			# selects the values from cells  of the current column
            z2.append(val)
        z1.append(z2)
	# creating an array that will be used for plotting
    data = np.array(z1)

    with open(loc_save + "index" + patient_id + ".fastq_VJ_usage.csv") as f:
        data2 = pd.read_csv(f, header=None, nrows=1)
        J_segs = data2.iloc[0, sjc2:]
    x1 = []
    for key in J_segs:
        x1.append(key.replace(".1", ""))

    y1 = []
    for value in V_segs:
        y1.append(value)

    column_names = x1
    row_names = y1

    fig = plt.figure()
    ax = Axes3D(fig)
	# Work out matrix dimensions
    lx = len(data[0])  
    ly = len(data[:, 0])
	# Set up a mesh of positions
    xpos = np.arange(0, lx, 1)  
    ypos = np.arange(0, ly, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'g']
	# Convert positions to 1D array
    xpos = xpos.flatten()  
    ypos = ypos.flatten()
    zpos = np.zeros(lx * ly)

    x_data, y_data = np.meshgrid(np.arange(columns), np.arange(rows))
    hopa = np.zeros(columns)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = data.flatten()
    for i in range(rows):
        ax.bar3d(x_data[i], y_data[i], hopa, 1, 1, data[i], alpha=0.5, color=colors[i] * columns)  
	
	# setting up the labels
    ax.w_xaxis.set_ticklabels(column_names)
    ax.w_yaxis.set_ticklabels(row_names)
    ax.set_xlabel('J_Segments')
    ax.set_ylabel('V_Segments')
    ax.set_zlabel('Raw Count')

	# setting up the x axis
    ticksx = np.arange(0.5, columns, 1)
    plt.xticks(ticksx, column_names)

	# setting up the y axis
    ticksy = np.arange(0.5, rows, 1)
    plt.yticks(ticksy, row_names)
    plt.suptitle('VJ Usage Raw Graph', fontsize=16)

	# saving the VJ Usage Raw Graph
    ax.view_init(elev=40., azim=-62)
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_TCRG VJ Usage Raw Graph.png", dpi=200)  
    plt.close(fig)

    # Pie Chart top 10 from top 500 Reads
    print("Plot Pie Chart top 10 from top 500 Reads")

	# selecting tsv file
    tsv_file = loc + region + "/index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.tsv"
    csv_table = pd.read_table(tsv_file, sep='\t')
    csv_table.to_csv(loc_save + 'index' + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv", index=False)

	# selecting V and J segs
    with open(loc_save + "index" + patient_id + ".fastq_read_summary_merged_top10_searchtop500.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
        V_segs = fastq_data.iloc[1:, 7]
        J_segs = fastq_data.iloc[1:, 8]
        total_reads = fastq_data.iloc[1:, 9]

    
    x=[]
    y=[]

    for i in range(1, len(V_segs) + 1):
        name = V_segs[i] + "-" + J_segs[i]
        x.append(name)
        y.append(total_reads[i])
    
    colors = ['gold', 'lawngreen', 'aqua', 'orange', "black", "magenta"
        , "dodgerblue", "maroon", "violet", "pink"]

	# plotting
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.set_position([0.5, 0.1, 0.5, 0.8])
    patches, texts = plt.pie(y, colors=colors, shadow=True, )
    labels = ['{0} - {1:.2f} %'.format(i, float(j)) for i, j in zip(x, y)]
    sort_legend = True
    plt.title("Top 10 Merged Sequences")
    plt.legend(patches, labels, loc='best', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)

    plt.axis('equal')
    plt.savefig("./plots/" + patient_id[0:3] + "_LymphoVisual/" + patient_id + "_TCRG Pie Chart.png", dpi=300)
    plt.close(fig)


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
		
		# select a patient
        for patient_id in lista_id:
            for folder in onlyfiles:
                if folder.startswith("index" + patient_id):
                    print("Working on: " + folder)

                    # PARAMETERS for IGH_FR1

					# location of the tsv files
                    loc = "."  
					# location where to save the csv files
                    loc_save = "./CSV_files_IGH_FR1/"  
                    region = "/IGH_FR1_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary_family.tsv"
					# dictionary with all possible combinations of V-J segs
                    dic = {"VH1-J1": [], "VH1-J2": [], "VH1-J3": [], "VH1-J4": [], "VH1-J5": [], "VH1-J6": [], "VH1-none": [],
                           "VH2-J1": [], "VH2-J2": [], "VH2-J3": [], "VH2-J4": [], "VH2-J5": [], "VH2-J6": [], "V2-none":  [],
                           "VH3-J1": [], "VH3-J2": [], "VH3-J3": [], "VH3-J4": [], "VH3-J5": [], "VH3-J6": [], "VH3-none": [],
                           "VH4-J1": [], "VH4-J2": [], "VH4-J3": [], "VH4-J4": [], "VH4-J5": [], "VH4-J6": [], "VH4-none": [],
                           "VH5-J1": [], "VH5-J2": [], "VH5-J3": [], "VH5-J4": [], "VH5-J5": [], "VH5-J6": [], "VH5-none": [],
                           "VH6-J1": [], "VH6-J2": [], "VH6-J3": [], "VH6-J4": [], "VH6-J5": [], "VH6-J6": [], "VH6-none": [],
                           "VH7-J1": [], "VH7-J2": [], "VH7-J3": [], "VH7-J4": [], "VH7-J5": [], "VH7-J6": [], "VH7-none": [],

                           }
                    # IGH_FR1 parameters for selecting the V,J and total_reads values
                    svr = 1					
                    svc = 7  
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # IGH_FR1 VJ Usage Graph parametars
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage_family.tsv"
                    svc2 = 8
                    sjc2 = 9
                    rows = 8
                    columns = 7
					# parameters for percentage graph
                    srr1 = 9  
                    srr2 = 16
					# parameters for raw count graph
                    srr3 = 1 
                    srr4 = 8

                    try:
                        IGH_FR1(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                                region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                                srr3, srr4)

                    except Exception:
                        pass

    mypath = "./IGH_FR2_output_v1.3.3_automated"
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
                    
					# PARAMETERS for IGH_FR2

					# location of the tsv files
                    loc = "." 
					# location where to save the csv files
                    loc_save = "./CSV_files_IGH_FR2/"  
                    region = "/IGH_FR2_output_v1.3.3_automated"

                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary_family.tsv"
					# dictionary with all possible combinations of V-J segs
                    dic = {"VH1-J1": [], "VH1-J2": [], "VH1-J3": [], "VH1-J4": [], "VH1-J5": [], "VH1-J6": [], "VH1-none": [],
                           "VH2-J1": [], "VH2-J2": [], "VH2-J3": [], "VH2-J4": [], "VH2-J5": [], "VH2-J6": [], "V2-none":  [],
                           "VH3-J1": [], "VH3-J2": [], "VH3-J3": [], "VH3-J4": [], "VH3-J5": [], "VH3-J6": [], "VH3-none": [],
                           "VH4-J1": [], "VH4-J2": [], "VH4-J3": [], "VH4-J4": [], "VH4-J5": [], "VH4-J6": [], "VH4-none": [],
                           "VH5-J1": [], "VH5-J2": [], "VH5-J3": [], "VH5-J4": [], "VH5-J5": [], "VH5-J6": [], "VH5-none": [],
                           "VH6-J1": [], "VH6-J2": [], "VH6-J3": [], "VH6-J4": [], "VH6-J5": [], "VH6-J6": [], "VH6-none": [],
                           "VH7-J1": [], "VH7-J2": [], "VH7-J3": [], "VH7-J4": [], "VH7-J5": [], "VH7-J6": [], "VH7-none": [],

                           }
					# IGH_FR2 parameters for selecting the V,J and total_reads values
                    svr = 1  
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # IGH_FR1 VJ Usage Graph parametars
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage_family.tsv"
                    svc2 = 8
                    sjc2 = 9
                    rows = 8
                    columns = 7
					# parameters for percentage graph
                    srr1 = 9  
                    srr2 = 16
					# parameters for raw count graph
                    srr3 = 1  
                    srr4 = 8

                    try:
                        IGH_FR2(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                                region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                                srr3, srr4)
                    except Exception:
                        pass

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
                    
					# PARAMETERS for IGK !

                    # location of the tsv files
                    loc = "."  
					# location where to save the csv files
                    loc_save = "./CSV_files_IGK/"  
                    region = "/IGK_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary_family.tsv"
					# dictionary with all possible combinations of V-J segs
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
					# IGK parameters for selecting the V,J and total_reads values
                    svr = 1  
                    svc = 7 
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # IGK VJ Usage Graph parametars
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage_family.tsv"
                    svc2 = 8
                    sjc2 = 9
                    rows = 9
                    columns = 7
					# parameters for percentage graph
                    srr1 = 9  
                    srr2 = 16
					# parameters for RAW count graph
                    srr3 = 1  
                    srr4 = 8

                    try:
                        IGK(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                            region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                            srr3, srr4)
                    except Exception:
                        pass

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

                    # PARAMETERS for TCRG
					
					# location of the tsv files
                    loc = "."  
					# location where to save the csv files
                    loc_save = "./CSV_files_TCRG/"  
                    region = "/TCRG_output_v1.3.3_automated"
                    folder_name = "./plots/" + patient_id[0:3] + "_LymphoVisual"
                    region_graph = "/index" + patient_id + ".fastq_read_summary.tsv"
					# dictionary with all possible combinations of V-J segs
                    dic = {"Vg2-Jg1/2": [], "Vg2-JgP": [], "Vg2-JgP1": [], "Vg2-JgP2": [], "Vg2-none": [],
                           "Vg3-Jg1/2": [], "Vg3-JgP": [], "Vg3-JgP1": [], "Vg3-JgP2": [], "Vg3-none": [],
                           "Vg4-Jg1/2": [], "Vg4-JgP": [], "Vg4-JgP1": [], "Vg4-JgP2": [], "Vg4-none": [],
                           "Vg5-Jg1/2": [], "Vg5-JgP": [], "Vg5-JgP1": [], "Vg5-JgP2": [], "Vg5-none": [],
                           "Vg8-Jg1/2": [], "Vg8-JgP": [], "Vg8-JgP1": [], "Vg8-JgP2": [], "Vg8-none": [],
                           "Vg9-Jg1/2": [], "Vg9-JgP": [], "Vg9-JgP1": [], "Vg9-JgP2": [], "Vg9-none": [],
                           "Vg10-Jg1/2": [], "Vg10-JgP": [], "Vg10-JgP1": [], "Vg10-JgP2": [], "Vg10-none": [],
                           "Vg11-Jg1/2": [], "Vg11-JgP": [], "Vg11-JgP1": [], "Vg11-JgP2": [], "Vg11-none": [],
                           "none-Jg1/2": [], "none-JgP": [], "none-JgP1": [], "none-JgP2": [], "none-none": []
                           }
					# TCRG parameters for selecting the V,J and total_reads values
                    svr = 1 
                    svc = 7  
                    sjr = 1
                    sjc = 8
                    srr = 1
                    src = 9

                    # TCRG VJ Usage Graph parametars
                    region_graph2 = "/index" + patient_id + ".fastq_VJ_usage.tsv"
                    svc2 = 7
                    sjc2 = 8
                    rows = 9
                    columns = 5
					# parameters for percentage graph
                    srr1 = 8  # 
                    srr2 = 13
					# parameters for RAW count graph
                    srr3 = 1  
                    srr4 = 6

                    TCRG(folder_name, loc, loc_save, region, patient_id, region_graph, dic, svr, svc, sjr, sjc, srr, src,
                         region_graph2, svc2, sjc2, rows, columns, srr1, srr2,
                         srr3, srr4)

                    break




if __name__ == "__main__":
    main(sys.argv[1:])

